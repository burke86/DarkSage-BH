#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"


void reincorporate_gas(int p, int centralgal, double dt)
{
    double move_factor;
    double sat_rad = get_satellite_radius(p, centralgal);
    const double t_dyn = 0.1 / sqrt(Hubble_sqr_z(Halo[Gal[centralgal].HaloNr].SnapNum));
    
    // reincorporate the fountaining gas
    if(Gal[p].FountainGas > 0.0)
    {
        if(dt >= Gal[p].FountainTime)
        {
            Gal[p].HotGas += Gal[p].FountainGas;
            Gal[p].MetalsHotGas += Gal[p].MetalsFountainGas;
            Gal[p].FountainGas = 0.0;
            Gal[p].MetalsFountainGas = 0.0;
            Gal[p].FountainTime = 0.0;
        }
        else
        {
            move_factor = dt / Gal[p].FountainTime;
            Gal[p].HotGas += (Gal[p].FountainGas * move_factor);
            Gal[p].MetalsHotGas += (Gal[p].MetalsFountainGas * move_factor);
            Gal[p].FountainGas -= (Gal[p].FountainGas * move_factor);
            Gal[p].MetalsFountainGas -= (Gal[p].MetalsFountainGas * move_factor);
            Gal[p].FountainTime -= dt;
        }
    }
    
    
    // recincorporate inflowing gas from the ejected reservoir
    if(Gal[p].EjectedMass > 0.0)
    {
        // use time-scale calculated in the code (calculated/modified below and when the halo grows)
        if(dt >= Gal[p].ReincTime)
        {
            Gal[p].FountainTime = (Gal[p].FountainGas * Gal[p].FountainTime + Gal[p].EjectedMass * t_dyn) / (Gal[p].FountainGas + Gal[p].EjectedMass);
            Gal[p].FountainGas += Gal[p].EjectedMass;
            Gal[p].MetalsFountainGas += Gal[p].MetalsEjectedMass;
            Gal[p].EjectedMass = 0.0;
            Gal[p].MetalsEjectedMass = 0.0;
            Gal[p].ReincTime = 0.0; // perhaps should be set to the dynamical time as a default instead, but should make no difference with zero mass
        }
        else
        {
            move_factor = dt / Gal[p].ReincTime;
            Gal[p].FountainTime = (Gal[p].FountainGas * Gal[p].FountainTime + move_factor*Gal[p].EjectedMass * t_dyn) / (Gal[p].FountainGas + move_factor*Gal[p].EjectedMass);
            Gal[p].FountainGas += (Gal[p].EjectedMass * move_factor);
            Gal[p].MetalsFountainGas += (Gal[p].MetalsEjectedMass * move_factor);
            Gal[p].EjectedMass -= (Gal[p].EjectedMass * move_factor);
            Gal[p].MetalsEjectedMass -= (Gal[p].MetalsEjectedMass * move_factor);
            Gal[p].ReincTime -= dt;
        }
    }
    
    
    // put outflowing gas into the ejected reservoir (or the central's fountain reservoir if appropriate)
    if(Gal[p].OutflowGas > 0.0)
    {
        // update ejected reservoir's reinc time OR central's fountain time first, then move mass accordingly from outflow -> ejected (or fountain of central)
        
        move_factor = dt / Gal[p].OutflowTime;
        const double eject_gas = Gal[p].OutflowGas * dmin(1.0, move_factor);
        
        if(p==centralgal || sat_rad > Gal[centralgal].Rvir)
        {
            // want to energy to calculate where the gas that is ejected gets to before it turns around
            // first need to work out its initial speed
            double outflow_kinetic, j_hot, discriminant, outflow_radial_speed;
            const double thermal = 0.5 * sqr(Gal[p].Vvir);
            outflow_kinetic = Gal[p].OutflowSpecificEnergy - Gal[p].EjectedPotential - thermal;
            j_hot = 2.0 * Gal[p].Vvir * Gal[p].CoolScaleRadius;
            discriminant = 2.0 * outflow_kinetic - sqr(j_hot/Gal[p].Rvir);
            assert(discriminant > 0.0);
            outflow_radial_speed = sqrt(discriminant);
            
            // iterate until the new turnaround radius is found based on the ejected gas's potential
            double R_min, R_max, R_guess, pot_guess, specific_energy;
            int iter;
            R_min = Gal[p].Rvir;
            R_max = 10 * Gal[p].Rvir;
            for(iter=0; iter<200; iter++)
            {
                R_guess = 0.5 * (R_min + R_max);
                pot_guess = NFW_potential(p, R_guess); // approximating an NFW potential here for efficiency
                specific_energy = pot_guess + 0.5*sqr(j_hot/R_guess) + thermal; // add thermal and tangential kinetic energy assuming angular momentum is conserved
                
                if(fabs((specific_energy-Gal[p].OutflowSpecificEnergy)/Gal[p].OutflowSpecificEnergy) < 1e-3)
                    break;
                
                if(specific_energy > Gal[p].OutflowSpecificEnergy)
                    R_max = R_guess;
                else
                    R_min = R_guess;
            }

            // assume time to reincorporate the freshly outflowing gas is that to get to the turnaround point and back, assuming its average speed is half its
            double reinc_time_new = 4.0 * (R_guess - Gal[p].Rvir) / outflow_radial_speed;
            Gal[p].ReincTime = (Gal[p].EjectedMass * Gal[p].ReincTime + eject_gas * reinc_time_new) / (Gal[p].EjectedMass + eject_gas);
            Gal[p].EjectedSpecificEnergy = (Gal[p].EjectedMass * Gal[p].EjectedSpecificEnergy + eject_gas * specific_energy) / (Gal[p].EjectedMass + eject_gas);
            
            // use time-scale calculated in the code that gets updated whenever mass is added to this reservoir
            if(dt >= Gal[p].OutflowTime)
            {
                Gal[p].EjectedMass += Gal[p].OutflowGas;
                Gal[p].MetalsEjectedMass += Gal[p].MetalsOutflowGas;
                Gal[p].OutflowGas = 0.0;
                Gal[p].MetalsOutflowGas = 0.0;
                Gal[p].OutflowTime = 0.0;
            }
            else
            {
                Gal[p].EjectedMass += eject_gas;
                Gal[p].MetalsEjectedMass += (Gal[p].MetalsOutflowGas * move_factor);
                Gal[p].OutflowGas -= eject_gas;
                Gal[p].MetalsOutflowGas -= (Gal[p].MetalsOutflowGas * move_factor);
                Gal[p].OutflowTime -= dt;
            }
        }
        else
        {
            // satellite internal to Rvir, so outflowing gas goes to central's fountain reservoir instead
            Gal[centralgal].FountainTime = (Gal[centralgal].FountainGas * Gal[centralgal].FountainTime + eject_gas * t_dyn) / (Gal[centralgal].FountainGas + eject_gas);
            
            if(dt >= Gal[p].OutflowTime)
            {
                Gal[centralgal].FountainGas += Gal[p].OutflowGas;
                Gal[centralgal].MetalsFountainGas += Gal[p].MetalsOutflowGas;
                Gal[p].OutflowGas = 0.0;
                Gal[p].MetalsOutflowGas = 0.0;
                Gal[p].OutflowTime = 0.0;
            }
            else
            {
                Gal[centralgal].FountainGas += eject_gas;
                Gal[centralgal].MetalsFountainGas += (Gal[p].MetalsOutflowGas * move_factor);
                Gal[p].OutflowGas -= eject_gas;
                Gal[p].MetalsOutflowGas -= (Gal[p].MetalsOutflowGas * move_factor);
                Gal[p].OutflowTime -= dt;
            }
            
        }
        
    }
    
}


void update_galactic_fountain_from_growth(int p)
{ // p should be a central galaxy
    
    // need to call this to make sure that relevant potential energies of reservoirs are updated
    update_disc_radii(p);
    
    double reinc_mass, metallicity;
    const double t_dyn = Gal[p].Rvir / Gal[p].Vvir;
    
    // quantity same for hot, fountain, and outflow reservoirs
    const double thermal = 0.5 * sqr(Gal[p].Vvir);
    const double potential_and_thermal = Gal[p].HotGasPotential + thermal;
    const double j_hot = 2.0 * Gal[p].Vvir * Gal[p].CoolScaleRadius;
    const double hot_specific_energy = potential_and_thermal + 0.5*sqr(j_hot)/Gal[p].R2_hot_av;
    
    if(Gal[p].OutflowGas > 0.0)
    {
        // what new average kinetic energy would be needed to maintain the outflow time-scale?
        const double outflow_distance = Gal[p].Rvir - sqrt(Gal[p].R2_hot_av);
        const double radial_speed_av = outflow_distance / Gal[p].OutflowTime;
        const double radial_speed_init = radial_speed_av + (Gal[p].EjectedPotential + 0.5*sqr(j_hot/Gal[p].Rvir) - Gal[p].HotGasPotential - 0.5*sqr(j_hot)/Gal[p].R2_hot_av) / (2.0 * radial_speed_av);
        
        if(radial_speed_init < 0.0)
        {
            // can't preserve the time-scale for the outflow, assume all the outflowing gas is now fountaining
            Gal[p].FountainTime = (Gal[p].FountainGas * Gal[p].FountainTime + Gal[p].OutflowGas * t_dyn) / (Gal[p].FountainGas + Gal[p].OutflowGas);
            Gal[p].FountainGas += Gal[p].OutflowGas;
            Gal[p].MetalsFountainGas += Gal[p].MetalsOutflowGas;
            Gal[p].OutflowGas = 0.0;
            Gal[p].MetalsOutflowGas = 0.0;
            Gal[p].OutflowTime = 0.0;
        }
        else
        {
            // use energy conservation, assuming outflow timescale remains, to calculate previously outflowing gas that no longer has enough energy to make it outside Rvir
            const double target_outflow_energy = potential_and_thermal + 0.5 * sqr(radial_speed_init);
            reinc_mass = Gal[p].OutflowGas * (target_outflow_energy - Gal[p].OutflowSpecificEnergy) / (target_outflow_energy - hot_specific_energy);
            assert(reinc_mass <= Gal[p].OutflowGas);
            metallicity = get_metallicity(Gal[p].OutflowGas, Gal[p].MetalsOutflowGas);
            Gal[p].FountainTime = (Gal[p].FountainGas * Gal[p].FountainTime + reinc_mass * t_dyn) / (Gal[p].FountainGas + reinc_mass);
            Gal[p].FountainGas += reinc_mass;
            Gal[p].MetalsFountainGas += (metallicity * reinc_mass);
            Gal[p].OutflowGas -= reinc_mass;
            Gal[p].MetalsOutflowGas -= (metallicity * reinc_mass);
            Gal[p].OutflowSpecificEnergy = target_outflow_energy; // this increases because it is the lower-energy material has been lost to the fountain reservoir
        }
    }
    
    if(Gal[p].EjectedMass > 0.0)
    {
        // what new average kinetic energy would be needed to maintain the reincorporation time-scale of the ejected reservoir?
        if(Gal[p].EjectedSpecificEnergy < Gal[p].EjectedPotential + thermal + 0.5*sqr(j_hot/Gal[p].Rvir))
        { // ejected specific energy from last snapshot is no longer enough for the gas to be ejected at all
            Gal[p].FountainTime = (Gal[p].FountainGas * Gal[p].FountainTime + Gal[p].EjectedMass * t_dyn) / (Gal[p].FountainGas + Gal[p].EjectedMass);
            Gal[p].FountainGas += Gal[p].EjectedMass;
            Gal[p].MetalsFountainGas += Gal[p].MetalsEjectedMass;
            Gal[p].EjectedMass = 0.0;
            Gal[p].MetalsEjectedMass = 0.0;
            Gal[p].ReincTime = 0.0;
        }
        else
        {
            double R_min, R_max, R_guess, energy_guess, Rratio;
            int iter;
            const int iter_max=200;
            R_min = Gal[p].Rvir;
            R_max = 10 * Gal[p].Rvir;
            for(iter=0; iter<iter_max; iter++)
            {
                R_guess = 0.5 * (R_min + R_max);
                energy_guess = NFW_potential(p, R_guess) + 0.5*sqr(j_hot/R_guess) + thermal;
                if(fabs((energy_guess-Gal[p].EjectedSpecificEnergy)/Gal[p].EjectedSpecificEnergy) < 1e-3)
                    break;
            }
            assert(iter<iter_max);
            Rratio = (R_guess-Gal[p].Rvir) / (R_guess-Gal[p].prevRvir);
            assert(Rratio < 1.0);
            Gal[p].ReincTime *= Rratio; // reduce reincorporation time based on the reduction in distance
        }

    }
    
    
    
    
}
