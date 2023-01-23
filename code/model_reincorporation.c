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
    // CONSIDER IF THIS SHOULD GO TO FOUNTAIN RESERVOIR FIRST
    if(Gal[p].EjectedMass > 0.0)
    {
        // use time-scale calculated in the code (calculated/modified below and when the halo grows)
        if(dt >= Gal[p].ReincTime)
        {
            Gal[p].HotGas += Gal[p].EjectedMass;
            Gal[p].MetalsHotGas += Gal[p].MetalsEjectedMass;
            Gal[p].EjectedMass = 0.0;
            Gal[p].MetalsEjectedMass = 0.0;
            Gal[p].ReincTime = 0.0; // perhaps should be set to the dynamical time as a default instead, but should make no difference with zero mass
        }
        else
        {
            move_factor = dt / Gal[p].ReincTime;
            Gal[p].HotGas += (Gal[p].EjectedMass * move_factor);
            Gal[p].MetalsHotGas += (Gal[p].MetalsFountainGas * move_factor);
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
            outflow_kinetic = Gal[p].OutflowSpecificEnergy - Gal[p].EjectedPotential - 0.5*sqr(Gal[p].Vvir);
            j_hot = 2.0 * Gal[p].Vvir * Gal[p].CoolScaleRadius;
            discriminant = 2.0 * outflow_kinetic_final - sqr(j_hot/Gal[p].Rvir);
            assert(discriminant > 0.0);
            outflow_radial_speed = sqrt(discriminant);
            
            // iterate until the new turnaround radius is found based on the ejected gas's potential
            double R_min, R_max, R_guess, pot_guess;
            int iter;
            R_min = Gal[p].Rvir;
            R_max = 10 * Gal[p].Rvir;
            for(iter=0; iter<200; iter++)
            {
                R_guess = 0.5 * (R_min + R_max);
                pot_guess = NFW_potential(p, R_guess); // approximating an NFW potential here for efficiency
                specific_energy = pot_guess + 0.5*sqr(j_hot/R_guess); // add tangential kinetic energy assuming angular momentum is conserved
                
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
            
            // use time-scale calculated in the code that gets updated whenever mass is added to this reservoir
            if(dt >= Gal[p].OutflowTime)
            {
                Gal[p].EjectedMass += Gal[p].OutflowGas;
                Gal[p].MetalsEjectedMass += Gal[p].MetalsOutflowMass;
                Gal[p].OutflowGas = 0.0;
                Gal[p].MetalsOutflowMass = 0.0;
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
            Gal[centralgal].FountainTime = (Gal[centralgal].FountainGas * Gal[centralgal].FountainTime + eject_gas * 0.1 / sqrt(Hubble_sqr_z(Halo[Gal[centralgal].HaloNr].SnapNum))) / (Gal[centralgal].FountainGas + eject_gas);
            
            if(dt >= Gal[p].OutflowTime)
            {
                Gal[centralgal].FountainMass += Gal[p].OutflowGas;
                Gal[centralgal].MetalsFountainMass += Gal[p].MetalsOutflowMass;
                Gal[p].OutflowGas = 0.0;
                Gal[p].MetalsOutflowMass = 0.0;
                Gal[p].OutflowTime = 0.0;
            }
            else
            {
                Gal[centralgal].FountainMass += eject_gas;
                Gal[centralgal].MetalsFountainMass += (Gal[p].MetalsOutflowGas * move_factor);
                Gal[p].OutflowGas -= eject_gas;
                Gal[p].MetalsOutflowGas -= (Gal[p].MetalsOutflowGas * move_factor);
                Gal[p].OutflowTime -= dt;
            }
            
        }
        
    }
    
}
