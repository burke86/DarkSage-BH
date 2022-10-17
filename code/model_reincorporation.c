#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



void reincorporate_gas(int centralgal, double dt, double time, int k_now)
{
    // use of "centralgal" here as a label is from old code, when this wasn't called for satelites.  This can now be called for any galaxy.
    int p = centralgal;
  double reincorporated, metallicity, eject_sum, eject_metals_sum, reincTime_fac, Vhot;
    int k, kk;
    
//    eject_sum = 0.0;
//    for(kk=0; kk<=N_AGE_BINS; kk++) eject_sum += Gal[p].EjectedMass_Reinc[kk];
//    assert((eject_sum <= 1.01 * Gal[p].EjectedMass) && (eject_sum >= 0.99 * Gal[p].EjectedMass));

    
  if(ReincorporationModel==0)
  {
      // SN velocity is 630km/s, and the condition for reincorporation is that the 
      // halo has an escape velocity greater than this, i.e. V_SN/sqrt(2) = 445.48km/s
      double Vcrit = 445.48 * ReIncorporationFactor;
      assert(Gal[p].HotGas >= Gal[p].MetalsHotGas);
      if(Gal[p].Vvir > Vcrit)
        reincorporated = ( Gal[p].Vvir / Vcrit - 1.0 ) * Gal[p].EjectedMass / (Gal[p].Rvir / Gal[p].Vvir) * dt;
      else
        reincorporated = 0.0;
  }
  else if(ReincorporationModel==1)
  {
      double t_reinc = 1.8e4/UnitTime_in_Megayears * (Hubble_h/Gal[p].Mvir);
      reincorporated = Gal[p].EjectedMass / t_reinc * dt;
  }
  else if (ReincorporationModel==2)
  {
      reincorporated = Gal[p].EjectedMass * (Gal[p].HotGasPotential - Gal[p].prevHotGasPotential) / (Gal[p].prevEjectedPotential - Gal[p].prevHotGasPotential) / STEPS;
  }
  else if ((ReincorporationModel==3) || (ReincorporationModel==4))
  {
//      double reincTime_fac = dmin( pow(Gal[p].prevRvir / Gal[p].Rvir, 2.0/STEPS) , 1.0) ;
      Vhot = (2 * Gal[p].Vvir * Gal[p].CoolScaleRadius) / sqrt(Gal[p].R2_hot_av);
      
      if(Gal[p].prevRvir - Gal[p].prevRhot > 0.0)
          reincTime_fac = pow(dmin(1.0, Gal[p].prevVhot/Vhot) * dmin( 1.0, cube(dmax(Gal[p].prevRvir - sqrt(Gal[p].R2_hot_av), 0.0) / (Gal[p].prevRvir - Gal[p].prevRhot)) ), 1.0/STEPS);
      else
          reincTime_fac = 1.0;
      
      if(!((reincTime_fac <= 1.0)*(reincTime_fac >= 0.0))) 
      {
          printf("reincTime_fac = %e\n", reincTime_fac);
          printf("Vhot, Gal[p].prevRvir, Gal[p].prevRhot = %e, %e, %e\n", Vhot, Gal[p].prevRvir, Gal[p].prevRhot);
      }
      assert(reincTime_fac <= 1.0);
      assert(reincTime_fac >= 0.0);
      Gal[p].ReincTime *= reincTime_fac;  // reduce reincorporation time to account for halo growth
      if(Gal[p].ReincTime <= 0.0 || dt >= Gal[p].ReincTime)
      {
          reincorporated = Gal[p].EjectedMass;
          Gal[p].ReincTime = 0.0;
      }
      else
      {
          reincorporated = Gal[p].EjectedMass * dt / Gal[p].ReincTime;
          Gal[p].ReincTime -= dt; // time to reincorporate what's left is reduced by the time that has elapsed
      }

//      if(Gal[p].ReincTime <= 0.0) Gal[p].ReincTime = 0.001 * dt; // set to small value for instant reincorporation in next sub-time-step
  }
  else
  {
      for(k=0; k<k_now; k++) assert(Gal[p].EjectedMass_Reinc[k] <= 0);
      
      // halo expansion should speed things up.  Mixing implicity happens between ejecta here too
      Vhot = (2 * Gal[p].Vvir * Gal[p].CoolScaleRadius) / sqrt(Gal[p].R2_hot_av);      
      double frac_drop = pow(Gal[p].Rvir / Gal[p].prevRvir, 2.0/STEPS) - 1.0;
      
      if (frac_drop >= 1.0)
      {
          for(k=k_now+1; k<=N_AGE_BINS; k++)
          {              
              Gal[p].EjectedMass_Reinc[k-1] += Gal[p].EjectedMass_Reinc[k];
              Gal[p].EjectedMass_Reinc[k] = 0.0;

              Gal[p].MetalsEjectedMass_Reinc[k-1] += Gal[p].MetalsEjectedMass_Reinc[k];
              Gal[p].MetalsEjectedMass_Reinc[k] = 0.0;
          }
      }
      else if (frac_drop > 0.0)
      {
          double drop_mass=0.0, drop_metals=0.0;
          for(k=k_now+1; k<=N_AGE_BINS; k++)
          {
              drop_mass = frac_drop * Gal[p].EjectedMass_Reinc[k];
              drop_metals = frac_drop * Gal[p].MetalsEjectedMass_Reinc[k];
              
              Gal[p].EjectedMass_Reinc[k] -= drop_mass;
              Gal[p].EjectedMass_Reinc[k-1] += drop_mass;

              Gal[p].MetalsEjectedMass_Reinc[k] -= drop_metals;
              Gal[p].MetalsEjectedMass_Reinc[k-1] += drop_metals;
          }
      }
      
      eject_sum=0.0;
      eject_metals_sum=0.0;
      for(k=k_now; k<=N_AGE_BINS; k++)
      {
//          if((Gal[p].EjectedMass_Reinc[k] < MIN_STARFORMATION) || (Gal[p].MetalsEjectedMass_Reinc[k] < MIN_STARFORMATION * BIG_BANG_METALLICITY))
//          {
//              Gal[p].EjectedMass_Reinc[k] = 0.0;
//              Gal[p].MetalsEjectedMass_Reinc[k] = 0.0;
//          }
          eject_sum += Gal[p].EjectedMass_Reinc[k];
          eject_metals_sum += Gal[p].MetalsEjectedMass_Reinc[k];
      }
      if(!((eject_sum <= 1.01 * Gal[p].EjectedMass) && (eject_sum >= 0.99 * Gal[p].EjectedMass))) 
      {
          printf("eject_sum, Gal[p].EjectedMass = %e, %e\n", eject_sum, Gal[p].EjectedMass);
          printf("eject_metals_sum, Gal[p].MetalsEjectedMass = %e, %e\n", eject_sum, Gal[p].EjectedMass);
      }
      assert((eject_sum <= 1.01 * Gal[p].EjectedMass) && (eject_sum >= 0.99 * Gal[p].EjectedMass));
      assert((eject_metals_sum <= 1.01 * Gal[p].MetalsEjectedMass) && (eject_metals_sum >= 0.99 * Gal[p].MetalsEjectedMass));
      Gal[p].EjectedMass = eject_sum;
      Gal[p].MetalsEjectedMass = eject_metals_sum;

      
      reincorporated = Gal[p].EjectedMass_Reinc[k_now] * dt / (time - AgeBinEdge[k_now+1]);
      if(reincorporated < MIN_STARFORMATION) reincorporated = 0.0;
      if(reincorporated > Gal[p].EjectedMass_Reinc[k_now]) reincorporated = Gal[p].EjectedMass_Reinc[k_now];
      if(!(reincorporated >= 0)) printf("reincorporated, Gal[p].EjectedMass_Reinc[k_now] = %e, %e\n", reincorporated, Gal[p].EjectedMass_Reinc[k_now]);
      assert(reincorporated >= 0);
  }

    
  if(reincorporated > 0.0)
  {
    check_ejected(centralgal);
    if(reincorporated < Gal[p].EjectedMass)
    {
        if(ReincorporationModel<5)
        {
            metallicity = get_metallicity(Gal[p].EjectedMass, Gal[p].MetalsEjectedMass);
            assert(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass);
        }
        else
        {
            metallicity = get_metallicity(Gal[p].EjectedMass_Reinc[k_now], Gal[p].MetalsEjectedMass_Reinc[k_now]);
            assert(Gal[p].EjectedMass_Reinc[k_now] >= Gal[p].MetalsEjectedMass_Reinc[k_now]);
        }
        
        Gal[p].EjectedMass -= reincorporated;
        Gal[p].MetalsEjectedMass -= metallicity * reincorporated;
        Gal[p].HotGas += reincorporated;
        Gal[p].MetalsHotGas += metallicity * reincorporated;
        
        if(ReincorporationModel==5)
        {
            Gal[p].EjectedMass_Reinc[k_now] -= reincorporated;
            Gal[p].MetalsEjectedMass_Reinc[k_now] -= metallicity * reincorporated;
        }

    }
    else
    {
        Gal[p].HotGas += Gal[p].EjectedMass;
        Gal[p].MetalsHotGas += Gal[p].MetalsEjectedMass;
        Gal[p].EjectedMass = 0.0;
        Gal[p].MetalsEjectedMass = 0.0;
        
        if(ReincorporationModel==5)
        {
            Gal[p].EjectedMass_Reinc[k_now] = 0.0;
            Gal[p].MetalsEjectedMass_Reinc[k_now] = 0.0;
        }
    }
  }
  assert(Gal[p].HotGas >= Gal[p].MetalsHotGas);

}


void update_reincorporation_time(int p, double new_ejected_mass, double time, int k_now, double metallicity)
{
    if(ReincorporationModel<5)
        Gal[p].ReincTime = (Gal[p].EjectedMass * Gal[p].ReincTime + new_ejected_mass * Gal[p].ReincTimeFresh) / (Gal[p].EjectedMass + new_ejected_mass);
    else
    {
        // find the age bin that corresponds to ReincTimeFresh in the future
        int k=0;
        for(k=k_now; k<N_AGE_BINS; k++)
        {
            if(Gal[p].ReincTimeFresh <= time - AgeBinEdge[k+1]) break;
        }
        Gal[p].EjectedMass_Reinc[k] += new_ejected_mass;
        Gal[p].MetalsEjectedMass_Reinc[k] += metallicity * new_ejected_mass;
        
    }
}
