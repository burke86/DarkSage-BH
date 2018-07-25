#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



void init_galaxy(int p, int halonr)
{
    int j, step;
    double SpinMag;
    
    if(halonr != Halo[halonr].FirstHaloInFOFgroup)
    {
        printf("Hah?\n");
        ABORT(1);
    }
    
    SpinMag = sqrt(sqr(Halo[halonr].Spin[0]) + sqr(Halo[halonr].Spin[1]) + sqr(Halo[halonr].Spin[2]));
    
    Gal[p].Type = 0;
    
    Gal[p].GalaxyNr = GalaxyCounter;
    GalaxyCounter++;
    
    Gal[p].HaloNr = halonr;
    Gal[p].MostBoundID = Halo[halonr].MostBoundID;
    Gal[p].SnapNum = Halo[halonr].SnapNum - 1;
    
    Gal[p].mergeType = 0;
    Gal[p].mergeIntoID = -1;
    Gal[p].mergeIntoSnapNum = -1;
    Gal[p].dT = -1.0;
    
    for(j = 0; j < 3; j++)
    {
        Gal[p].Pos[j] = Halo[halonr].Pos[j];
        Gal[p].Vel[j] = Halo[halonr].Vel[j];
        if(SpinMag>0)
        {
            Gal[p].SpinStars[j] = Halo[halonr].Spin[j] / SpinMag;
            Gal[p].SpinGas[j] = Halo[halonr].Spin[j] / SpinMag;
            Gal[p].SpinHot[j] = Halo[halonr].Spin[j] / SpinMag;
        }
        else
        {
            Gal[p].SpinStars[j] = 0.0;
            Gal[p].SpinGas[j] = 0.0;
            Gal[p].SpinHot[j] = 0.0;
        }
        Gal[p].SpinClassicalBulge[j] = 0.0;
        Gal[p].SpinSecularBulge[j] = 0.0;
    }
    
    Gal[p].Len = Halo[halonr].Len;
    Gal[p].LenMax = Halo[halonr].Len;
    Gal[p].Vmax = Halo[halonr].Vmax;
    Gal[p].Vvir = get_virial_velocity(halonr, p);
    Gal[p].Mvir = get_virial_mass(halonr, p);
    Gal[p].Rvir = get_virial_radius(halonr, p);
    
    Gal[p].deltaMvir = 0.0;
    
    Gal[p].ColdGas = 0.0;
    Gal[p].StellarMass = 0.0;
    Gal[p].ClassicalBulgeMass = 0.0;
    Gal[p].SecularBulgeMass = 0.0;
    Gal[p].HotGas = 0.0;
    Gal[p].EjectedMass = 0.0;
    Gal[p].BlackHoleMass = 0.0;
    Gal[p].ICS = 0.0;
    
    Gal[p].StarsFromH2 = 0.0;
    Gal[p].StarsInstability = 0.0;
    Gal[p].StarsMergeBurst = 0.0;
    
    Gal[p].MetalsColdGas = 0.0;
    Gal[p].MetalsStellarMass = 0.0;
    Gal[p].ClassicalMetalsBulgeMass = 0.0;
    Gal[p].SecularMetalsBulgeMass = 0.0;
    Gal[p].MetalsHotGas = 0.0;
    Gal[p].MetalsEjectedMass = 0.0;
    Gal[p].MetalsICS = 0.0;
    
    Gal[p].AccretedGasMass = 0.0;
    Gal[p].EjectedSNGasMass = 0.0;
    Gal[p].EjectedQuasarGasMass = 0.0;
    Gal[p].MaxStrippedGas = 0.0;
    
    Gal[p].TotInstabEvents = 0;
    Gal[p].TotInstabEventsGas = 0;
    Gal[p].TotInstabEventsStar = 0;
    Gal[p].TotInstabAnnuliGas = 0;
    Gal[p].TotInstabAnnuliStar = 0;
    Gal[p].FirstUnstableGas = 0;
    Gal[p].FirstUnstableStar = 0;
    
    Gal[p].DiscRadii[0] = 0.0;
    for(j=0; j<N_BINS; j++)
    {
        if(Gal[p].Vvir>0.0)
        Gal[p].DiscRadii[j+1] = DiscBinEdge[j+1] / Gal[p].Vvir;
        else
        Gal[p].DiscRadii[j+1] = DiscBinEdge[j+1]; // This is essentially just a place-holder for problem galaxies
        Gal[p].DiscGas[j] = 0.0;
        Gal[p].DiscStars[j] = 0.0;
        Gal[p].DiscGasMetals[j] = 0.0;
        Gal[p].DiscStarsMetals[j] = 0.0;
        Gal[p].DiscSFR[j] = 0.0;
        Gal[p].TotSinkGas[j] = 0.0;
        Gal[p].TotSinkStar[j] = 0.0;
    }
    
    for(step = 0; step < STEPS; step++)
    {
        Gal[p].SfrFromH2[step] = 0.0;
        Gal[p].SfrInstab[step] = 0.0;
        Gal[p].SfrMerge[step] = 0.0;
//        Gal[p].SfrDiskColdGas[step] = 0.0;
//        Gal[p].SfrDiskColdGasMetals[step] = 0.0;
//        Gal[p].SfrBulgeColdGas[step] = 0.0;
//        Gal[p].SfrBulgeColdGasMetals[step] = 0.0;
    }
    
    Gal[p].DiskScaleRadius = get_disk_radius(halonr, p);
    Gal[p].StellarDiscScaleRadius = 1.0*Gal[p].DiskScaleRadius;
    Gal[p].CoolScaleRadius = 1.0*Gal[p].DiskScaleRadius;
    Gal[p].MergTime = 999.9;
    Gal[p].Cooling = 0.0;
    Gal[p].Heating = 0.0;
    Gal[p].r_heat = 0.0;
    Gal[p].LastMajorMerger = -1.0;
    Gal[p].LastMinorMerger = -1.0;
    Gal[p].OutflowRate = 0.0;
    
    Gal[p].infallMvir = -1.0;  //infall properties
    Gal[p].infallVvir = -1.0;
    Gal[p].infallVmax = -1.0;
    
}



double get_disk_radius(int halonr, int p)
{
    // See Mo, Mao & White (1998) eq12, and using a Bullock style lambda.
    double SpinMagnitude, SpinParameter, radius;
    
    SpinMagnitude = sqrt(Halo[halonr].Spin[0] * Halo[halonr].Spin[0] +
                         Halo[halonr].Spin[1] * Halo[halonr].Spin[1] + Halo[halonr].Spin[2] * Halo[halonr].Spin[2]);
    
    if(SpinMagnitude==0 && Gal[p].DiskScaleRadius>=0) return Gal[p].DiskScaleRadius;
    
    SpinParameter = SpinMagnitude / (1.414 * Gal[p].Vvir * Gal[p].Rvir);
    
    radius = (SpinParameter / 1.414) * Gal[p].Rvir;
    
    if(radius <= 0.0 || radius!=radius || radius==INFINITY)
    {
        if(!(SpinMagnitude>0)) printf("ERROR: Halo with SpinMagntiude of 0\n");
        radius = 0.0;
    }
    assert(radius>0);
    
    return radius;
    
}



double get_metallicity(double gas, double metals)
{
    double metallicity;
    
//    if(metals>gas)
//    printf("get_metallicity report: metals, gas/stars = %e, %e\n", metals, gas);
    
    if(gas > 0.0 && metals > 0.0)
    {
        metallicity = metals / gas;
        if(metallicity < 1.0)
            return metallicity;
        else
            return 1.0;
    }
    else
        return 0.0;
    
}



double dmax(double x, double y)
{
    if(x > y)
        return x;
    else
        return y;
}


double sqr(double x)
{
    return x*x;
}


double cube(double x)
{
    return x*x*x;
}


double exp_f(double x)
{
    // Perform a faster (but less precise) exponential
    x = 1.0 + x / 256.0;
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x; x *= x; x *= x;
    return x;
}

double get_virial_mass(int halonr, int p)
{
    if(halonr == Halo[halonr].FirstHaloInFOFgroup && Halo[halonr].Mvir > 0.0)
        return Halo[halonr].Mvir;   /* take spherical overdensity mass estimate */
    else if(Halo[halonr].Len>0)
        return Halo[halonr].Len * PartMass;
    else if(p!=-1)
        return Gal[p].StellarMass + Gal[p].ColdGas + Gal[p].HotGas + Gal[p].BlackHoleMass + Gal[p].ICS;
    else
        return 0.0; // This completes the logic but shouldn't happen
}



double get_virial_velocity(int halonr, int p)
{
    double Rvir;
    
    Rvir = get_virial_radius(halonr, p);
    
    if(Rvir > 0.0)
        return sqrt(G * get_virial_mass(halonr, p) / Rvir);
    else
        return 0.0;
}



double get_virial_radius(int halonr, int p)
{
    // return Halo[halonr].Rvir;  // Used for Bolshoi
    
    double zplus1, hubble_of_z_sq, rhocrit, fac;
    
    zplus1 = 1 + ZZ[Halo[halonr].SnapNum];
    hubble_of_z_sq =
    Hubble * Hubble *(Omega * zplus1 * zplus1 * zplus1 + (1 - Omega - OmegaLambda) * zplus1 * zplus1 +
                      OmegaLambda);
    
    rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * G);
    fac = 1 / (200 * 4 * M_PI / 3.0 * rhocrit);
    
    return cbrt(get_virial_mass(halonr, p) * fac);
}


double get_disc_gas(int p)
{
    double DiscGasSum, DiscMetalsSum;
    int l, fix;
    
    DiscGasSum = DiscMetalsSum = 0.0;
    fix = 1;
    for(l=N_BINS-1; l>=0; l--)
    {
        if(Gal[p].DiscGas[l] < 1e-20) // Negligible gas content is negligible
        {
            Gal[p].DiscGas[l] = 0.0;
            Gal[p].DiscGasMetals[l] = 0.0;
        }
        
        DiscGasSum += Gal[p].DiscGas[l];
        DiscMetalsSum += Gal[p].DiscGasMetals[l];
    }
    
    if((Gal[p].ColdGas < 1e-15 && Gal[p].ColdGas!=0.0) || (DiscGasSum < 1e-15 && DiscGasSum!=0.0) || (Gal[p].ColdGas<=0 && Gal[p].MetalsColdGas>0.0))
    {
        Gal[p].ColdGas = 0.0;
        Gal[p].MetalsColdGas = 0.0;
        for(l=0; l<N_BINS; l++)
        {
            Gal[p].DiscGas[l] = 0.0;
            Gal[p].DiscGasMetals[l] = 0.0;
        }
        DiscGasSum = 0.0;
        DiscMetalsSum = 0.0;
        fix = 0;
    }
    
    if(Gal[p].MetalsColdGas<=0.0 && fix)
    {
        for(l=0; l<N_BINS; l++)
            Gal[p].DiscGasMetals[l] = 0.0;
        DiscMetalsSum = 0.0;
        Gal[p].MetalsColdGas = 0.0;
    }
    
    if(DiscGasSum>1.001*Gal[p].ColdGas || DiscGasSum<Gal[p].ColdGas/1.001)
    {
//        printf("get_disc_gas report ... DiscSum, ColdGas =  %e, %e\n", DiscGasSum, Gal[p].ColdGas);
//        printf("get_disc_gas report ... MetalsSum, ColdMetals =  %e, %e\n", DiscMetalsSum, Gal[p].MetalsColdGas);
        
        Gal[p].ColdGas = DiscGasSum; // Prevent small errors from blowing up
        Gal[p].MetalsColdGas = DiscMetalsSum;
        
    }
    return DiscGasSum;
}

double get_disc_stars(int p)
{
    double DiscStarSum, DiscAndBulge, dumped_mass;
    int l;
    
    DiscStarSum = 0.0;
    dumped_mass = 0.0;
    for(l=N_BINS-1; l>=0; l--)
    {
        if(Gal[p].DiscStars[l] < 1e-11) // This would be less than a single star in an annulus.  Not likely.
        {
            dumped_mass += Gal[p].DiscStars[l];
            Gal[p].DiscStars[l] = 0.0;
            Gal[p].DiscStarsMetals[l] = 0.0;
        }
        DiscStarSum += Gal[p].DiscStars[l];
    }
    
    DiscAndBulge = DiscStarSum + Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass;
    
    if(DiscAndBulge>1.001*Gal[p].StellarMass || DiscAndBulge<Gal[p].StellarMass/1.001)
    {
//        printf("get_disc_stars report: DiscAndBulge, StellarMass, dumped_mass = %e, %e, %e\n", DiscAndBulge, Gal[p].StellarMass, dumped_mass);
        Gal[p].StellarMass = DiscAndBulge; // Prevent small errors from blowing up
        if(Gal[p].StellarMass <= Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass)
        {
            for(l=0; l<N_BINS; l++)
                Gal[p].DiscStars[l] = 0.0;
            DiscStarSum = 0.0;
            Gal[p].StellarMass = Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass;
        }
    }
    
    return DiscStarSum;
}

void check_channel_stars(int p)
{
    double ChannelFrac;
    
    if(Gal[p].StellarMass>0.0)
    {
        ChannelFrac = Gal[p].StellarMass / (Gal[p].StarsFromH2+Gal[p].StarsInstability+Gal[p].StarsMergeBurst);
//        if(ChannelFrac>1.01 || ChannelFrac<0.99) printf("1. ChannelFrac, StellarMass = %e, %e\n", ChannelFrac, Gal[p].StellarMass);
        Gal[p].StarsFromH2 *= ChannelFrac;
        Gal[p].StarsInstability *= ChannelFrac;
        Gal[p].StarsMergeBurst *= ChannelFrac;
    }
    else
    {
        Gal[p].StarsFromH2 = 0.0;
        Gal[p].StarsInstability = 0.0;
        Gal[p].StarsMergeBurst = 0.0;
    }
}

double get_disc_ang_mom(int p, int type)
{
    // type=0 for gas, type=1 for stars
    double J_sum;
    int l;
    
    J_sum = 0.0;
    if(type==0)
    {
        for(l=N_BINS-1; l>=0; l--)
        J_sum += Gal[p].DiscGas[l] * (DiscBinEdge[l]+DiscBinEdge[l+1])/2.0;
    }
    else if(type==1)
    {
        for(l=N_BINS-1; l>=0; l--)
        J_sum += Gal[p].DiscStars[l] * (DiscBinEdge[l]+DiscBinEdge[l+1])/2.0;
    }
    
    return J_sum;
}

void check_ejected(int p)
{
    if(!(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass))
    {
//        printf("ejected mass, metals = %e, %e\n", Gal[p].EjectedMass, Gal[p].MetalsEjectedMass);
        if(Gal[p].EjectedMass <= 1e-10 || Gal[p].MetalsEjectedMass <= 1e-10)
        {
            Gal[p].EjectedMass = 0.0;
            Gal[p].MetalsEjectedMass = 0.0;
        }
    }
}


void update_disc_radii(int p)
{
    // Calculate the radius corresponding to an annulus edge for a given galaxy.  Calculation is iterative given a generic rotation curve format.
    int i, j, j_max;
    double left, right, tol, r_try, j_try, dif, v_try;
    double M_D, M_int, M_DM, M_CB, M_SB, M_ICS, M_hot;
    double z, a, b, c_DM, c, r_2, X, M_DM_tot, rho_const;
    double a_CB, M_CB_inf, a_SB, M_SB_inf, a_ICS, M_ICS_inf;
    double f_support, BTT, v_max;
    
    update_stellardisc_scaleradius(p); // need this at the start, as disc scale radii are part of this calculation

    tol = 1e-3;
    j_max = 100;
    
    // Determine the distribution of dark matter in the halo =====
    M_DM_tot = Gal[p].Mvir - Gal[p].HotGas - Gal[p].ColdGas - Gal[p].StellarMass - Gal[p].ICS - Gal[p].BlackHoleMass; // One may want to include Ejected Gas in this too
    
    if(M_DM_tot < 0.0) M_DM_tot = 0.0;
    
    X = log10(Gal[p].StellarMass/Gal[p].Mvir);
    
    z = ZZ[Gal[p].SnapNum];
    if(z>5.0) z=5.0;
    a = 0.520 + (0.905-0.520)*exp_f(-0.617*pow(z,1.21)); // Dutton & Maccio 2014
    b = -0.101 + 0.026*z; // Dutton & Maccio 2014
    c_DM = pow(10.0, a+b*log10(Gal[p].Mvir*UnitMass_in_g/(SOLAR_MASS*1e12))); // Dutton & Maccio 2014
    if(Gal[p].Type==0 && Gal[p].StellarMass>0 && Gal[p].Mvir>0)
        c = c_DM * (1.0 + 3e-5*exp_f(3.4*(X+4.5))); // Di Cintio et al 2014b
    else
        c = 1.0*c_DM; // Should only happen for satellite-satellite mergers, where X cannot be trusted
    r_2 = Gal[p].Rvir / c; // Di Cintio et al 2014b
    rho_const = M_DM_tot / (log((Gal[p].Rvir+r_2)/r_2) - Gal[p].Rvir/(Gal[p].Rvir+r_2));
    // ===========================================================
    
    // Determine distribution for bulge and ICS ==================
    a_SB = 0.2 * Gal[p].StellarDiscScaleRadius / (1.0 + sqrt(0.5)); // Fisher & Drory (2008)
    M_SB_inf = Gal[p].SecularBulgeMass * sqr((Gal[p].Rvir+a_SB)/Gal[p].Rvir);
    
    a_CB = pow(10.0, (log10(Gal[p].ClassicalBulgeMass*UnitMass_in_g/SOLAR_MASS/Hubble_h)-10.21)/1.13) * (CM_PER_MPC/1e3) / UnitLength_in_cm * Hubble_h; // Sofue 2015
    M_CB_inf = Gal[p].ClassicalBulgeMass * sqr((Gal[p].Rvir+a_CB)/Gal[p].Rvir);
    
    a_ICS = 0.0;
    if(Gal[p].ClassicalBulgeMass>0.0)
        a_ICS = 13.0 * a_CB; // Gonzalez et al (2005)
    else if(a_SB>0.0)
        a_ICS = 13.0 * a_SB; // Feeding Fisher & Drory (2008) relation into Gonzalez et al (2005)
    else if(Gal[p].ICS > 0.0)
    {
        printf("Issue with ICS size, ICS = %e\n", Gal[p].ICS);
    }
    M_ICS_inf = Gal[p].ICS * sqr((Gal[p].Rvir+a_ICS)/Gal[p].Rvir);
    // ===========================================================
    
    M_D = 0.0;
    left = 0.0;
    
    BTT = (Gal[p].ClassicalBulgeMass+Gal[p].SecularBulgeMass)/(Gal[p].StellarMass+Gal[p].ColdGas);

    if(Gal[p].Vmax > Gal[p].Vvir)
        v_max = 2.0*Gal[p].Vmax;
    else if(Gal[p].Vvir > 0.0)
        v_max = 2.0*Gal[p].Vvir;
    else
    {
        v_max = 500.0; // Randomly large value here...
        printf("Set v_max to %.1f\n\n", v_max);
    }
    
    if(Gal[p].Mvir>0.0 && BTT<1.0)
    {
        double two_inv_Vvir = 2.0/Gal[p].Vvir;
        double inv_r_2 = 1.0/r_2;
        double hot_fraction = Gal[p].HotGas/Gal[p].Rvir;
        double exponent_support = -3.0*(1.0-BTT)/Gal[p].StellarDiscScaleRadius;
        
        for(i=1; i<N_BINS+1; i++)
        {
            right = DiscBinEdge[i] * two_inv_Vvir;
            if(right<Gal[p].Rvir) right = Gal[p].Rvir;
            if(right<8.0*left) right = 8.0*left;
            M_D += Gal[p].DiscStars[i-1] + Gal[p].DiscGas[i-1];
            
            double inv_DiscBinEdge = 1.0/DiscBinEdge[i];
            
            for(j=0; j<j_max; j++)
            {
                r_try = (left+right)*0.5;
                
                if(Gal[p].Mvir>0.0)
                    M_DM = rho_const * (log((r_try+r_2)*inv_r_2) - r_try/(r_try+r_2));
                else
                    M_DM = 0.0;
                M_SB = M_SB_inf * sqr(r_try/(r_try + a_SB));
                M_CB = M_CB_inf * sqr(r_try/(r_try + a_CB));
                M_ICS = M_ICS_inf * sqr(r_try/(r_try + a_ICS));
                M_hot = hot_fraction * r_try;
                M_int = M_DM + M_D + M_CB + M_SB + M_ICS + M_hot + Gal[p].BlackHoleMass;
                
                
                f_support = 1.0 - exp_f(exponent_support*r_try); // Fraction of support in rotation
                v_try = sqrt(G*M_int*f_support/r_try);
                if(v_try<v_max || v_max <= 0)
                    j_try = r_try*v_try;
                else
                    j_try = r_try*v_max;
                dif = j_try*inv_DiscBinEdge - 1.0;
                
                if(j==j_max-1) printf("Iterations maxed out in update_disc_radii\n");
                
                // Found correct r (within tolerance)
                if(fabs(dif) <= tol || (right-left)/right <= tol)
                    break;
                
                // Reset boundaries for next guess
                if(dif>0)
                    right = r_try;
                else
                    left = r_try;
            }
            assert(r_try==r_try && r_try!=INFINITY);
            Gal[p].DiscRadii[i] = r_try;
            left = r_try;
        }
    }
    
    update_stellardisc_scaleradius(p);
}


void update_stellardisc_scaleradius(int p)
{
    int i, j;
    double SMcum_norm[N_BINS+1];
    double DiscStarSum = get_disc_stars(p);
    
    SMcum_norm[0] = 0.0;
    // Can't update the disc size if there are no stars in the disc
    if(DiscStarSum>0.0)
    {
        for(i=1; i<N_BINS+1; i++)
        {
            if(i==0)
                SMcum_norm[i] = (Gal[p].DiscStars[i-1]/DiscStarSum);
            else
                SMcum_norm[i] = (Gal[p].DiscStars[i-1]/DiscStarSum) + SMcum_norm[i-1];
            if(SMcum_norm[i] >= 0.5)
                break;
        }

        for(j=i; j<N_BINS+1; j++)
        {
            SMcum_norm[j] = (Gal[p].DiscStars[j-1]/DiscStarSum) + SMcum_norm[j-1];
            if(SMcum_norm[j] >= 0.9)
                break;
        }

        // These are exponential scale radii that correspond to the actual r50 and r90 values of the disc
        double Rscale50 = (Gal[p].DiscRadii[i]*(0.5-SMcum_norm[i-1]) + Gal[p].DiscRadii[i-1]*(SMcum_norm[i]-0.5)) / (SMcum_norm[i]-SMcum_norm[i-1]) * 0.59581;
        double Rscale90 = (Gal[p].DiscRadii[j]*(0.9-SMcum_norm[j-1]) + Gal[p].DiscRadii[j-1]*(SMcum_norm[j]-0.9)) / (SMcum_norm[j]-SMcum_norm[j-1]) * 0.25709;
        
        // Take a weighted average of those scale radii for the disc's actual value.
        Gal[p].StellarDiscScaleRadius = (Rscale50 + RadiusWeight*Rscale90) / (RadiusWeight+1.0); 
    }
    
    if(Gal[p].StellarDiscScaleRadius<=0.0)
    {
        printf("BUG: StellarDiscScaleRadius reassigned from %e to DiskScale Radius\n", Gal[p].StellarDiscScaleRadius);
        Gal[p].StellarDiscScaleRadius = 1.0 * Gal[p].DiskScaleRadius; // Some functions still need a number for the scale radius even if there aren't any stars actually in the disc.
    }
}
