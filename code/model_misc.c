#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



void init_galaxy(int p, int halonr)
{
    int j, k, step;
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
    Gal[p].mergeIntoGalaxyNr = -1;
    Gal[p].RootID = -1;
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
    
    Gal[p].HaloScaleRadius = pow(10.0, 0.52 + (-0.101 + 0.026*ZZ[Halo[halonr].SnapNum])*log10(Gal[p].Mvir*UnitMass_in_g/(SOLAR_MASS*1e12)) ); // Initialisation based on Dutton and Maccio (2014). Gets updated as part of update_disc_radii()
    
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
        Gal[p].Potential[j+1] = NFW_potential(p, Gal[p].DiscRadii[j+1]);
        Gal[p].VelDispStars[j] = 0.0;
    }
    Gal[p].Potential[0] = Gal[p].Potential[1];

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
    Gal[p].GasDiscScaleRadius = 1.0*Gal[p].DiskScaleRadius;
    Gal[p].CoolScaleRadius = 1.0*Gal[p].DiskScaleRadius;
    Gal[p].MergTime = 999.9;
    Gal[p].Cooling = 0.0;
    Gal[p].Heating = 0.0;
    Gal[p].MaxRadioModeAccretionRate = 0.0;
    Gal[p].LastMajorMerger = -1.0;
    Gal[p].LastMinorMerger = -1.0;
    Gal[p].SNreheatRate = 0.0;
    Gal[p].SNejectRate = 0.0;
    Gal[p].c_beta = MIN_C_BETA;
    Gal[p].R2_hot_av = sqr(Gal[p].Rvir*0.5);
    
    Gal[p].infallMvir = -1.0;  //infall properties
    Gal[p].infallVvir = -1.0;
    Gal[p].infallVmax = -1.0;
    
    for(k=0; k<N_AGE_BINS; k++)
    {
        Gal[p].ClassicalBulgeMassAge[k] = 0.0;
        Gal[p].SecularBulgeMassAge[k] = 0.0;
        Gal[p].ClassicalMetalsBulgeMassAge[k] = 0.0;
        Gal[p].SecularMetalsBulgeMassAge[k] = 0.0;
        Gal[p].ICS_Age[k] = 0.0;
        Gal[p].MetalsICS_Age[k] = 0.0;
        
        for(j=0; j<N_BINS; j++)
        {
            Gal[p].DiscStarsAge[j][k] = 0.0;
            Gal[p].DiscStarsMetalsAge[j][k] = 0.0;
            Gal[p].VelDispStarsAge[j][k] = 0.0;
        }
    }

//    update_disc_radii(p);
    
    Gal[p].EjectedPotential = 0.0;
    Gal[p].HotGasPotential = NFW_potential(p, 0.5*Gal[p].Rvir);
    Gal[p].prevHotGasPotential = 0.0;
    Gal[p].prevEjectedPotential = 0.0;
    Gal[p].ReincTime = 0.0;
    Gal[p].ReincTimeFresh = 0.0;
    
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
        printf("radius, SpinMagnitude, Vvir, Rvir = %e, %e, %e %e\n", radius, SpinMagnitude, Gal[p].Vvir, Gal[p].Rvir);
        printf("halo Mvirs = %e, %e, %e, %i\n", Halo[halonr].Mvir1, Halo[halonr].Mvir2, Halo[halonr].Mvir3, Halo[halonr].Len);
        radius = 0.0;
    }
    assert(radius>0);
    
    return radius;
    
}



double get_metallicity(double gas, double metals)
{
    double metallicity;
    
    if(metals>gas)
    printf("get_metallicity report: metals, gas/stars = %e, %e\n", metals, gas);
    
    if(gas > 0.0 && metals > 0.0)
    {
        metallicity = metals / gas;
        if(metallicity < 1.0)
            return dmax(metallicity, BIG_BANG_METALLICITY);
        else
            return 1.0;
    }
    else
        return 0.0;
    
}


double dmin(double x, double y)
{
    if(x < y)
        return x;
    else
        return y;
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
    x = 1.0 + x * 0.00390625; // 0.00390625 = 1/256
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x; x *= x; x *= x;
    return x;
}

double get_virial_mass(int halonr, int p)
{
    if(halonr == Halo[halonr].FirstHaloInFOFgroup && Halo[halonr].Mvir1 > 0.0 && MvirDefinition==1)
    {
        return Halo[halonr].Mvir1;
    }
    else if(halonr == Halo[halonr].FirstHaloInFOFgroup && Halo[halonr].Mvir2 > 0.0 && MvirDefinition==2)
    {
        return Halo[halonr].Mvir2;
    }
    else if(halonr == Halo[halonr].FirstHaloInFOFgroup && Halo[halonr].Mvir3 > 0.0 && MvirDefinition==3)
    {
        return Halo[halonr].Mvir3;
    }
    else if(Halo[halonr].Len>0)
    {
        return Halo[halonr].Len * PartMass;
    }
    else if(p!=-1)
    {
        return Gal[p].StellarMass + Gal[p].ColdGas + Gal[p].HotGas + Gal[p].BlackHoleMass + Gal[p].ICS;
    }
    else
        return 0.0;
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
    
    double hubble_of_z_sq, rhocrit, fac;
    
    hubble_of_z_sq = Hubble_sqr_z(Halo[halonr].SnapNum);
    rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * G);
    fac = 1 / (200 * 4 * M_PI / 3.0 * rhocrit);
    
    return cbrt(get_virial_mass(halonr, p) * fac);
}

double Hubble_sqr_z(int snap)
{
    // Calculate the square of the Hubble parameter at input snapshot
    double zplus1 = 1 + ZZ[snap];
    return sqr(Hubble) * (Omega * cube(zplus1) + (1 - Omega - OmegaLambda) * sqr(zplus1) + OmegaLambda);
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
        else if(Gal[p].DiscGasMetals[l] < 1e-20)
            Gal[p].DiscGasMetals[l] = 0.0;
        
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
    
    if(DiscGasSum>1.001*Gal[p].ColdGas || DiscGasSum<Gal[p].ColdGas*0.999)
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
    double DiscStarSum, DiscAndBulge, AgeSum; //, dumped_mass
    int l, k;
    
    DiscStarSum = 0.0;
//    dumped_mass = 0.0;
    for(l=N_BINS-1; l>=0; l--)
    {
        if(Gal[p].DiscStars[l]<1e-11) // This would be less than a single star in an annulus.  Not likely.
        {
//            dumped_mass += Gal[p].DiscStars[l];
            Gal[p].DiscStars[l] = 0.0;
            Gal[p].DiscStarsMetals[l] = 0.0;
            
            if(AgeStructOut>0)
            {
                for(k=0; k<N_AGE_BINS; k++)
                {
                    Gal[p].DiscStarsAge[l][k] = 0.0;
                    Gal[p].DiscStarsMetalsAge[l][k] = 0.0;
                }
            }
        }
        else if(AgeStructOut>0)
        {
            AgeSum = 0.0;
            for(k=0; k<N_AGE_BINS; k++) AgeSum += Gal[p].DiscStarsAge[l][k];
            if(AgeSum>1.001*Gal[p].DiscStars[l] || AgeSum<0.999*Gal[p].DiscStars[l])
            {
                printf("get_disc_stars age report: l, AgeSum, Gal[p].DiscStars[l]: %i, %e, %e\n", l, AgeSum, Gal[p].DiscStars[l]);
                assert(AgeSum<=1.01*Gal[p].DiscStars[l] && AgeSum>=0.99*Gal[p].DiscStars[l]);
                Gal[p].DiscStars[l] = 1.0*AgeSum;
            }
        }
        DiscStarSum += Gal[p].DiscStars[l];
    }
    
    DiscAndBulge = DiscStarSum + Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass;
    
    if(DiscAndBulge>1.001*Gal[p].StellarMass || DiscAndBulge<0.999*Gal[p].StellarMass)
    {
//        printf("get_disc_stars report: DiscAndBulge, StellarMass, dumped_mass = %e, %e, %e\n", DiscAndBulge, Gal[p].StellarMass, dumped_mass);
        Gal[p].StellarMass = DiscAndBulge; // Prevent small errors from blowing up
        if(Gal[p].StellarMass <= Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass)
        {
            for(l=0; l<N_BINS; l++)
            {
                Gal[p].DiscStars[l] = 0.0;
                Gal[p].DiscStarsMetals[l] = 0.0;
            
                if(AgeStructOut>0)
                {
                    for(k=0; k<N_AGE_BINS; k++)
                    {
                        Gal[p].DiscStarsAge[l][k] = 0.0;
                        Gal[p].DiscStarsMetalsAge[l][k] = 0.0;
                    }
                }
            }
            DiscStarSum = 0.0;
            Gal[p].StellarMass = Gal[p].ClassicalBulgeMass + Gal[p].SecularBulgeMass;
        }
    }
    
    return DiscStarSum;
}

void check_channel_stars(int p)
{
    if(DelayedFeedbackOn>=1) return;
    
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
    
    // Check the age breakdown of bulges
    if(AgeStructOut>0 && Gal[p].StellarMass>0.0)
    {
        int k;
        double SecularSum=0.0, ClassicalSum=0.0;
        for(k=0; k<N_AGE_BINS; k++)
        {
            SecularSum += Gal[p].SecularBulgeMassAge[k];
            ClassicalSum += Gal[p].ClassicalBulgeMassAge[k];
        }
        
        if(SecularSum>1.01*Gal[p].SecularBulgeMass || SecularSum<0.99*Gal[p].SecularBulgeMass)
            printf("Secular bulge mass not consistent with age bins: SecularSum, Gal[p].SecularBulgeMass = %e, %e\n", SecularSum, Gal[p].SecularBulgeMass);

        if(ClassicalSum>1.01*Gal[p].ClassicalBulgeMass || ClassicalSum<0.99*Gal[p].ClassicalBulgeMass)
            printf("Classical bulge mass not consistent with age bins: ClassicalSum, Gal[p].ClassicalBulgeMass = %e, %e\n", ClassicalSum, Gal[p].ClassicalBulgeMass);
        
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
        J_sum += Gal[p].DiscGas[l] * (DiscBinEdge[l]+DiscBinEdge[l+1])*0.5;
    }
    else if(type==1)
    {
        for(l=N_BINS-1; l>=0; l--)
        J_sum += Gal[p].DiscStars[l] * (DiscBinEdge[l]+DiscBinEdge[l+1])*0.5;
    }
    
    return J_sum;
}

void check_ejected(int p)
{
    if(!(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass))
    {
//        printf("ejected mass, metals = %e, %e\n", Gal[p].EjectedMass, Gal[p].MetalsEjectedMass);
        if(Gal[p].EjectedMass <= 1e-10 || Gal[p].MetalsEjectedMass <= 1e-10*BIG_BANG_METALLICITY)
        {
            Gal[p].EjectedMass = 0.0;
            Gal[p].MetalsEjectedMass = 0.0;
        }
    }
}


static inline double v2_disc_cterm(const double x, const double c)
{
    const double cx = c*x;
    const double inv_cx = 1.0/cx;
    return (c + 4.8*c*exp(-0.35*cx - 3.5*inv_cx)) / (cx + inv_cx*inv_cx + 2*sqrt(inv_cx));
}

void update_disc_radii(int p)
{
    // Calculate the radius corresponding to an annulus edge for a given galaxy.  Calculation is iterative given a generic rotation curve format.
    int i, k;
    double M_D, M_int, M_DM, M_CB, M_SB, M_ICS, M_hot;
    double z, a, b, c_DM, c, r_2, X, M_DM_tot, rho_const;
    double a_CB, M_CB_inf, a_SB, M_SB_inf, a_ICS, M_ICS_inf;
    double f_support, BTT, v_max;
    double v2_spherical, v2_sdisc, v2_gdisc;
    
    const double inv_Rvir = 1.0/Gal[p].Rvir;
    
    update_stellardisc_scaleradius(p); // need this at the start, as disc scale radii are part of this calculation
    update_gasdisc_scaleradius(p);
    
    // Determine the distribution of dark matter in the halo =====
    M_DM_tot = Gal[p].Mvir - Gal[p].HotGas - Gal[p].ColdGas - Gal[p].StellarMass - Gal[p].ICS - Gal[p].BlackHoleMass; // One may want to include Ejected Gas in this too
//    double baryon_fraction = (Gal[p].HotGas + Gal[p].ColdGas + Gal[p].StellarMass + Gal[p].ICS + Gal[p].BlackHoleMass) / Gal[p].Mvir;
    
    if(M_DM_tot < 0.0) M_DM_tot = 0.0;
    
    X = log10(Gal[p].StellarMass/Gal[p].Mvir);
    
    z = ZZ[Halo[Gal[p].HaloNr].SnapNum];
    if(z>5.0) z=5.0;
    a = 0.520 + (0.905-0.520)*exp_f(-0.617*pow(z,1.21)); // Dutton & Maccio 2014
    b = -0.101 + 0.026*z; // Dutton & Maccio 2014
    c_DM = pow(10.0, a+b*log10(Gal[p].Mvir*UnitMass_in_g/(SOLAR_MASS*1e12))); // Dutton & Maccio 2014
    if(Gal[p].Type==0 && Gal[p].StellarMass>0 && Gal[p].Mvir>0)
        c = c_DM * (1.0 + 3e-5*exp_f(3.4*(X+4.5))); // Di Cintio et al 2014b
    else
        c = 1.0*c_DM; // Should only happen for satellite-satellite mergers, where X cannot be trusted
    r_2 = Gal[p].Rvir / c; // Di Cintio et al 2014b
    Gal[p].HaloScaleRadius = r_2; // NEW ADDITION IN 2021
    rho_const = M_DM_tot / (log((Gal[p].Rvir+r_2)/r_2) - Gal[p].Rvir/(Gal[p].Rvir+r_2));
    assert(rho_const>=0.0);
    // ===========================================================
    
    // Determine distribution for bulge and ICS ==================
    a_SB = 0.2 * Gal[p].StellarDiscScaleRadius / (1.0 + sqrt(0.5)); // Fisher & Drory (2008)
    M_SB_inf = Gal[p].SecularBulgeMass * sqr((Gal[p].Rvir+a_SB)*inv_Rvir);
    
    a_CB = pow(10.0, (log10(Gal[p].ClassicalBulgeMass*UnitMass_in_g/SOLAR_MASS/Hubble_h)-10.21)/1.13) * (CM_PER_MPC/1e3) / UnitLength_in_cm * Hubble_h; // Sofue 2015
    M_CB_inf = Gal[p].ClassicalBulgeMass * sqr((Gal[p].Rvir+a_CB)*inv_Rvir);
    
    a_ICS = 0.0;
    if(Gal[p].ClassicalBulgeMass>0.0)
        a_ICS = 13.0 * a_CB; // Gonzalez et al (2005)
    else if(a_SB>0.0)
        a_ICS = 13.0 * a_SB; // Feeding Fisher & Drory (2008) relation into Gonzalez et al (2005)
    else if(Gal[p].ICS > 0.0)
    {
        printf("Issue with ICS size, ICS = %e\n", Gal[p].ICS);
    }
    M_ICS_inf = Gal[p].ICS * sqr((Gal[p].Rvir+a_ICS)*inv_Rvir);
    // ===========================================================
    
    M_D = 0.0;
//    left = 0.0;
    
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
        const double hot_fraction = Gal[p].HotGas * inv_Rvir; // when assuming a singular isothermal sphere
        const double exponent_support = -3.0*(1.0-BTT)/Gal[p].StellarDiscScaleRadius;
        const int NUM_R_BINS=51;

        // when assuming a beta profile
        const double c_beta = Gal[p].c_beta;
        const double cb_term = 1.0/(1.0 - c_beta * atan(1.0/c_beta));
        const double hot_stuff = Gal[p].HotGas * cb_term;
        if(!(c_beta>=0))
        {
            printf("c_beta, z, SnapNum = %e, %e, %i\n", c_beta, z, Gal[p].SnapNum);
        }
        assert(c_beta>=0);
        assert(hot_stuff>=0);

        
        // set up array of radii and j from which to interpolate
        double analytic_j[NUM_R_BINS], analytic_r[NUM_R_BINS], analytic_fsupport[NUM_R_BINS], analytic_potential[NUM_R_BINS];
        analytic_j[0] = 0.0;
        analytic_r[0] = 0.0;
        analytic_potential[NUM_R_BINS-1] = 0.0; // consider making this -G*Gal[p].Mvir/Gal[p].Rvir, which would aassume that analytic_r was always <Rvir.  This assumption likely wouldn't be strictly true.  In a way, it doesn't matter, provided potentials are always used for taking differences (what else would they be useful for?)
        double r_jmax = 10.0*DiscBinEdge[N_BINS]/Gal[p].Vvir;
        if(r_jmax<Gal[p].Rvir) r_jmax = 1.0*Gal[p].Rvir;
//        if(baryon_fraction>1.0) r_jmax *= (100*baryon_fraction);
        double r = r_jmax;
        const double inv_ExponentBin = 1.0/ExponentBin;
        const double c_sdisc = Gal[p].Rvir / Gal[p].StellarDiscScaleRadius;
        const double c_gdisc = Gal[p].Rvir / Gal[p].GasDiscScaleRadius;
        const double GM_sdisc_r = G * (Gal[p].StellarMass - Gal[p].ClassicalBulgeMass - Gal[p].SecularBulgeMass) * inv_Rvir;
        const double GM_gdisc_r = G * Gal[p].ColdGas * inv_Rvir;
        double vrot, rrat, RonRvir, AvHotPotential, reincTime, intgd1, intgd2;
        int kmax;
        double Rhot = sqrt(Gal[p].R2_hot_av);
        M_DM = 1.0; // random initialisation to trigger if statement
        
        for(i=NUM_R_BINS-1; i>0; i--)
        {
            analytic_r[i] = r;
            
            // Numerical-resolution issues can arise if rrat is too small
            rrat = r/r_2;
            if(rrat<1e-8)
            {
                i += 1;
                break;  
            }
            
            RonRvir = r * inv_Rvir;
            
            // M_DM profile can momentarily dip negative at centre. Just setting all DM to be zero from that point inward
            if(M_DM>0.0)
                M_DM = rho_const * (log(rrat+1) - rrat/(rrat+1));
            else 
                M_DM = 0.0;
            M_SB = M_SB_inf * sqr(r/(r + a_SB));
            M_CB = M_CB_inf * sqr(r/(r + a_CB));
            M_ICS = M_ICS_inf * sqr(r/(r + a_ICS));
            if(HotGasProfileType==0)
                M_hot = hot_fraction * r;
            else
                M_hot = hot_stuff * (RonRvir - c_beta * atan(RonRvir/c_beta));
            M_int = M_DM + M_CB + M_SB + M_ICS + M_hot + Gal[p].BlackHoleMass;

            v2_spherical = G * M_int / r;
            v2_sdisc = GM_sdisc_r * v2_disc_cterm(RonRvir, c_sdisc);
            v2_gdisc = GM_gdisc_r * v2_disc_cterm(RonRvir, c_gdisc);
            
            f_support = 1.0 - exp(exponent_support*r); // Fraction of support in rotation
            if(f_support<1e-5) f_support = 1e-5; // Minimum safety net
            analytic_fsupport[i] = f_support; // needed later for calculating stored potential
            vrot = sqrt(f_support * (v2_spherical+v2_sdisc+v2_gdisc));
       
            analytic_j[i] = vrot * r;
            if(i<NUM_R_BINS-1)
            {
                assert(analytic_j[i]<analytic_j[i+1]);
                analytic_potential[i] = analytic_potential[i+1] - (analytic_r[i+1]-analytic_r[i]) * 0.5 * (sqr(analytic_j[i+1])/(analytic_fsupport[i+1]*cube(analytic_r[i+1])) + sqr(vrot)/(f_support*r)); // numerically integrating from "infinity" (the largest analytic_r) down to zero, step by step
            }
            
            if((i==NUM_R_BINS-1) && analytic_j[i]<DiscBinEdge[N_BINS]) // Need to expand range for interpolation!
            {
                r *= 10.0; 
                i += 1;
            }
            else
                r *= inv_ExponentBin;
   
        }
        analytic_potential[0] = analytic_potential[1]; // assumes the potential flattens at the centre, as internal mass goes to zero as radius goes to zero


        
        gsl_interp_accel *acc = gsl_interp_accel_alloc();
//        printf("i=%i\,", i);
        if(i>0) // reducing bins to avoid numerical issues with oversampling centre
        {
            const int NUM_R_BINS_REDUCED = NUM_R_BINS - i;
            double analytic_j_reduced[NUM_R_BINS_REDUCED], analytic_r_reduced[NUM_R_BINS_REDUCED], analytic_potential_reduced[NUM_R_BINS_REDUCED];
            for(k=0; k<NUM_R_BINS_REDUCED; k++)
            {
                analytic_j_reduced[k] = 1.0*analytic_j[k+i];
                analytic_r_reduced[k] = 1.0*analytic_r[k+i];
                analytic_potential_reduced[k] = 1.0*analytic_potential[k+i];
            }
            gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, NUM_R_BINS_REDUCED);
            gsl_spline_init(spline, analytic_j_reduced, analytic_r_reduced, NUM_R_BINS_REDUCED);

            gsl_spline *spline2 = gsl_spline_alloc(gsl_interp_cspline, NUM_R_BINS_REDUCED);
            gsl_spline_init(spline2, analytic_r_reduced, analytic_potential_reduced, NUM_R_BINS_REDUCED);
            
            for(k=1; k<N_BINS+1; k++)
            {
                if(DiscBinEdge[k]<analytic_j_reduced[0]) printf("If this assert statement is triggered, it's probably because the analytic_r and analytic_j arrays don't go deep enough.  Changing the rrat threshold for the break statement in the previous loop or increasing NUM_R_BINS might help.  First, check that everything in the parameter file is accurate, especially quantities like particle mass.");
                assert(DiscBinEdge[k] >= analytic_j_reduced[0]);
                assert(DiscBinEdge[k] <= analytic_j_reduced[NUM_R_BINS_REDUCED-1]);
                Gal[p].DiscRadii[k] = gsl_spline_eval(spline, DiscBinEdge[k], acc);
                Gal[p].Potential[k] = gsl_spline_eval(spline2, Gal[p].DiscRadii[k], acc);
            }
            
            AvHotPotential = 0.0;
            for(k=0; k<NUM_R_BINS_REDUCED-1; k++)
            {
                AvHotPotential += 0.5*(analytic_potential_reduced[k] + analytic_potential_reduced[k+1]) * cb_term * ((analytic_r_reduced[k+1]*inv_Rvir - c_beta*atan(analytic_r_reduced[k+1]*inv_Rvir/c_beta)) - (analytic_r_reduced[k]*inv_Rvir - c_beta*atan(analytic_r_reduced[k]*inv_Rvir/c_beta)));
                assert(AvHotPotential<=0);
                if(analytic_r_reduced[k+1] > Gal[p].Rvir) break;
            }
            
//            Gal[p].HotGasPotential = gsl_spline_eval(spline2, 0.5*Gal[p].Rvir, acc);
            Gal[p].HotGasPotential = AvHotPotential;
            Gal[p].EjectedPotential = gsl_spline_eval(spline2, 1.0*Gal[p].Rvir, acc);
            
            // calculate reioncorpotation time for freshly ejected gas, based on the current potential
            if(ReincorpotationModel==3)
            {
                reincTime = 0.0;
                kmax = k+1;
                intgd1 = 1.0 / sqrt(2*(Gal[p].EjectedPotential - analytic_potential_reduced[kmax]));
                for(k=kmax; k>0; k--)
                {
                    intgd2 = 1.0 / sqrt(2*(Gal[p].EjectedPotential - analytic_potential_reduced[k-1]));
                    reincTime += (analytic_r_reduced[k] - analytic_r_reduced[k-1]) * 0.5 * (intgd1 + intgd2);
                    if(analytic_r_reduced[k-1] <= Rhot) break;
                    intgd1 = 0.0 + intgd2;
                }
                Gal[p].ReincTimeFresh = reincTime;
            }
            
            gsl_spline_free (spline);
            gsl_spline_free (spline2);
            gsl_interp_accel_free (acc);

        }
        else
        {
            gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, NUM_R_BINS);
            gsl_spline_init(spline, analytic_j, analytic_r, NUM_R_BINS);

            gsl_spline *spline2 = gsl_spline_alloc(gsl_interp_cspline, NUM_R_BINS);
            gsl_spline_init(spline2, analytic_r, analytic_potential, NUM_R_BINS);
            
            for(k=1; k<N_BINS+1; k++)
            {
                assert(DiscBinEdge[k] >= analytic_j[0]);
                assert(DiscBinEdge[k] <= analytic_j[NUM_R_BINS-1]);
                Gal[p].DiscRadii[k] = gsl_spline_eval(spline, DiscBinEdge[k], acc);
                Gal[p].Potential[k] = gsl_spline_eval(spline2, Gal[p].DiscRadii[k], acc);
            }
            
            AvHotPotential = 0.0;
            for(k=0; k<NUM_R_BINS-1; k++)
            {
                AvHotPotential += 0.5*(analytic_potential[k] + analytic_potential[k+1]) * cb_term * ((analytic_r[k+1]*inv_Rvir - c_beta*atan(analytic_r[k+1]*inv_Rvir/c_beta)) - (analytic_r[k]*inv_Rvir - c_beta*atan(analytic_r[k]*inv_Rvir/c_beta)));
                assert(AvHotPotential<=0);
                if(analytic_r[k+1] > Gal[p].Rvir) break;
            }
            Gal[p].HotGasPotential = AvHotPotential;
            
            
//            Gal[p].HotGasPotential = dmax(analytic_potential[0], gsl_spline_eval(spline2, 0.5*Gal[p].Rvir, acc));
            Gal[p].EjectedPotential = dmax(0.9999*Gal[p].HotGasPotential, gsl_spline_eval(spline2, Gal[p].Rvir, acc)); // ensures the ejected potential mass is always a little higher (less negative, i.e. closer to zero) than hot (which is should be by definition).  This is only needed in niche instances where analytic_r doesn't probe deep enough
            
//            if(Gal[p].EjectedPotential < Gal[p].HotGasPotential)
//                for(k=0; k<NUM_R_BINS; k++) printf("r, Potential = %e, %e\n", analytic_r[k], analytic_potential[k]);

            // calculate reioncorpotation time for freshly ejected gas, based on the current potential
            if(ReincorpotationModel==3)
            {
                reincTime = 0.0;
                kmax = k+1;
                intgd1 = 1.0 / sqrt(2*(Gal[p].EjectedPotential - analytic_potential[kmax]));
                for(k=kmax; k>0; k--)
                {
                    intgd2 = 1.0 / sqrt(2*(Gal[p].EjectedPotential - analytic_potential[k-1]));
                    reincTime += (analytic_r[k] - analytic_r[k-1]) * 0.5 * (intgd1 + intgd2);
                    if(analytic_r[k-1] <= Rhot) break;
                    intgd1 = 0.0 + intgd2;
                }
                Gal[p].ReincTimeFresh = reincTime;
            }
            
            
            gsl_spline_free (spline);
            gsl_spline_free (spline2);
            gsl_interp_accel_free (acc);
        }

    }
    
    // calculate potential energy as a function of radius (in the disc plane)
//    double potential_sum = 0.0;
//    i = NUM_R_BINS-1;
//    for(k=N_BINS; k>=0; k--)
//    {
//        while(analytic_r[i])
//            potential_sum -= 
//    }
    
    if(Gal[p].Potential[0]>=0) Gal[p].Potential[0] = Gal[p].Potential[1];
    
    
    update_stellardisc_scaleradius(p);
    // if other functionality is added for gas disc scale radius, update it here too
    
    
    if(Gal[p].HotGasPotential>=0)
    {
        printf("i = %i\n", i);
        printf("k = %i\n", k);
        printf("BTT = %e\n", BTT);
    }
    assert(Gal[p].HotGasPotential<0);
    
    
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
//        printf("BUG: StellarDiscScaleRadius reassigned from %e to DiskScale Radius -- %e, %e\n", Gal[p].StellarDiscScaleRadius, DiscStarSum, Gal[p].DiskScaleRadius);
        Gal[p].StellarDiscScaleRadius = 1.0 * Gal[p].DiskScaleRadius; // Some functions still need a number for the scale radius even if there aren't any stars actually in the disc.
    }
}


void update_gasdisc_scaleradius(int p)
{
    int i, j;
    double GMcum_norm[N_BINS+1];
    double DiscGasSum = get_disc_gas(p);
    
    GMcum_norm[0] = 0.0;
    // Can't update the disc size if there are no stars in the disc
    if(DiscGasSum>0.0)
    {
        for(i=1; i<N_BINS+1; i++)
        {
            if(i==0)
            GMcum_norm[i] = (Gal[p].DiscGas[i-1]/DiscGasSum);
            else
            GMcum_norm[i] = (Gal[p].DiscGas[i-1]/DiscGasSum) + GMcum_norm[i-1];
            if(GMcum_norm[i] >= 0.5)
            break;
        }
        
        for(j=i; j<N_BINS+1; j++)
        {
            GMcum_norm[j] = (Gal[p].DiscGas[j-1]/DiscGasSum) + GMcum_norm[j-1];
            if(GMcum_norm[j] >= 0.9)
                break;
        }
        
        // These are exponential scale radii that correspond to the actual r50 and r90 values of the disc
        double Rscale50 = (Gal[p].DiscRadii[i]*(0.5-GMcum_norm[i-1]) + Gal[p].DiscRadii[i-1]*(GMcum_norm[i]-0.5)) / (GMcum_norm[i]-GMcum_norm[i-1]) * 0.59581;
        double Rscale90 = (Gal[p].DiscRadii[j]*(0.9-GMcum_norm[j-1]) + Gal[p].DiscRadii[j-1]*(GMcum_norm[j]-0.9)) / (GMcum_norm[j]-GMcum_norm[j-1]) * 0.25709;
        
        // Take a weighted average of those scale radii for the disc's actual value.
        Gal[p].GasDiscScaleRadius = (Rscale50 + RadiusWeight*Rscale90) / (RadiusWeight+1.0); 
    }
    
    if(Gal[p].GasDiscScaleRadius<=0.0)
    {
//        printf("BUG: GasDiscScaleRadius reassigned from %e to DiskScale Radius\n", Gal[p].GasDiscScaleRadius);
        Gal[p].GasDiscScaleRadius = 1.0 * Gal[p].DiskScaleRadius; // Some functions still need a number for the scale radius even if there isn't any gas actually in the disc.
    }
}


double NFW_potential(int p, double r)
{
    double radius_sum = Gal[p].Rvir + Gal[p].HaloScaleRadius;
    double pot_energy = -Gal[p].Mvir * log(1.0 + r/Gal[p].HaloScaleRadius) / (r * (log(radius_sum/Gal[p].HaloScaleRadius) - Gal[p].Rvir/radius_sum) );
    assert(pot_energy<=0);
    return pot_energy;
}


int get_stellar_age_bin_index(double time)
{
    int k_now;
    for(k_now=0; k_now<N_AGE_BINS; k_now++)
    {
        if(time<=AgeBinEdge[k_now+1])
            break; // calculate which time bin 
    }
    return k_now;
}


void get_RecycleFraction_and_NumSNperMass(double t0, double t1, double *stellar_output)
{
    double m0, m1;//, piecewise_int, frac_already_lost, norm_multiplier;
    
    // note the swapping of labels 0 and 1, as min time is max mass and vice versa. Note that mass is in solar masses, not Dark Sage internal units!!  m0 and m1 refer to the range of initial star masses over which the IMF is integrated.
    
    assert(t1>t0); // no point running this otherwise
    
    m0 = 1.0 / pow(t1 * UnitTime_in_s / SEC_PER_MEGAYEAR * 1e-4 / Hubble_h, 0.4);
    
    // built-in assumption that stars > 50 solar collapse rapidly to black holes without returning gas to the ISM or causing feedback
    if(m0>50.0) 
    {
        stellar_output[0] = 0.0;
        stellar_output[1] = 0.0;
        return;
    }
    
    // built-in assumption that subsolar stars don't cause any outflows of gas or feedback
    if(m0<1.0) 
        m0 = 1.0; 
    
    if(t0>0)
    {
        m1 = 1.0 / pow(t0 * UnitTime_in_s / SEC_PER_MEGAYEAR * 1e-4 / Hubble_h, 0.4);
        if(m1>50.0) m1 = 50.0; 
        if(m1<1.0)
        {
            stellar_output[0] = 0.0;
            stellar_output[1] = 0.0;
            return;
        }
    }
    else
        m1 = 50.0;
    
    if(!(m1>m0))
        printf("m0, m1, t0, t1 = %e, %e, %e, %e\n", m0, m1, t0, t1);
    assert(m1>m0);
    
    double denom = (integrate_m_IMF(0.1,m1) + integrate_mremnant_IMF(m1,100.0));
    stellar_output[0] = (integrate_m_IMF(m0,m1) - integrate_mremnant_IMF(m0,m1)) / denom;
    stellar_output[1] = get_numSN_perMass(m0,m1) / (denom * SOLAR_MASS * Hubble_h) * UnitMass_in_g; // convert to internal units (inverse mass)
    return;
}


double integrate_m_IMF(double m0, double m1)
{
    // Definite integral of m*IMF between m0 and m1. Mass unit is solar (not the Dark Sage internal mass unit).
    if(m1<1)
        printf("Have not added functionality for integrate_m_IMF to take m1<1 yet! Either it is time to add this to the code, or there is a bug causing this to happen.");
    assert(m1>=1);
    
    if(m1>100)
        m1 = 100.0;
    
    if(m0<1)
    {
        if(m0>0.1)
            printf("Have not added functionality for integrate_m_IMF to take 0.1<m0<1 yet!  Resetting to 0.1. Either it is time to add this to the code, or there is a bug causing this to happen: m0=%e\n", m0);
        return 0.405 - 0.7945925666666667 * (pow(m1, -0.3) - 1.0);
    }
    else
        return -0.7945925666666667 * (pow(m1, -0.3) - pow(m0, -0.3));
}

double indef_integral_mremnant_IMF(double m, int piece)
{
    // analytic indefinite integral of mremnant*IMF. Mass unit is solar (not the Dark Sage internal mass unit).    
    
    if(m>=50 && piece==4)
        return -0.7945925666666667 * pow(m, -0.3);
    else if(m>=1 && m<=7 && piece==1)
        return -0.08141517683076922*pow(m,-1.3) - 0.0667457756*pow(m,-0.3);
    else if(m>=7 && m<=8 && piece==2)
        return -0.07683098894615384*pow(m,-1.3) - 0.08661058976666666*pow(m,-0.3);
    else if(m>=8 && m<=50 && piece==3)
        return -0.2567145215384615*pow(m,-1.3);
    else if(m<1)
        printf("Have not added functionality for indef_integral_mremnant_IMF to take m<1 yet! Either it is time to add this to the code, or there is a bug causing this to happen.");
    else
        printf("Improper use of indef_integral_mremnant_IMF()!");
    return 0.0; // place-holder; should probably intentionally crash the code here
}

double integrate_mremnant_IMF(double m0, double m1)
{
    double piecewise_int = 0.0;
    // another if statement here for when m0<1?  Or would this never come up?
    if(m0<1.0) piecewise_int += integrate_m_IMF(dmax(0.1,m0), dmin(1.0,m1));
    if(m1<=1.0) return piecewise_int;
    if(m0<7.0) piecewise_int += (indef_integral_mremnant_IMF(dmin(7.0,m1),1) - indef_integral_mremnant_IMF(dmax(1.0,m0),1));
    if(m1<=7.0) return piecewise_int;
    if(m0<8.0) piecewise_int += (indef_integral_mremnant_IMF(dmin(8.0,m1),2) - indef_integral_mremnant_IMF(dmax(7.0,m0),2));
    if(m1<=8.0) return piecewise_int;
    if(m0<50.0) piecewise_int += (indef_integral_mremnant_IMF(dmin(50.0,m1),3) - indef_integral_mremnant_IMF(dmax(8.0,m0),3));
    if(m1<=50.0) return piecewise_int;
    piecewise_int += integrate_m_IMF(dmax(50.0,m0), dmin(100.0,m1));
    return piecewise_int;
}


double get_numSN_perMass(double m0, double m1)
{
    // Intergrate the IMF between mass range to calculate the number of supernovae in that time per unit remaining mass of the stellar population from which they originate
    
    // starts the same as get_recycle_fraction
    double piecewise_int;
    if(m0<1.0) m0 = 1.0;
    if(m1>50.0) m1 = 50.0;
    assert(m1>m0);
    
    piecewise_int = 0.0;
    if(m0<8.0) piecewise_int += (0.01194933634822933 * Ratio_Ia_II * (pow(dmax(1.0,m0),-1.3) - pow(dmin(8.0,m1),-1.3)) ) ;
    if(m1<=8.0) return piecewise_int;
    piecewise_int += (0.18336751538461538 * (pow(dmax(8.0,m0),-1.3) - pow(dmin(50.0,m1),-1.3)) );
    return piecewise_int;
    
}


// Not sure this is acutally useful with my design!
//double m_returned(double m_initial)
//{
//    // calculate the mass returned to the ISM from a star of initial mass m_initial.  Assumes both input and output are in solar masses (not Dark Sage internal mass units).
//    if(m_initial>=1 && m_initial<=7)
//        return 0.444 + 0.084 * m_initial;
//    else if(m_initial>7 && m_initial<8)
//        return 0.419 + 0.109*m_initial;
//    else if(m_initial>=8 && m_initial<=50)
//        return 1.4;
//    else
//        return 0.0;
//}


double get_satellite_potential(int p, int centralgal)
{
    if(p==centralgal) return 0.0; // not a satellite in this case
    
    double r = get_satellite_radius(p, centralgal);
    double Potential, frac_low, HotRad;
//    int i;
    
    HotRad = 0.5*Gal[centralgal].Rvir;
    
    if(r >= Gal[centralgal].Rvir)
        return Gal[centralgal].EjectedPotential;
    
    else if(r <= Gal[centralgal].DiscRadii[1])
        return Gal[centralgal].Potential[1]; // First 2 Potential entries should be the same
        
    else if(r < Gal[centralgal].DiscRadii[N_BINS-1])
    {
//            printf("Interp with N: r, DiscRadii0, DiscRadiiLast, Potential0, PotentialLast = %e, %e, %e, %e, %e\n", r, Gal[centralgal].DiscRadii[0], Gal[centralgal].DiscRadii[N_BINS], Gal[centralgal].Potential[0], Gal[centralgal].Potential[N_BINS]);
//            for(i=1; i<N_BINS; i++) printf("i, Radius, Potential = %i, %e, %e\n", i, Gal[centralgal].DiscRadii[i], Gal[centralgal].Potential[i]);
        
        // interpolate satellite's potential using its position and the central's potential profile
        gsl_interp_accel *acc = gsl_interp_accel_alloc();
        gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, N_BINS);
        gsl_spline_init(spline, Gal[centralgal].DiscRadii, Gal[centralgal].Potential, N_BINS);
        Potential = gsl_spline_eval(spline, r, acc);
        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);
        return Potential;
    }
    
    else if(r < Gal[centralgal].DiscRadii[N_BINS])
    {
        // implicitly assuming that r is closer to at least 2 of DiscRadii[N_BINS], DiscRadii[N_BINS-1] and 0.5*Rvir than Rvir (this will still work if that isn't true, it just wouldn't be as accurate an interpolation as it might have otherwise been in that case)
        
        if(r >= HotRad && HotRad >= Gal[centralgal].DiscRadii[N_BINS-1])
        {
            frac_low = (r - HotRad) / (Gal[centralgal].DiscRadii[N_BINS] - HotRad);
            return frac_low*Gal[centralgal].Potential[N_BINS] + (1-frac_low)*Gal[centralgal].HotGasPotential;
        }
        else if(r < HotRad && HotRad < Gal[centralgal].DiscRadii[N_BINS])
        {
            frac_low = (r - Gal[centralgal].DiscRadii[N_BINS-1]) / (HotRad - Gal[centralgal].DiscRadii[N_BINS-1]);
            return frac_low*Gal[centralgal].HotGasPotential + (1-frac_low)*Gal[centralgal].Potential[N_BINS-1];
        }
        else
        {
            frac_low = (r - Gal[centralgal].DiscRadii[N_BINS-1]) / (Gal[centralgal].DiscRadii[N_BINS] - Gal[centralgal].DiscRadii[N_BINS-1]);
            return frac_low*Gal[centralgal].Potential[N_BINS] + (1-frac_low)*Gal[centralgal].Potential[N_BINS-1];
        }
    }
    else
    {
        if(r < HotRad)
        {
            frac_low = (r - Gal[centralgal].DiscRadii[N_BINS]) / (HotRad - Gal[centralgal].DiscRadii[N_BINS]);
            return frac_low*Gal[centralgal].HotGasPotential + (1-frac_low)*Gal[centralgal].Potential[N_BINS];
        }
        else if(HotRad >= Gal[centralgal].DiscRadii[N_BINS])
        {
            frac_low = (r - HotRad) / HotRad; // Assumption here that Rvir-HotRad = HotRad, i.e. HotRad = Rvir/2 (which is currently hard-coded)
            return frac_low*Gal[centralgal].EjectedPotential + (1-frac_low)*Gal[centralgal].HotGasPotential;
        }
        else
        {
            frac_low = (r - Gal[centralgal].DiscRadii[N_BINS]) / (Gal[centralgal].Rvir - Gal[centralgal].DiscRadii[N_BINS]);
            return frac_low*Gal[centralgal].EjectedPotential + (1-frac_low)*Gal[centralgal].Potential[N_BINS];
        }

    }
    
    
    
//    else // I probably don't need to call an interpolation function here, but not obvious that doing the interpolation manually would be any faster
//    {
//        if(Gal[centralgal].DiscRadii[N_BINS] < 0.5*Gal[centralgal].Rvir)
//        {
//            printf("Interp with 3\n");
//            double radius_arr[3], potential_arr[3];
//            radius_arr[0] = 1.0*Gal[centralgal].DiscRadii[N_BINS];
//            radius_arr[1] = 0.5*Gal[centralgal].Rvir;
//            radius_arr[2] = 1.0*Gal[centralgal].Rvir;
//            potential_arr[0] = 1.0*Gal[centralgal].Potential[N_BINS];
//            potential_arr[1] = 1.0*Gal[centralgal].HotGasPotential;
//            potential_arr[2] = 1.0*Gal[centralgal].EjectedPotential;
//            
//            gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, 3);
//            gsl_spline_init(spline, radius_arr, potential_arr, 3);
//            Potential = gsl_spline_eval(spline, r, acc);
//            gsl_spline_free (spline);
//            gsl_interp_accel_free (acc);
//            return Potential;
//            
//        }
//        else
//        {
//            printf("Interp with 2\n");
//            double radius_arr[2], potential_arr[2];
//            radius_arr[0] = 1.0*Gal[centralgal].DiscRadii[N_BINS];
//            radius_arr[1] = 1.0*Gal[centralgal].Rvir;
//            potential_arr[0] = 1.0*Gal[centralgal].Potential[N_BINS];
//            potential_arr[1] = 1.0*Gal[centralgal].EjectedPotential;
//            
//            gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, 2);
//            gsl_spline_init(spline, radius_arr, potential_arr, 2);
//            Potential = gsl_spline_eval(spline, r, acc);
//            gsl_spline_free (spline);
//            gsl_interp_accel_free (acc);
//            return Potential;
//        }
//    }
    
}

double get_satellite_radius(int p, int centralgal)
{
    int i;
    double dx;
    double r2 = 0.0;
    for(i=0; i<3; i++)
    {
        dx = fabs(Gal[p].Pos[i] - Gal[centralgal].Pos[i]);
        if(dx>HalfBoxLen) dx -= BoxLen;
        r2 += sqr(dx);
    }
    return sqrt(r2) * AA[Gal[p].SnapNum]; // returns in physical units, not comoving
}


double get_satellite_mass(int p)
{
    // 'virial mass' should always include baryons, but tidal stripping can cause this to fall below the baryon mass
    return dmax(Gal[p].Mvir, Gal[p].StellarMass + Gal[p].ColdGas + Gal[p].HotGas + Gal[p].BlackHoleMass);
}


double get_Mhost_internal(int p, int centralgal, double dr)
{
    // Get the mass of the host halo internal to a satellite
    // dr serves as a way of calculating internal mass at a radius close to but not exactly equal to that of a satellite (useful for calculating mass gradients near satellites, for example)
    
    if(p==centralgal)
    {
        printf("get_Mhost_internal() inappropriately called for a central\n");
        return 0.0;
    }
    
    double Rvir_host = Gal[centralgal].Rvir;
    double Mhost;
    double SatelliteRadius = get_satellite_radius(p, centralgal) + dr;
    double SatelliteMass = get_satellite_mass(p); 
    int i;
    
    // calculate internal mass of halo from satellite's position
    if(SatelliteRadius < Rvir_host)
    {
        double M_DM_tot, X, z, a, b, c_DM, c, r_2, rho_const, M_DM;
        M_DM_tot = Gal[centralgal].Mvir - Gal[centralgal].HotGas - Gal[centralgal].ColdGas - Gal[centralgal].StellarMass - Gal[centralgal].BlackHoleMass - SatelliteMass; // doesn't properly account for mass of other satellites -- effectively treats them like dark matter
        if(M_DM_tot > 0.0) 
        {
            X = log10(Gal[centralgal].StellarMass/Gal[centralgal].Mvir);
            z = ZZ[Gal[centralgal].SnapNum];
            if(z>5.0) z=5.0;
            a = 0.520 + (0.905-0.520)*exp(-0.617*pow(z,1.21)); // Dutton & Maccio 2014
            b = -0.101 + 0.026*z; // Dutton & Maccio 2014
            c_DM = pow(10.0, a+b*log10(Gal[centralgal].Mvir*UnitMass_in_g/(SOLAR_MASS*1e12))); // Dutton & Maccio 2014
            c = c_DM * (1.0 + 3e-5*exp(3.4*(X+4.5))); // Di Cintio et al 2014b
            r_2 = Gal[centralgal].Rvir / c; // Di Cintio et al 2014b
            rho_const = M_DM_tot / (log((Gal[centralgal].Rvir+r_2)/r_2) - Gal[centralgal].Rvir/(Gal[centralgal].Rvir+r_2));
            M_DM = rho_const * (log((SatelliteRadius+r_2)/r_2) - SatelliteRadius/(SatelliteRadius+r_2));
        }
        else
            M_DM = 0.0;
        
        double M_hot;
        if(HotGasProfileType==0)
            M_hot = Gal[centralgal].HotGas * SatelliteRadius / Rvir_host;
        else
        {
            const double c_beta = Gal[centralgal].c_beta;
            const double cb_term = 1.0/(1.0 - c_beta * atan(1.0/c_beta));
            const double hot_stuff = Gal[centralgal].HotGas * cb_term;
            const double RonRvir = SatelliteRadius / Rvir_host;
            M_hot = hot_stuff * (RonRvir - c_beta * atan(RonRvir/c_beta));
        }
        
        const double a_SB = 0.2 * Gal[centralgal].StellarDiscScaleRadius / (1.0 + sqrt(0.5)); // Fisher & Drory (2008)
        const double M_iBulge = Gal[centralgal].SecularBulgeMass * sqr((Rvir_host+a_SB)/Rvir_host) * sqr(SatelliteRadius/(SatelliteRadius + a_SB));
        
        const double a_CB = pow(10.0, (log10(Gal[centralgal].ClassicalBulgeMass*UnitMass_in_g/SOLAR_MASS/Hubble_h)-10.21)/1.13) * (CM_PER_MPC/1e3) / UnitLength_in_cm * Hubble_h; // Sofue 2015
        const double M_mBulge = Gal[centralgal].ClassicalBulgeMass * sqr((Rvir_host+a_CB)/Rvir_host) * sqr(SatelliteRadius/(SatelliteRadius + a_CB));
        
        double a_ICS = 0.0;
        if(Gal[p].ClassicalBulgeMass>0.0)
            a_ICS = 13.0 * a_CB; // Gonzalez et al (2005)
        else if(a_SB>0.0)
            a_ICS = 13.0 * a_SB; // Feeding Fisher & Drory (2008) relation into Gonzalez et al (2005)      
        const double M_ICS = Gal[centralgal].ICS * sqr((Rvir_host+a_ICS)/Rvir_host) * sqr(SatelliteRadius/(SatelliteRadius + a_ICS));
        
        // Add mass from the disc
        double M_disc = 0.0;
        for(i=0; i<N_BINS; i++)
        {
            if(Gal[centralgal].DiscRadii[i+1] <= SatelliteRadius)
                M_disc += (Gal[centralgal].DiscGas[i] + Gal[centralgal].DiscStars[i]);
            else
            {
                M_disc += ((Gal[centralgal].DiscGas[i] + Gal[centralgal].DiscStars[i]) * sqr((SatelliteRadius - Gal[centralgal].DiscRadii[i])/(Gal[centralgal].DiscRadii[i+1]-Gal[centralgal].DiscRadii[i])));
                break;
            }
        }
        
        Mhost = dmin(Gal[centralgal].Mvir, M_DM + M_hot + M_disc + M_iBulge + M_mBulge + M_ICS + Gal[centralgal].BlackHoleMass);
    }
    else
        Mhost = Gal[centralgal].Mvir;
    
    return Mhost;
}


void rotate(double *pos, double axis[3], double angle)
{ // rotate a 3-vector pos by angle about axis
    double dot = pos[0]*axis[0] + pos[1]*axis[1] + pos[2]*axis[2];
    
    double cross[3];
    cross[0] = axis[1]*pos[2] - axis[2]*pos[1];
    cross[1] = axis[2]*pos[0] - axis[0]*pos[2];
    cross[2] = axis[0]*pos[1] - axis[1]*pos[0];
    
    double cosa = cos(angle);
    double sina = sin(angle);
    
    int i;
    for(i=0; i<3; i++) pos[i] = pos[i]*cosa + cross[i]*sina + axis[i]*dot*(1-cosa);
}
