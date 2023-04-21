#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#ifdef MPI
#include <mpi.h>
#endif

#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"


char bufz0[1000];
int exitfail = 1;

struct sigaction saveaction_XCPU;
volatile sig_atomic_t gotXCPU = 0;



void termination_handler(int signum)
{
    gotXCPU = 1;
    sigaction(SIGXCPU, &saveaction_XCPU, NULL);
    if(saveaction_XCPU.sa_handler != NULL)
    (*saveaction_XCPU.sa_handler) (signum);
}



void myexit(int signum)
{
#ifdef MPI
    printf("Task: %d\tnode: %s\tis exiting.\n\n\n", ThisTask, ThisNode);
#else
    printf("Exiting\n\n\n");
#endif
    exit(signum);
}



void bye()
{
#ifdef MPI
    MPI_Finalize();
    free(ThisNode);
#endif
    
    if(exitfail)
    {
        unlink(bufz0);
#ifdef MPI
        if(ThisTask == 0 && gotXCPU == 1)
        printf("Received XCPU, exiting. But we'll be back.\n");
#endif
    }
}



int main(int argc, char **argv)
{
    int filenr, tree, halonr, i, k;
    double k_float;
    struct sigaction current_XCPU;
    
    struct stat filestatus;
    FILE *fd;
#ifdef MPI
    time_t start, current;
#endif
    
#ifdef MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);
    
    ThisNode = malloc(MPI_MAX_PROCESSOR_NAME * sizeof(char));
    
    MPI_Get_processor_name(ThisNode, &nodeNameLen);
    if (nodeNameLen >= MPI_MAX_PROCESSOR_NAME) 
    {
        printf("Node name string not long enough!...\n");
        ABORT(0);
    }
    
#endif
    
    if(argc != 2)
    {
        printf("\n  usage: DARK SAGE <parameterfile>\n\n");
        ABORT(1);
    }
    
    atexit(bye);
    
    sigaction(SIGXCPU, NULL, &saveaction_XCPU);
    current_XCPU = saveaction_XCPU;
    current_XCPU.sa_handler = termination_handler;
    sigaction(SIGXCPU, &current_XCPU, NULL);
    
    read_parameter_file(argv[1]);
    init();
    srand(pow(getpid() % 100, 2.0)); // Seed with last 3 digits will be too similar when numbers approach 1000
    
    // Define the specific-angular-momentum bins used to collect disc mass
    DiscBinEdge[0] = 0.0;
    for(i=1; i<N_BINS+1; i++)
    {
        DiscBinEdge[i] = FirstBin*(CM_PER_MPC/UnitLength_in_cm/1e3)/(UnitVelocity_in_cm_per_s/1e5) *pow(ExponentBin, i-1);
    }
    
    // Define age bins (in terms of look-back time) for stellar content
    // Default bin edges match the simulation alist. Otherwise the alist will be interpolated
    int LastSnap = ListOutputSnaps[NOUT-1];
    if(N_AGE_BINS==LastSnap)
    {
        for(i=0; i<N_AGE_BINS+1; i++)
            AgeBinEdge[i] = time_to_present(ZZ[LastSnap-i]);
    }
    else
    {
        double age_snap_ratio = (double)(LastSnap)/(double)N_AGE_BINS;
        AgeBinEdge[0] = time_to_present(ZZ[LastSnap]);
        for(i=1; i<N_AGE_BINS; i++)
        {
            k_float = (double)LastSnap - (double)i*age_snap_ratio;
            k = (int)k_float;
            AgeBinEdge[i] = time_to_present(ZZ[k]);
        }
        AgeBinEdge[N_AGE_BINS] = time_to_present(ZZ[0]);
    }
            
    // Set counts for prograde and retrograde satellite collisions
    RetroCount = 0;
    ProCount = 0;
    
    // Determine the total returned mass fraction from a population of stars
    double StellarOutput[2];
    get_RecycleFraction_and_NumSNperMass(0, 1e20, StellarOutput); // arbitrarily large number for upper bound on time
//    printf("FinalRecycleFraction, SNperMassFormed = %e, %e\n", StellarOutput[0], StellarOutput[1]);
    FinalRecycleFraction = 1.0*StellarOutput[0];
    
    // when using the instantaneous recycling approximation
    // these terms are otherwise regularly updated in core_build_model.c
    if(DelayedFeedbackOn==0)
    {
        RecycleFraction = FinalRecycleFraction;
        SNperMassFormed = 1.0*StellarOutput[1];
    }
        
    HalfBoxLen = 0.5 * BoxLen; // useful to have a field of half the box length for calculating galaxy--galaxy distances
    
    // conversion factors for percentage mass radii to exponential scale radii
    DiscScalePercentConversion[0] = 1.0 / 0.531812;
    DiscScalePercentConversion[1] = 1.0 / 0.824388;
    DiscScalePercentConversion[2] = 1.0 / 1.09735;
    DiscScalePercentConversion[3] = 1.0 / 1.37642;
    DiscScalePercentConversion[4] = 1.0 / 1.67835;
    DiscScalePercentConversion[5] = 1.0 / 2.02231;
    DiscScalePercentConversion[6] = 1.0 / 2.43922;
    DiscScalePercentConversion[7] = 1.0 / 2.99431;
    DiscScalePercentConversion[8] = 1.0 / 3.88927;
    
    for(i=0; i<9; i++)
        DiscScalePercentValues[i] = 0.1 * (i+1);
    
    // calculate the UV background at each snapshot
    // using the FG09 background in the Lyman--Werner band
    double z_arr[45], UVLW_arr[31], GammaHI_arr[45];
    for(i=0; i<31; i++) z_arr[i] = 0.2*i;
    UVLW_arr[0] = 1.38e-3;
    UVLW_arr[1] = 2.38e-3;
    UVLW_arr[2] = 3.81e-3;
    UVLW_arr[3] = 5.81e-3;
    UVLW_arr[4] = 8.41e-3;
    UVLW_arr[5] = 1.15e-2;
    UVLW_arr[6] = 1.52e-2;
    UVLW_arr[7] = 1.95e-2;
    UVLW_arr[8] = 2.42e-2;
    UVLW_arr[9] = 2.95e-2;
    UVLW_arr[10] = 3.51e-2;
    UVLW_arr[11] = 4.17e-2;
    UVLW_arr[12] = 4.88e-2;
    UVLW_arr[13] = 5.63e-2;
    UVLW_arr[14] = 6.34e-2;
    UVLW_arr[15] = 7.19e-2;
    UVLW_arr[16] = 8.04e-2;
    UVLW_arr[17] = 8.89e-2;
    UVLW_arr[18] = 9.96e-2;
    UVLW_arr[19] = 1.08e-1;
    UVLW_arr[20] = 1.17e-1;
    UVLW_arr[21] = 1.26e-1;
    UVLW_arr[22] = 1.35e-1;
    UVLW_arr[23] = 1.43e-1;
    UVLW_arr[24] = 1.51e-1;
    UVLW_arr[25] = 1.58e-1;
    UVLW_arr[26] = 1.63e-1;
    UVLW_arr[27] = 1.69e-1;
    UVLW_arr[28] = 1.75e-1;
    UVLW_arr[29] = 1.79e-1;
    UVLW_arr[30] = 1.84e-1;
    
    // same again but for ionizing background, not the molecule-dissociating background
    double z_qHI[29], qHI[29];
    for(i=0; i<29; i++) z_qHI[i] = 0.25*i;
    z_qHI[0] = 0.158;
    z_qHI[1] = 0.311;
    z_qHI[2] = 0.569;
    z_qHI[3] = 0.929;
    z_qHI[4] = 1.371;
    z_qHI[5] = 1.841;
    z_qHI[6] = 2.260;
    z_qHI[7] = 2.574;
    z_qHI[8] = 2.768;
    z_qHI[9] = 2.852;
    z_qHI[10] = 2.839;
    z_qHI[11] = 2.762;
    z_qHI[12] = 2.642;
    z_qHI[13] = 2.511;
    z_qHI[14] = 2.384;
    z_qHI[15] = 2.272;
    z_qHI[16] = 2.171;
    z_qHI[17] = 2.083;
    z_qHI[18] = 2.002;
    z_qHI[19] = 1.921;
    z_qHI[20] = 1.833;
    z_qHI[21] = 1.745;
    z_qHI[22] = 1.661;
    z_qHI[23] = 1.573;
    z_qHI[24] = 1.487;
    z_qHI[25] = 1.399;
    z_qHI[26] = 1.305;
    z_qHI[27] = 1.216;
    z_qHI[28] = 1.127;
    
    int snap, i_arr;
    for(snap=0; snap<MAXSNAPS; snap++)
    {
        for(i_arr=0; i_arr<31; i_arr++)
            if(z_arr[i_arr] >= ZZ[snap]) break;
        
        if(i_arr<30)
            UVB_z[snap] = UVLW_arr[i_arr] + (UVLW_arr[i_arr+1] - UVLW_arr[i_arr]) * (ZZ[snap]-z_arr[i_arr]) / (z_arr[i_arr+1]-z_arr[i_arr]);
        else
            UVB_z[snap] = UVLW_arr[30];
    }
    
    // UV flux per SFR at distance squared in internal units
    // Value of 96.21 kpc^2/(Msun/yr) comes from Diemer, Stevens et al. (2018).
    UVMW_perSFRdensity = 96.21 * sqr(Hubble_h * 1e-3 * CM_PER_MPC / UnitLength_in_cm) / (SOLAR_MASS / UnitMass_in_g) * (SEC_PER_YEAR / UnitTime_in_s);
    Sigma_R1_fac = 50.0 * (Hubble_h * SOLAR_MASS / UnitMass_in_g) / sqr(Hubble_h * 1e-6 * CM_PER_MPC / UnitLength_in_cm);
    
    
#ifdef MPI
    // A small delay so that processors don't use the same file
    //    printf("Small delay for processors\n");
    time(&start);
    do
    time(&current);
    while(difftime(current, start) < 5.0 * ThisTask);
#endif
    
#ifdef MPI
    for(filenr = FirstFile+ThisTask; filenr <= LastFile; filenr += NTask)
#else
    for(filenr = FirstFile; filenr <= LastFile; filenr++)
#endif
    {
        sprintf(bufz0, "%s/%s.%d", SimulationDir, TreeName, filenr);
        
        // Sleep for case of running with MPI without actually enabling MPI, so processors don't do the same job!
        const unsigned int sleep_time = (10000ULL*(getpid() % 100));
        //printf("random_sleep_time = %u ns pid = %zu\n", random_sleep_time,getpid());
        usleep(sleep_time);
        if(!(fd = fopen(bufz0, "r")))
        {
            printf("-- missing tree %s ... skipping\n", bufz0);
            continue;  // tree file does not exist, move along
        }
        else
        fclose(fd);
        
        sprintf(bufz0, "%s/%s_z%1.3f_%d", OutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[0]], filenr);
        if(stat(bufz0, &filestatus) == 0)
        {
            printf("-- output for tree %s already exists ... skipping\n", bufz0);
            continue;  // output seems to already exist, dont overwrite, move along
        }
        
        if((fd = fopen(bufz0, "w"))) fclose(fd);
        
        FileNum = filenr;
        load_tree_table(filenr);
        
        
        for(tree = 0; tree < Ntrees; tree++)
        {
            
            assert(!gotXCPU);
            
            if(tree % 10000 == 0)
            {
#ifdef MPI
                printf("\ttask: %d\tnode: %s\tfile: %i\ttree: %i of %i\n", ThisTask, ThisNode, filenr, tree, Ntrees);
#else
                printf("\tfile: %i\ttree: %i of %i\n", filenr, tree, Ntrees);
#endif
                fflush(stdout);
            }
            
            TreeID = tree;
            load_tree(tree);
            
            gsl_rng_set(random_generator, filenr * 100000 + tree);
            NumGals = 0;
            GalaxyCounter = 0;
            for(halonr = 0; halonr < TreeNHalos[tree]; halonr++)
            {
                if(HaloAux[halonr].DoneFlag == 0)
                    construct_galaxies(halonr, tree);
//                if(Halo[halonr].SnapNum == Snaplistlen-1)
//                    assign_root_index(halonr);
            }
            
            
            save_galaxies(filenr, tree);
            free_galaxies_and_tree();
        }
        
        finalize_galaxy_file(filenr);
        free_tree_table();
        //      printf("\nPro v retro = %d, %d", ProCount, RetroCount);
        printf("\ndone file %d\n\n", filenr);
    }
    
    exitfail = 0;
    return 0;
}
