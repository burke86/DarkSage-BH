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
    
//    for(i=0; i<=N_AGE_BINS; i++) printf("i, AgeBinEdge[i], Gyr = %i, %e, %e\n", i, AgeBinEdge[i], AgeBinEdge[i]*1e-3*UnitTime_in_s/SEC_PER_MEGAYEAR/Hubble_h);
        
    // Set counts for prograde and retrograde satellite collisions
    RetroCount = 0;
    ProCount = 0;
    
    // Determine the total returned mass fraction from a population of stars
    double StellarOutput[2];
    get_RecycleFraction_and_NumSNperMass(0, 1e20, StellarOutput); // arbitrarily large number for upper bound on time
//    printf("FinalRecycleFraction, SNperMassFormed = %e, %e\n", StellarOutput[0], StellarOutput[1]);
    if(DelayedFeedbackOn>0)
        FinalRecycleFraction = 1.0*StellarOutput[0];
    else
        FinalRecycleFraction = 1.0 * RecycleFraction;
    
    // running a check on get_RecycleFraction_and_NumSNperMass() returning what I want it to
//    double test_mass = 1.0;
//    int tt;
//    double time_convert = 1.0 / (1e-3 * UnitTime_in_s / SEC_PER_MEGAYEAR / Hubble_h);
//    for(tt=0; tt<1500; tt++)
//    {
//        if(tt%30==0) continue;
//        get_RecycleFraction_and_NumSNperMass(tt*0.01*time_convert, (tt+1)*0.01*time_convert, StellarOutput);
////        get_RecycleFraction_and_NumSNperMass(0.0, 0.5*(AgeBinEdge[tt+1] - AgeBinEdge[tt]), StellarOutput);
////        printf("bin, instant age interval, return fraction = %i, %e, %e\n", tt, 0.5*(AgeBinEdge[tt+1]-AgeBinEdge[tt])/time_convert, StellarOutput[0]);
////        get_RecycleFraction_and_NumSNperMass(AgeBinEdge[tt], AgeBinEdge[tt+1], StellarOutput);
//        test_mass -= (StellarOutput[0]*test_mass);
//    }
//    printf("test_mass, 1-FinalRecycleFraction = %e, %e\n", test_mass, 1-FinalRecycleFraction);
//    assert(test_mass == 1-FinalRecycleFraction);
    
    
    // Used for SupernovaRecipeOn>3.  If DelayedFeedbackOn==1, this will be updated in core_build_model.c.  The value below assumes the instantaneous recycling (and therefore instantaneous feedback) approximation
    SNperMassFormed = 1.0*StellarOutput[1];
    
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
