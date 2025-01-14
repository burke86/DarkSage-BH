#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"

#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300



void read_parameter_file(char *fname)
{
  FILE *fd;
  char buf[400], buf1[400], buf2[400], buf3[400];
  int i, j, nt = 0, done;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int errorFlag = 0;

#ifdef MPI
  if(ThisTask == 0)
    printf("\nreading parameter file:\n\n");
#endif
  
  strcpy(tag[nt], "FileNameGalaxies");
  addr[nt] = FileNameGalaxies;
  id[nt++] = STRING;

  strcpy(tag[nt], "TreeName");
  addr[nt] = TreeName;
  id[nt++] = STRING;

  strcpy(tag[nt], "OutputDir");
  addr[nt] = OutputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "FirstFile");
  addr[nt] = &FirstFile;
  id[nt++] = INT;

  strcpy(tag[nt], "LastFile");
  addr[nt] = &LastFile;
  id[nt++] = INT;

  strcpy(tag[nt], "LastSnapShotNr");
  addr[nt] = &LastSnapShotNr;
  id[nt++] = INT;

  strcpy(tag[nt], "SimulationDir");
  addr[nt] = SimulationDir;
  id[nt++] = STRING;
  
  strcpy(tag[nt], "FileWithSnapList");
  addr[nt] = FileWithSnapList;
  id[nt++] = STRING;
    
  strcpy(tag[nt], "FirstBin");
  addr[nt] = &FirstBin;
  id[nt++] = DOUBLE;
    
  strcpy(tag[nt], "ExponentBin");
  addr[nt] = &ExponentBin;
  id[nt++] = DOUBLE;
    
  strcpy(tag[nt], "ThreshMajorMerger");
  addr[nt] = &ThreshMajorMerger;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "CoolingScaleSlope");
  addr[nt] = &CoolingScaleSlope;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "CoolingScaleConst");
  addr[nt] = &CoolingScaleConst;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
  addr[nt] = &UnitVelocity_in_cm_per_s;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "UnitLength_in_cm");
  addr[nt] = &UnitLength_in_cm;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "UnitMass_in_g");
  addr[nt] = &UnitMass_in_g;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Hubble_h");
  addr[nt] = &Hubble_h;
  id[nt++] = DOUBLE;
    
  strcpy(tag[nt], "BoxLen");
  addr[nt] = &BoxLen;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ReionizationOn");
  addr[nt] = &ReionizationOn;
  id[nt++] = INT;

  strcpy(tag[nt], "SupernovaRecipeOn");
  addr[nt] = &SupernovaRecipeOn;
  id[nt++] = INT;

  strcpy(tag[nt], "DiskInstabilityOn");
  addr[nt] = &DiskInstabilityOn;
  id[nt++] = INT;

  strcpy(tag[nt], "SFprescription");
  addr[nt] = &SFprescription;
  id[nt++] = INT;

  strcpy(tag[nt], "H2prescription");
  addr[nt] = &H2prescription;
  id[nt++] = INT;
    
  strcpy(tag[nt], "AGNrecipeOn");
  addr[nt] = &AGNrecipeOn;
  id[nt++] = INT;

  strcpy(tag[nt], "GasPrecessionOn");
  addr[nt] = &GasPrecessionOn;
  id[nt++] = INT;
    
  strcpy(tag[nt], "RamPressureOn");
  addr[nt] = &RamPressureOn;
  id[nt++] = INT;
    
  strcpy(tag[nt], "HotStripOn");
  addr[nt] = &HotStripOn;
  id[nt++] = INT;

  strcpy(tag[nt], "CoolingExponentialRadiusOn");
  addr[nt] = &CoolingExponentialRadiusOn;
  id[nt++] = INT;
    
  strcpy(tag[nt], "MvirDefinition");
  addr[nt] = &MvirDefinition;
  id[nt++] = INT;

  strcpy(tag[nt], "AgeStructOut");
  addr[nt] = &AgeStructOut;
  id[nt++] = INT;

  strcpy(tag[nt], "DelayedFeedbackOn");
  addr[nt] = &DelayedFeedbackOn;
  id[nt++] = INT;

  strcpy(tag[nt], "HotGasProfileType");
  addr[nt] = &HotGasProfileType;
  id[nt++] = INT;

  strcpy(tag[nt], "MergeTimeScaleForm");
  addr[nt] = &MergeTimeScaleForm;
  id[nt++] = INT;
    
  strcpy(tag[nt], "BaryonFrac");
  addr[nt] = &BaryonFrac;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Omega");
  addr[nt] = &Omega;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "OmegaLambda");
  addr[nt] = &OmegaLambda;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "PartMass");
  addr[nt] = &PartMass;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "EnergySN");
  addr[nt] = &EnergySN;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "SeedPowerLawIndex");
  addr[nt] = &SeedPowerLawIndex;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "SeedProb");
  addr[nt] = &SeedProb;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Yield");
  addr[nt] = &Yield;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "SfrEfficiency");
  addr[nt] = &SfrEfficiency;
  id[nt++] = DOUBLE;
        
  strcpy(tag[nt], "RadioModeEfficiency");
  addr[nt] = &RadioModeEfficiency;
  id[nt++] = DOUBLE;
    
  strcpy(tag[nt], "RadiativeEfficiency");
  addr[nt] = &RadiativeEfficiency;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt], "H2FractionFactor");
  addr[nt] = &H2FractionFactor;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt], "H2FractionExponent");
  addr[nt] = &H2FractionExponent;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt], "ClumpFactor");
  addr[nt] = &ClumpFactor;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt], "ClumpExponent");
  addr[nt] = &ClumpExponent;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt], "QTotMin");
  addr[nt] = &QTotMin;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt], "GasSinkRate");
  addr[nt] = &GasSinkRate;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt], "ThetaThresh");
  addr[nt] = &ThetaThresh;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt], "DegPerTdyn");
  addr[nt] = &DegPerTdyn;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt], "Reionization_z0");
  addr[nt] = &Reionization_z0;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Reionization_zr");
  addr[nt] = &Reionization_zr;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ThresholdSatDisruption");
  addr[nt] = &ThresholdSatDisruption;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt], "NumOutputs");
  addr[nt] = &NOUT;
  id[nt++] = INT;

  if((fd = fopen(fname, "r")))
  {
    while(!feof(fd))
    {
      *buf = 0;
      fgets(buf, 400, fd);
      if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
        continue;

      if(buf1[0] == '%' || buf1[0] == '-')
        continue;

      for(i = 0, j = -1; i < nt; i++)
        if(strcmp(buf1, tag[i]) == 0)
      {
        j = i;
        tag[i][0] = 0;
        break;
      }

      if(j >= 0)
      {
#ifdef MPI
        if(ThisTask == 0)
#endif
          printf("%35s\t%10s\n", buf1, buf2);

        switch (id[j])
        {
          case DOUBLE:
          *((double *) addr[j]) = atof(buf2);
          break;
          case STRING:
          strcpy(addr[j], buf2);
          break;
          case INT:
          *((int *) addr[j]) = atoi(buf2);
          break;
        }
      }
      else
      {
          printf("Error in file %s:   Tag '%s' not allowed or multiple defined (entire buffer = %s).\n", fname, buf1, buf);
        errorFlag = 1;
      }
    }
    fclose(fd);

    i = strlen(OutputDir);
    if(i > 0)
      if(OutputDir[i - 1] != '/')
      strcat(OutputDir, "/");
  }
  else
  {
    printf("Parameter file %s not found.\n", fname);
    errorFlag = 1;
  }


  for(i = 0; i < nt; i++)
  {
    if(*tag[i])
    {
      printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
      errorFlag = 1;
    }
  }

    assert(!errorFlag);
    printf("\n");
    
    assert(LastSnapShotNr+1 > 0 && LastSnapShotNr+1 < ABSOLUTEMAXSNAPS);
    MAXSNAPS = LastSnapShotNr + 1;
    
    if(!(NOUT == -1 || (NOUT > 0 && NOUT <= ABSOLUTEMAXSNAPS)))
    printf("NumOutputs must be -1 or between 1 and %i\n", ABSOLUTEMAXSNAPS);
    assert(NOUT == -1 || (NOUT > 0 && NOUT <= ABSOLUTEMAXSNAPS));
    
    // read in the output snapshot list
    if(NOUT == -1)
    {
        NOUT = MAXSNAPS;
        for (i=NOUT-1; i>=0; i--)
        ListOutputSnaps[i] = i;
        printf("all %i snapshots selected for output\n", NOUT);
    }
    else
    {
        printf("%i snapshots selected for output: ", NOUT);
        // reopen the parameter file
        fd = fopen(fname, "r");
        
        done = 0;
        while(!feof(fd) && !done)
        {
            // scan down to find the line with the snapshots
            fscanf(fd, "%s", buf);
            if(strcmp(buf, "->") == 0)
            {
                // read the snapshots into ListOutputSnaps
                for (i=0; i<NOUT; i++)
                {
                    fscanf(fd, "%d", &ListOutputSnaps[i]);
                    printf("%i ", ListOutputSnaps[i]);
                }
                done = 1;
            }
        }
        
        fclose(fd);
        assert(done);
        
        printf("\n");
    }


}
