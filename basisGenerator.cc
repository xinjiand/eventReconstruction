//&&&&&&&&&&&&&
//& C Headers &
//&&&&&&&&&&&&&
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

//&&&&&&&&&&&&&&&&&&
//& Custom Headers &
//&&&&&&&&&&&&&&&&&&
#include "parsing.hh"
#include "detectorParameters.hh"

//&&&&&&&&&&&&&
//& Constants &
//&&&&&&&&&&&&&
#define ENERGY 10000.0

//*************************************************************************************************

//&&&&&&&&
//& Main &
//&&&&&&&&
int main ()
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Declare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  char basisFileName[stringLen_g];    //
  char configFileName[stringLen_g];   //
  int configCheck, xyCellIndex;       //
  FILE * script, * basis, * transIn;  //
  
  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&
  
  //Check the configuration
  printf("This code will use the simple exponential with a 2.2 ns decay constant for the scintillator time response!\n");
  printf("\nIf this is correct enter an a non-zero integer otherwise enter 0\n");
  printf("> ");
  scanf("%d", &configCheck);
  
  if (!configCheck)
  {
    printf("The desired option is not coded. Sorry.\n");
    exit(1);
  }

  //Generate the light transport input file
  transIn = fopen("inputs/trans.dat", "w");
  if (transIn==NULL)
  {
    printf("Could not set the light transport input file!\n");
    exit(1);
  }
  
  printDetectorParameters(transIn);
  fclose(transIn);
  
  //Open the script to generate
  script = fopen("runBasis.sh", "w");
  if (script==NULL)
  {
    printf("Could not open the script! Exiting...\n");
    exit(1);
  }

  //Open the basis data file  
  initString(basisFileName, stringLen_g);
  sprintf(basisFileName, "reconstructionData/reconBasisData_%d_%d%d%d_%d_%d.dat", NCELL,
    LATTICE_FLAG, MIRROR_FLAG, GUIDE_FLAG, TIME_BIN_SIZE, N_TIME_BIN_PRINT);
    
  basis = fopen(basisFileName, "w");
  if (basis==NULL)
  {
    printf("Could not open the basis file! Exiting...\n");
    exit(1);
  }
  
  //Initialize other parameters
  if (LATTICE_FLAG==1) { xyCellIndex = int(NCELL)/2 + 1; }  //Use the center channel in an air-gap lattice
  else { xyCellIndex = 1; }                                 //Use the corner channel in non air-gap lattices
                                                            //to increase the number of side band PMTs
  
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Output the Basis Preamble &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  printDetectorParameters(basis);
  fclose(basis);
  
  
  //&&&&&&&&&&&&&&&&&&&&&&&
  //& Generate the Script &
  //&&&&&&&&&&&&&&&&&&&&&&&
  
  fprintf(script,
    "#This script runs the light transport code and the basis transform code to generated the needed\n");
  fprintf(script, "#waveforms for the reconstruction.\n");
  fprintf(script, "echo \"Start Time:\" ;date\n");

  //Set the configuration file for the light transport
  initString(configFileName, stringLen_g);   
  sprintf(configFileName, "inputs/transConfig.txt");

  fprintf(script, "rm %s\n", configFileName);
  fprintf(script, "echo \"%d\" > %s\n", NCELL, configFileName);
  fprintf(script, "echo \"%d\" >> %s\n", LATTICE_FLAG, configFileName);
  fprintf(script, "echo \"%d\" >> %s\n", MIRROR_FLAG, configFileName);
  fprintf(script, "echo \"%d\" >> %s\n", GUIDE_FLAG, configFileName);
  fprintf(script, "echo \"%d\" >> %s\n", PMT_COUPLING_FLAG, configFileName);
  
  fprintf(script, "echo \"'outputs/trans-pmt.dat'\" >> %s\n", configFileName);
  fprintf(script, "echo \"'outputs/trans-time.dat'\" >> %s\n", configFileName);

  //Call the light transport
  for (int i=1; i<=NCELL; i++)
  {
    fprintf(script, "./trans 0 %lf 0.0 %d %d %d %s\n", ENERGY, xyCellIndex, xyCellIndex, i,
      configFileName);
    fprintf(script, "./addBasis %lf 0.0 %d %d %d\n", ENERGY, xyCellIndex, xyCellIndex, i);
    if (i==1) { fprintf(script, "echo \"Finished with the first cell\"\n"); }
    else { fprintf(script, "echo \"Finished with %d cells\"\n", i); }
  }
	fprintf(script, "echo \"End Time:\" ;date\n");
	
  chmod("runBasis.sh", 0777); //Make the script executable

  //&&&&&&&
  //& End &
  //&&&&&&&
  
  fclose(script);
  printf("Now run the script runBasis.sh\n");
  
  return 0;
}

//*************************************************************************************************
