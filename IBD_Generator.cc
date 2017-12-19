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

//&&&&&&&&&&&&&&&&&&&&&&&&&&
//& GNU Scientific Library &
//&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <gsl/gsl_rng.h>

//&&&&&&&&&&&&&&&&&&
//& Custom Headers &
//&&&&&&&&&&&&&&&&&&
#include "parsing.hh"
#include "randomNumberGenerator.hh"
#include "reconLib.hh"

//&&&&&&&&&&&&&
//& Constants &
//&&&&&&&&&&&&&
#define VERBOSE_LEVEL 0

//&&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Declarations &
//&&&&&&&&&&&&&&&&&&&&&&&&&

void determineEventStartPoint(int distribFlag, vector_t * delta, int genFlag, gsl_rng * rGen);
void printDeposit(deposit_t dep, int depNumber, FILE * f, int particleID);
void deleteEvent(double energy[][NCELL][NCELL], double timeFirst[][NCELL][NCELL],
  double timeCent[][NCELL][NCELL]);
int addDeposit(deposit_t dep, double energy[][NCELL][NCELL], double timeFirst[][NCELL][NCELL], 
  double timeCent[][NCELL][NCELL]);

//*************************************************************************************************

//&&&&&&&&
//& Main &
//&&&&&&&&
int main ()
{
  //TO DO: Implement script splitting

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Declare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  //Varibles for program configuration
  char positronFileName[stringLen_g];     //TO DO: Go through and comment these
  char neutronFileName[stringLen_g];      //
  char scriptDirName[stringLen_g];        //
  char scriptName[stringLen_g];           //
  char configFileName[stringLen_g];       //
  char positronScriptName[stringLen_g];   //
  char neutronScriptName[stringLen_g];    //
  char eventDirName[stringLen_g];         //
  char eventFileName[stringLen_g];        //
  char string0[stringLen_g];              //
  char string1[stringLen_g];              //
  char string2[stringLen_g];              //
  char string3[stringLen_g];              //
  char readPositron[stringLen_g];         //
  char readNeutron[stringLen_g];          //
  int genFlag, seed, seedFlag;            //
  int nCells, latticeFlag, mirrorFlag;    //
  int guideFlag, distribFlag;             //
  int quenchingFlag, scinTimeResponse;    //
  int configCheck, nEnergy;               //

  //Variables for reading the geant data
  int check, nLines, positronEventNum;    //
  int neutronEventNum;                    //
  int nEntries, eventCnt=0, geantSeed;    //
  int startParticleID;                    //
  int nDep, nPositronDep, nNeutronDep;    //
  int firstAlpha;                         //
  double initialEnergyPositron;           //
  double initialEnergyNeutron;            //
  double startThetaPositron;              //
  double startPhiPositron;                //
  double startThetaNeutron;               //
  double startPhiNeutron;                 //
  double startTimePositron;               //
  double startTimeNeutron;                //
  double energyIBD, deltaT, E;            //
  vector_t startPosPositron;              //
  vector_t startPosNeutron;               //
  deposit_t * depPositron, * depNeutron;  // 
  deposit_t * tempList;                   //
  
  //Variables for generating events
  double alpha, gamma;                            //
  double cAlpha, sAlpha;                          //
  double cBeta, sBeta;                            //
  double cGamma, sGamma;                          //
  double rotMatrix[3][3];                         //
  double u, v, w, mag;                            //
  vector_t * delta = new vector_t;                //
  vector_t positronInitDir, neutronInitDir;       //
  double positronEnergy[NCELL][NCELL][NCELL];     //
  double positronTimeFirst[NCELL][NCELL][NCELL];  //
  double positronTimeCent[NCELL][NCELL][NCELL];   //
  double neutronEnergy[NCELL][NCELL][NCELL];      //
  double neutronTimeFirst[NCELL][NCELL][NCELL];   //
  double neutronTimeCent[NCELL][NCELL][NCELL];    //
  int nDepOutside, nInAfterOut;                   //
  int nPositronCells, nNeutronCells, depCnt;      //
  gsl_rng *  rGen;                                //
  const gsl_rng_type * T;                         //
  
  //Files
  DIR * dir;                                    //
  FILE * config, * list, * positron, * neutron; //
  FILE * script, * out, * event, * transIn;     //
  FILE * positronScript, * neutronScript;       //
   

  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&

  //Clear the cell lists
  deleteEvent(positronEnergy, positronTimeFirst, positronTimeCent);
  deleteEvent(neutronEnergy, neutronTimeFirst, neutronTimeCent);
  
  //Read the program control variables
  config = fopen("inputs/readConfig.txt", "r");
  if (config==NULL)
  {
    printf("Could not open the file readConfig.txt! Exiting...\n");
    exit(1);
  }

  if (getNumEntries(config)!=9)
  {
    printf("The file inputs/readConfig.txt is improperly formatted! Exiting...\n");
    exit(1);
  }
  
  for (int i=0; i<9; i++)
  {
    initString(string0, stringLen_g);
    initString(string1, stringLen_g);
    getEntry(config, string0, stringLen_g);  
    
    if (i==0)       { sscanf(string0, "%s %d", string1, &genFlag);          }
    else if (i==1)  { sscanf(string0, "%s %d", string1, &seedFlag);         }
    else if (i==2)  { sscanf(string0, "%s %d", string1, &nCells);           }
    else if (i==3)  { sscanf(string0, "%s %d", string1, &latticeFlag);      }
    else if (i==4)  { sscanf(string0, "%s %d", string1, &mirrorFlag);       }
    else if (i==5)  { sscanf(string0, "%s %d", string1, &guideFlag);        }
    else if (i==6)  { sscanf(string0, "%s %d", string1, &distribFlag);      }
    else if (i==7)  { sscanf(string0, "%s %d", string1, &quenchingFlag);    }
    else            { sscanf(string0, "%s %d", string1, &scinTimeResponse); }
  }
  
  if (distribFlag==2)
  {
    printf("IBD events are internal to the detector! Edit inputs/readConfig.txt and re-run.\n");
    exit(1);
  }
  
  fclose(config);

  //Initialize the random number generator
  if (genFlag!=0)
  {
    if (genFlag==1)       T = gsl_rng_default;
    else if (genFlag==2)  T = gsl_rng_ranlxs0;
    else if (genFlag==3)  T = gsl_rng_ranlxs1;
    else if (genFlag==4)  T = gsl_rng_ranlxs2;
    else if (genFlag==5)  T = gsl_rng_ranlxd1;
    else if (genFlag==6)  T = gsl_rng_ranlxd2;
    else if (genFlag==7)  T = gsl_rng_ranlux;
    else if (genFlag==8)  T = gsl_rng_ranlux389;
    else if (genFlag==9)  T = gsl_rng_cmrg;
    else if (genFlag==10) T = gsl_rng_mrg;
    else if (genFlag==11) T = gsl_rng_taus;
    else if (genFlag==12) T = gsl_rng_taus2;
    else if (genFlag==13) T = gsl_rng_gfsr4;
    else 
    {
      printf("You entered an illegal value!\n");
      printf("Generator set to GSL default.\n");
      T = gsl_rng_default;
    }

    rGen = gsl_rng_alloc(T);
  }  
  
  if (seedFlag==0) {seed=1013;}
  else {seed=time(NULL);}
  rngSeed(genFlag, seed, rGen);

  //Open the output file with the script names written to it
  out = fopen("./inputs/scriptGenerator.in", "w");
  if (out==NULL)
  {
    printf("Could not open the file inputs/scriptGenerator.in! Exiting...\n");
    exit(1);
  }


  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Confirm the Configuration &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  if (NCELL!=nCells || LATTICE_FLAG!=latticeFlag || MIRROR_FLAG!=mirrorFlag ||
    GUIDE_FLAG!=guideFlag)
  {
    printf("There is a mismatch between the macros defined in detectorParameters.hh and in readConfig.txt!\n");
    printf("Adjust this before running. Exiting...\n");
    exit(1);
  }
  
  //Tell user the configuration
  printf("Dear user this is the configuration that you inputted:\n");
  printf("  The seed is %d\n", int(seed));

  //Detector size
  printf("  The detector size will be %dx%dx%d cells.\n", NCELL, NCELL, NCELL);

  //Event distribution
  if (distribFlag==0)
  {
    printf("  The events will be generated in the CENTER of the detector.\n");
  } 
  else if (distribFlag==1)
  {
    printf("  The events will be generated UNIFORM in the detector.\n");
  }
  else if (distribFlag==3)
  {
    printf("  The GEANT4 deposits will be used.\n");
    printf("  Ensure that the Geant4 and light transport geometries are the same!\n");     
  }
  else
  {
    printf("  The option %d does not exist for distribFlag! Exiting...\n", distribFlag);   
    exit(1);     
  }

  //Mirrors
  if (mirrorFlag) { printf("  Mirrors are ON.\n"); }
  else { printf("  Mirrors are OFF.\n"); }

  //Guides
  if (guideFlag==1) { printf("  Light guides are ON with %d channel into 1 guide.\n", guideFlag); }
  else if (guideFlag!=0)
  {
    printf("  Light guides are ON with %d channels into 1 guides.\n", guideFlag);
    printf("\nWARNING! The reconstruction will not work with these parameters!\n");
  }
  else
  {
    printf("  Light guides are OFF.\n");
    printf("\nWARNING! The reconstruction will not work with these parameters!\n");
    guideFlag=0;  //To ensure a proper call to the light transport function
  }

  //Quenching
  if (quenchingFlag==0) { printf("  Assume legacy quenching method\n"); }
  else { printf("  Assume all deposits are given in MeVee\n"); }
  
  //Scintillator response
  if (scinTimeResponse==0)
  {
    printf("  Scintillator has a time response of a pure exponential with a 2.2 ns decay constant (EJ-254 datasheet)\n");
  }
  else
  {
    //TO DO: Give the correct reference when the paper is publised
    printf("  Scintillator has the time responses from M. J. I. Balmer, et al., NIM, A788, 146 (2015)\n");  
  }
  
  //Reconstruction program
  if (RECON_PROGRAM==0)
  {
    printf("\nThe program chargeReconstruction.f will be used. Note this is a charge only reconstruction\n");
  }
  else
  {
    printf("\nThe program reconstruct.cc will be used. Set options in inputs/reconConfig.txt.\n");  
  }

  //confirm configuration
  printf("\nIf this configuration is correct enter an a non-zero integer otherwise enter 0\n");
  printf("> ");
  scanf("%d", &configCheck);
  
  if (!configCheck)
  {
    printf("Edit the configuration in inputs/readConfig.txt. Exiting...\n");
    fclose(out);
    return 0;
  }
  

  //&&&&&&&&&&&&&&&&&&&&&&&
  //& Read the Data Files &
  //&&&&&&&&&&&&&&&&&&&&&&&

  list = fopen("./inputs/IBD.in", "r");
  if (list==NULL)
  {
    printf("Cannot open the file! Exiting...\n");
    exit(1); 
  }
  nEntries = getNumEntries(list);
  
  if (nEntries%2!=0)
  {
    printf("Please list a positron file then a neutron file. Exiting...\n");
    exit(1);
  }
   

  //&&&&&&&&&&&&&&&&&
  //& Read the Data &
  //&&&&&&&&&&&&&&&&&

  for (int i=0; i<nEntries/2; i++)
  {
    //Get the positron file name and open
    initString(positronFileName, stringLen_g);
    getEntry(list, positronFileName, stringLen_g);

    positron = fopen(positronFileName, "r");
    if (positron==NULL)
    {
      printf("Could not open %s! Skipping to the next pair in the list.\n", positronFileName);
      getEntry(list, neutronFileName, stringLen_g);
      continue;
    }

    //Get the neutron file name and open
    initString(neutronFileName, stringLen_g);
    getEntry(list, neutronFileName, stringLen_g);

    neutron = fopen(neutronFileName, "r");
    if (neutron==NULL)
    {
      printf("Could not open %s! Skipping to the next pair in the list.\n", neutronFileName);
      continue;    
    }
    
    //Check script directory
    copyString(string0, positronFileName, strlen(positronFileName));
    removePath(string0, 0); //Get the energy of the incident neutrino
    
    for (int si=0; si<int(strlen(string0)); si++)
    {
      if (string0[si]=='_') string0[si] = ' ';
    }
    
    initString(string1, stringLen_g);
    initString(string2, stringLen_g);    
    initString(string3, stringLen_g);    
    sscanf(string0, "%s %s %s", string1, string2, string3);
    
    initString(string0, stringLen_g);
    for (int si=0; si<int(strlen(string2))-3; si++)
    {
      string0[si] = string2[si];
    }
    
    sscanf(string0, "%d", &nEnergy);
    initString(scriptDirName, stringLen_g);
    sprintf(scriptDirName, "simulationScripts/IBD_%05dkeV", nEnergy+1800);

    dir = opendir(scriptDirName);
    if (dir) { closedir(dir); } //Directory exists -> close and move on
    else if (errno==ENOENT)     //Directory does not exist -> create
    {
      check = mkdir(scriptDirName, 0777);
      if (check)
      {
        printf("Error in opening script directory! Skipping to the next file in the list.\n");
        fclose(positron);
        fclose(neutron);        
        continue;
      }
    }
    else                        //opendir() failed for some other reason
    {
      printf("Error in opening script directory! Skipping to the next file in the list.\n");
      fclose(positron);
      fclose(neutron);
      continue;
    }
    strcat(scriptDirName, "/");

    //Open the bash script and name the script similar to the data file name
    copyString(scriptName, scriptDirName, strlen(scriptDirName));
    initString(string0, stringLen_g);
    sprintf(string0, "IBD_%05dkeV_000000000_Run.sh", nEnergy+1800);
    strcat(scriptName, string0);  
        
    fprintf(out, "%s\n", scriptName);
    script = fopen(scriptName, "w");
    if (script==NULL)
    {
      printf("Could not open %s! Skipping to the next file in the list\n", scriptName);
      fclose(positron);
      fclose(neutron);
      continue;
    }

    chmod(scriptName, 0777); //Make the script executable

    //Generate the script
    
    //Start with clean output files
    fprintf(script, "rm outputs/lightTransPositron.out\n");
    fprintf(script, "rm outputs/lightTransNeutron.out\n");
    if (RECON_PROGRAM==0) fprintf(script, "rm outputs/reconstruction.out\n"); //If using Bruce's
                                                                              //reconstrution
                                                                              
    //Generate the light transport configuration file    
    initString(configFileName, stringLen_g);   
    sprintf(configFileName, "inputs/transConfigPositron.txt");
    fprintf(script, "rm %s\n", configFileName);
    
    fprintf(script, "echo \"%d\" > %s\n", NCELL, configFileName);
    fprintf(script, "echo \"%d\" >> %s\n", LATTICE_FLAG, configFileName);
    fprintf(script, "echo \"%d\" >> %s\n", MIRROR_FLAG, configFileName);
    fprintf(script, "echo \"%d\" >> %s\n", GUIDE_FLAG, configFileName);
    
    fprintf(script, "echo \"'outputs/trans-pmt-positron.dat'\" >> %s\n", configFileName);
    fprintf(script, "echo \"'outputs/trans-time-positron.dat'\" >> %s\n", configFileName);
    
    initString(configFileName, stringLen_g);   
    sprintf(configFileName, "inputs/transConfigNeutron.txt");
    fprintf(script, "rm %s\n", configFileName);
    
    fprintf(script, "echo \"%d\" > %s\n", NCELL, configFileName);
    fprintf(script, "echo \"%d\" >> %s\n", LATTICE_FLAG, configFileName);
    fprintf(script, "echo \"%d\" >> %s\n", MIRROR_FLAG, configFileName);
    fprintf(script, "echo \"%d\" >> %s\n", GUIDE_FLAG, configFileName);
    
    fprintf(script, "echo \"'outputs/trans-pmt-neutron.dat'\" >> %s\n", configFileName);
    fprintf(script, "echo \"'outputs/trans-time-neutron.dat'\" >> %s\n", configFileName);    
    
    //Generate the scripts for the positron and the neutron seperatly
    copyString(positronScriptName, scriptDirName, strlen(scriptDirName));
    copyString(neutronScriptName, scriptDirName, strlen(scriptDirName));

    initString(string0, stringLen_g);
    sprintf(string0, "IBD_%05dkeV_positron_000000000_Run.sh", nEnergy+1800);
    initString(string1, stringLen_g);  
    sprintf(string1, "IBD_%05dkeV_neutron_000000000_Run.sh", nEnergy+1800);
    
    strcat(positronScriptName, string0);      
    strcat(neutronScriptName, string1);
    
    positronScript = fopen(positronScriptName, "w");
    if (positronScript==NULL)
    {
      printf("Could not open %s! Skipping to the next file in the list\n", positronScriptName);
      fclose(positron);
      fclose(neutron);
      fclose(script);
      continue;    
    }
    
    neutronScript = fopen(neutronScriptName, "w");
    if (neutronScript==NULL)
    {
      printf("Could not open %s! Skipping to the next file in the list\n", neutronScriptName);
      fclose(positron);
      fclose(neutron);
      fclose(script);
      fclose(positronScript);
      continue;    
    }

    chmod(positronScriptName, 0777);
    chmod(neutronScriptName, 0777);
        
    fprintf(script, "./%s\n", positronScriptName);
    fprintf(script, "./%s\n", neutronScriptName);
    fprintf(script, "./combineIBD_Events 0\n");
    fprintf(script, "echo \"%s processed!\"\n", scriptName);
    
    fclose(script);

    //Open the event file directory
    initString(eventDirName, stringLen_g);
    sprintf(eventDirName, "eventFiles/IBD_%05dkeV", nEnergy+1800);

    dir = opendir(eventDirName);
    if (dir) { closedir(dir); } //Directory exists -> close and move on
    else if (errno==ENOENT)     //Directory does not exist -> create
    {
      check = mkdir(eventDirName, 0777);
      if (check)
      {
        printf("Error in opening script directory! Skipping to the next file in the list.\n");
        fclose(positron);
        fclose(neutron);
        fclose(positronScript);
        fclose(neutronScript);
        continue;
      }
    }
    else                        //opendir() failed for some other reason
    {
      printf("Error in opening script directory! Skipping to the next file in the list.\n");
      fclose(positron);
      fclose(neutron);
      fclose(positronScript);
      fclose(neutronScript);
      continue;
    }
    strcat(eventDirName, "/");

    //Open the event file
    copyString(eventFileName, eventDirName, strlen(eventDirName));
    initString(string0, stringLen_g);
    sprintf(string0, "IBD_%05dkeV_000000000_Event.dat", nEnergy+1800); //For now assume just one file
    strcat(eventFileName, string0); 

    event = fopen(eventFileName, "w");
    if (event==NULL)
    {
      printf("Could not open %s! Skipping to the next file in the list\n", string0);
      fclose(positron);
      fclose(neutron);
      fclose(positronScript);
      fclose(neutronScript);
      continue;
    }
    if (VERBOSE_LEVEL==0) printf("%s has ", eventFileName);

    //Check if this is a neutron data file and if so it is for boron or lithium loading
    int boronOrLithium=0;
    if (strstr(neutronFileName, "boron")) { boronOrLithium = 1; }
    else if (strstr(neutronFileName, "lithium")) { boronOrLithium = 2; }
 
    //Loop over the data and generate the script
    initString(readPositron, stringLen_g);
    initString(readNeutron, stringLen_g);    
    while (fgets(readPositron, stringLen_g, positron)!=NULL && 
      fgets(readNeutron, stringLen_g, neutron)!=NULL)
    {
      //Check the event headers
      if (strstr(readPositron, "#Event")==NULL || strstr(readNeutron, "#Event")==NULL)
      {
        printf("Error at event header!\n");
        printf("  Positron: %s", readPositron);
        printf("  Neutron: %s", readNeutron);
        continue;
      }
            
      //Read the positron event header
      initString(readPositron, stringLen_g);
      fgets(readPositron, stringLen_g, positron);
      sscanf(readPositron, "%d", &nLines);
      for (int j=0; j<nLines; j++)
      {
        initString(readPositron, stringLen_g);
        fgets(readPositron, stringLen_g, positron);
        if (j==0) { sscanf(readPositron, "%d", &positronEventNum); }
        else if (j==1) { sscanf(readPositron, "%d", &geantSeed); }
      }
      
      //Read the neutron event header
      initString(readNeutron, stringLen_g);
      fgets(readNeutron, stringLen_g, neutron);
      sscanf(readNeutron, "%d", &nLines);
      for (int j=0; j<nLines; j++)
      {
        initString(readNeutron, stringLen_g);
        fgets(readNeutron, stringLen_g, neutron);
        if (j==0) { sscanf(readNeutron, "%d", &neutronEventNum); }
        //****Using the positron "geantSeed", but really it doesn't matter this this has not been
        //****implemented yet. 2015-01-31
        //TO DO: Place both the positron and the neutron seed in the event file
      }      
      
      //Write the IBD event header
      if (positronEventNum!=neutronEventNum)
      {
        printf("The positron and neutron event numbers are not aligned if you care about that.\n");
      }
      
      initString(string0, stringLen_g);
      sprintf(string0, "#Event\n2\n  %d\n  %d\n", positronEventNum, geantSeed); //Just use the
      fputs(string0, event);                                                    //positron event
                                                                                //number because who 
                                                                                //cares
      //Check the Geant4 properties headers
      initString(readPositron, stringLen_g);
      initString(readNeutron, stringLen_g);
      fgets(readPositron, stringLen_g, positron);
      fgets(readNeutron, stringLen_g, neutron);
      if (strstr(readPositron, "#Geant4_properties")==NULL || 
        strstr(readNeutron, "#Geant4_properties")==NULL)
      {
        printf("Error at Geant4 properties header!\n");
        printf("  Positron: %s", readPositron);
        printf("  Neutron: %s", readNeutron);
        continue;
      } 
 
      //Read the positron Geant4 properties header
      initString(readPositron, stringLen_g);
      fgets(readPositron, stringLen_g, positron);
      sscanf(readPositron, "%d", &nLines);
      for (int j=0; j<nLines; j++)
      {
        initString(readPositron, stringLen_g);
        fgets(readPositron, stringLen_g, positron);
        if (j==0) { sscanf(readPositron, "%d", &startParticleID); }
        else if (j==1) { sscanf(readPositron, "%lf", &initialEnergyPositron); }
        else if (j==2) { sscanf(readPositron, "%lf %lf %lf", &startPosPositron.x,
          &startPosPositron.y, &startPosPositron.z); }
        else if (j==3) { sscanf(readPositron, "%lf %lf", &startThetaPositron, &startPhiPositron); }
        else if (j==4) { sscanf(readPositron, "%lf", &startTimePositron); }
      }
 
      //Read the neutron Geant4 properties header
      initString(readNeutron, stringLen_g);
      fgets(readNeutron, stringLen_g, neutron);
      sscanf(readNeutron, "%d", &nLines);
      for (int j=0; j<nLines; j++)
      {
        initString(readNeutron, stringLen_g);
        fgets(readNeutron, stringLen_g, neutron);
        if (j==0) { sscanf(readNeutron, "%d", &startParticleID); }
        else if (j==1) { sscanf(readNeutron, "%lf", &initialEnergyNeutron); }
        else if (j==2) { sscanf(readNeutron, "%lf %lf %lf", &startPosNeutron.x, &startPosNeutron.y, 
          &startPosNeutron.z); }
        else if (j==3) { sscanf(readNeutron, "%lf %lf", &startThetaNeutron, &startPhiNeutron); }
        else if (j==4) { sscanf(readNeutron, "%lf", &startTimeNeutron); }
      }
      
      //Write the IBD Geant4 properties header
      initString(string0, stringLen_g);
      sprintf(string0, "#Geant4_properties\n5\n  100\n");
      fputs(string0, event);
      
      energyIBD = initialEnergyPositron + 1800.0;
      initString(string0, stringLen_g);
      sprintf(string0, "  %lf\n", energyIBD);
      fputs(string0, event);
      
      if ( startPosPositron.x!=startPosNeutron.x || startPosPositron.y!=startPosNeutron.y || 
        startPosPositron.z!=startPosNeutron.z )
      {
        printf("The starting positions do not match. Using the positron initial position\n");
      }
      initString(string0, stringLen_g);
      sprintf(string0, "  %lf %lf %lf\n", startPosPositron.x, startPosPositron.y,
        startPosPositron.z);
      fputs(string0, event);
      
      if ( startThetaPositron!=startThetaNeutron || startPhiPositron!=startPhiNeutron )
      {
        printf("The starting directions are not the same. Using the positron initial direction\n");
      }
      initString(string0, stringLen_g);
      sprintf(string0, "  %lf %lf\n", startThetaPositron, startPhiPositron);
      fputs(string0, event);
      
      if (startTimePositron!=startTimeNeutron)
      {
        printf("The start times are not the same. Using the positron start time and offsetting");
        printf(" the neutron times to match\n");
        deltaT = startTimePositron - startTimeNeutron;
      }
      else { deltaT=0.0; }
      initString(string0, stringLen_g);
      sprintf(string0, "  %lf\n", startTimePositron);
      fputs(string0, event);
      
      //Check event deposits header
      initString(readPositron, stringLen_g);
      initString(readNeutron, stringLen_g);      
      fgets(readPositron, stringLen_g, positron);
      fgets(readNeutron, stringLen_g, neutron);
      if (strstr(readPositron, "#Energy_deposits")==NULL || 
        strstr(readNeutron, "#Energy_deposits")==NULL)
      {
        printf("Error at energy properties header!\n");
        printf("  Positron: %s", readPositron);
        printf("  Neutron: %s", readNeutron);
        continue;
      }       
      
      //Write the energy deposits header
      initString(readPositron, stringLen_g);
      initString(readNeutron, stringLen_g); 
      fgets(readPositron, stringLen_g, positron);      
      fgets(readNeutron, stringLen_g, neutron);      

      sscanf(readPositron, "%d", &nPositronDep);
      sscanf(readNeutron, "%d", &nNeutronDep);
      nDep = nPositronDep + nNeutronDep;

      initString(string0, stringLen_g);
      sprintf(string0, "#Energy_deposits\n%d %d %d\n", nDep, nPositronDep, nNeutronDep);
      fputs(string0, event);
            
      //Allocate the deposit arrays
      depPositron = new deposit_t [nPositronDep];
      for (int j=0; j<nPositronDep; j++)
      {
        depPositron[j].id=0;
        depPositron[j].x=0.0;
        depPositron[j].y=0.0;
        depPositron[j].z=0.0;
        depPositron[j].t=0.0;
        depPositron[j].energy=0.0;
      }
      
      depNeutron = new deposit_t [nNeutronDep];
      for (int j=0; j<nNeutronDep; j++)
      {
        depNeutron[j].id=0;
        depNeutron[j].x=0.0;
        depNeutron[j].y=0.0;
        depNeutron[j].z=0.0;
        depNeutron[j].t=0.0;
        depNeutron[j].energy=0.0;
      }      
      
      //Read and modify the positron energy deposits
      for (int j=0; j<nPositronDep; j++)
      {
        initString(readPositron, stringLen_g);
        fgets(readPositron, stringLen_g, positron);
        fputs(readPositron, event);
        
        sscanf(readPositron, "%d %lf %lf %lf %lf %lf %lf", 
          &depPositron[j].id, &depPositron[j].x, &depPositron[j].y, &depPositron[j].z, &depPositron[j].t, 
            &E, &depPositron[j].energy);
            
        depPositron[j].x = depPositron[j].x/10.0; //Convert from mm to cm
        depPositron[j].y = depPositron[j].y/10.0;
        depPositron[j].z = depPositron[j].z/10.0;
      }
      
      //Set the deposit id's so that the light transport generates the scintillation light with
      //the correct time dependence
      if (scinTimeResponse==0)  //Pure exponential with a 2.2 ns decay time
      {
        for (int depNum=0; depNum<nPositronDep; depNum++) { depPositron[depNum].id = 0; }
      }
      else                      //Uses the gamma, fast neutron, or thermal neutron time response
      {                         //from "Comparative analysis of pulse shape discrimination..."
        for (int depNum=0; depNum<nPositronDep; depNum++)
        {
          if      (depPositron[depNum].id<=3)  { depPositron[depNum].id = 1; }        //Electrons/gammas
          else if (depPositron[depNum].id==4)  { depPositron[depNum].id = 2; }        //Proton scatters
          else if (depPositron[depNum].id==8)  { depPositron[depNum].id = 3; }        //Neutron capture
          else                                 { depPositron[depNum].energy = 0.0; }  //Axe deposits from other than
        }                                                                           //e-, e+, gamma, proton, or alpha
      }
 
      //Read and modify the neutron energy deposits
      firstAlpha=1;
      for (int j=0; j<nNeutronDep; j++)
      {
        initString(readNeutron, stringLen_g);
        fgets(readNeutron, stringLen_g, neutron);
        sscanf(readNeutron, "%d %lf %lf %lf %lf %lf %lf", 
          &depNeutron[j].id, &depNeutron[j].x, &depNeutron[j].y, &depNeutron[j].z, &depNeutron[j].t, &E, 
            &depNeutron[j].energy);
        depNeutron[j].t += deltaT;

        if (quenchingFlag==0) //Adjust energy to reflect quenching in the scintillator  
        {        
          if (depNeutron[j].id>3  || depNeutron[j].id==0) //Zero energy deposits not an electron,
          {                                               //positron, or gamma.
            depNeutron[j].energy = 0.0;
          }
                                                          
          if (depNeutron[j].id==8 && firstAlpha==1) //If the particle is an alpha set it to the
          {                                         //electron equivalent 
            if      (boronOrLithium==1) { depNeutron[j].energy = 0.07; }
            else if (boronOrLithium==2) { depNeutron[j].energy = 0.40; }
            firstAlpha=0;
          }
        }
 
        initString(string0, stringLen_g);
        sprintf(string0, "  %d %lf %lf %lf %lf %lf %lf\n", depNeutron[j].id, depNeutron[j].x,
          depNeutron[j].y, depNeutron[j].z, depNeutron[j].t, E, depNeutron[j].energy);
        fputs(string0, event);
        
        depNeutron[j].x = depNeutron[j].x/10.0; //Conver from mm to cm
        depNeutron[j].y = depNeutron[j].y/10.0;
        depNeutron[j].z = depNeutron[j].z/10.0;
      }

      //Set the deposit id's so that the light transport generates the scintillation light with
      //the correct time dependence
      if (scinTimeResponse==0)  //Pure exponential with a 2.2 ns decay time
      {
        for (int depNum=0; depNum<nNeutronDep; depNum++) { depNeutron[depNum].id = 0; }
      }
      else                      //Uses the gamma, fast neutron, or thermal neutron time response
      {                         //from "Comparative analysis of pulse shape discrimination..."
        for (int depNum=0; depNum<nNeutronDep; depNum++)
        {
          if      (depNeutron[depNum].id<=3)  { depNeutron[depNum].id = 1; }        //Electrons/gammas
          else if (depNeutron[depNum].id==4)  { depNeutron[depNum].id = 2; }        //Proton scatters
          else if (depNeutron[depNum].id==8)  { depNeutron[depNum].id = 3; }        //Neutron capture
          else                                { depNeutron[depNum].energy = 0.0; }  //Axe deposits from other than
        }                                                                           //e-, e+, gamma, proton, or alpha
      }
      
      //Remove the zero energy
      int aboveZeroCnt=0;
      for (int j=0; j<nNeutronDep; j++)
      {
        if (depNeutron[j].energy>0.0) aboveZeroCnt++;
      }

      tempList = new deposit_t [aboveZeroCnt];
      aboveZeroCnt=0;
      for (int j=0; j<nNeutronDep; j++)
      {
        if (depNeutron[j].energy>0.0)
        {
          tempList[aboveZeroCnt].id = depNeutron[j].id;
          tempList[aboveZeroCnt].x = depNeutron[j].x;
          tempList[aboveZeroCnt].y = depNeutron[j].y;
          tempList[aboveZeroCnt].z = depNeutron[j].z;
          tempList[aboveZeroCnt].t = depNeutron[j].t;
          tempList[aboveZeroCnt].energy = depNeutron[j].energy;
          
          aboveZeroCnt++;
        }
      }
      
      delete [] depNeutron;
      nNeutronDep = aboveZeroCnt;
      depNeutron = new deposit_t [nNeutronDep];
      for (int j=0; j<nNeutronDep; j++)
      {
        depNeutron[j].id = tempList[j].id;
        depNeutron[j].x = tempList[j].x;
        depNeutron[j].y = tempList[j].y;
        depNeutron[j].z = tempList[j].z;
        depNeutron[j].t = tempList[j].t;
        depNeutron[j].energy = tempList[j].energy;      
      }
      
      delete [] tempList;

      for (int j=1; j<nNeutronDep; j++) //Offset the neutron times
      {
        depNeutron[j].t = depNeutron[j].t - depNeutron[0].t;
      }
      depNeutron[0].t = 0.0;

      //Set the positron event
      
      //Set the offset
      delta->x=0.0; delta->y=0.0; delta->z=0.0;
      determineEventStartPoint(distribFlag, delta, genFlag, rGen);
            
      //Define the random Euler angles
      
      //TO DO: Implement the differential cross section for emission of the positron
      
      alpha = 2.0*PI*rng(genFlag, rGen);    //Assume that the positron direction is isotropic with
      cBeta = 2.0*rng(genFlag, rGen) - 1.0; //respect to the incident neutrino direction. This is 
      gamma = 2.0*PI*rng(genFlag, rGen);    //not strictly true (see P. Vogel and J. F. Beacom,
                                            //Phys. Rev. D, 60, 053003 (1999)), but it is a good
                                            //enough approximation for now. 

      //Set the starting positrion in light transport properties
      fputs("#Light_transport_properties\n3\n", event);      
      
      initString(string0, stringLen_g);
      sprintf(string0, "  %lf %lf %lf\n", delta->x+startPosPositron.x, delta->y+startPosPositron.y, 
        delta->z+startPosPositron.z);
      fputs(string0, event);

      //Define the rotation matrix      
      cAlpha = cos(alpha);
      sAlpha = sin(alpha);
      sBeta = sqrt(1.0 - cBeta*cBeta);
      cGamma = cos(gamma);
      sGamma = sin(gamma);

      rotMatrix[0][0] = cGamma*cBeta*cAlpha - sGamma*sAlpha;
      rotMatrix[0][1] = cGamma*cBeta*sAlpha + sGamma*cAlpha;
      rotMatrix[0][2] = -cGamma*sBeta;
      rotMatrix[1][0] = -sGamma*cBeta*cAlpha - cGamma*sAlpha;
      rotMatrix[1][1] = -sGamma*cBeta*sAlpha + cGamma*cAlpha;
      rotMatrix[1][2] = sGamma*sBeta;
      rotMatrix[2][0] = sBeta*cAlpha;
      rotMatrix[2][1] = sBeta*sAlpha;
      rotMatrix[2][2] = cBeta;            
            
      //Apply the transformation and add to the detector
      nDepOutside=0;
      nInAfterOut=0;
      depCnt=0;
      for (int j=0; j<nPositronDep; j++)
      {
        u = rotMatrix[0][0]*depPositron[j].x + rotMatrix[0][1]*depPositron[j].y + 
          rotMatrix[0][2]*depPositron[j].z;
        v = rotMatrix[1][0]*depPositron[j].x + rotMatrix[1][1]*depPositron[j].y + 
          rotMatrix[1][2]*depPositron[j].z;
        w = rotMatrix[2][0]*depPositron[j].x + rotMatrix[2][1]*depPositron[j].y + 
          rotMatrix[2][2]*depPositron[j].z;             

        u = u + delta->x;  v = v + delta->y;  w = w + delta->z;
        depPositron[j].x=u;  depPositron[j].y=v;  depPositron[j].z=w;
        
        check = addDeposit(depPositron[j], positronEnergy, positronTimeFirst, positronTimeCent);
        
        if (check==0) { nDepOutside++; }
        else
        {
          if (depPositron[j].energy>E_MIN)
          {
            if (depCnt==0)
            {
              fprintf(positronScript, "echo \"#Event %d\" >> outputs/lightTransPositron.out\n",
                positronEventNum);
            } 

            printDeposit(depPositron[j], depCnt, positronScript, 3);
            depCnt++;
            if (nDepOutside!=0) nInAfterOut++;     
                            
          } //End check if the deposit is above threshold
        } //End check if inside the detector
      } //End loop over positron deposits
      
      if (VERBOSE_LEVEL==1 && nDepOutside!=0 && nInAfterOut!=0)
      {
        printf("Positron Event %d:\n", positronEventNum);
        printf("  nDepOutside=%d\n", nDepOutside);
        printf("  nInAfterOut=%d\n", nInAfterOut);
      }
      
      if (depCnt!=0)
      {
        fprintf(positronScript, "NUM_LINES=$(cat outputs/trans-pmt-positron.dat | wc -l )\n");
        fprintf(positronScript, "echo $NUM_LINES >> outputs/lightTransPositron.out\n");
        fprintf(positronScript, "cat outputs/trans-pmt-positron.dat >> outputs/lightTransPositron.out\n");
        
        fprintf(positronScript, "NUM_LINES=$(cat outputs/trans-time-positron.dat | wc -l )\n");
        fprintf(positronScript, "echo $NUM_LINES >> outputs/lightTransPositron.out\n");          
        fprintf(positronScript, "cat outputs/trans-time-positron.dat >> outputs/lightTransPositron.out\n");     
      }
      else
      {
        fprintf(positronScript, "echo \"#Event %d\" >> outputs/lightTransPositron.out\n",
          positronEventNum);
        fprintf(positronScript, "echo \"0\" >> outputs/lightTransPositron.out\n");    
        fprintf(positronScript, "echo \"0\" >> outputs/lightTransPositron.out\n");
      }
      
      if (RECON_PROGRAM==0) //If using Bruce's reconstruction then call here
      {
        fprintf(positronScript, "./chargeReconstruction %d %d %d %d\n", NCELL, LATTICE_FLAG,
          MIRROR_FLAG, positronEventNum);
      }
          
      //Set the neutron event     
      
      //Get the neutron's direction. Assume that the neutrino is incident from the z-direction and
      //require that the neutron's direction satisfies conservation of momentum. 
      positronInitDir.x = cos(startPhiPositron*(PI/180.0))*sin(startThetaPositron*(PI/180.0));
      positronInitDir.y = sin(startPhiPositron*(PI/180.0))*sin(startThetaPositron*(PI/180.0));
      positronInitDir.z = cos(startThetaPositron*(PI/180.0));

      u = rotMatrix[0][0]*positronInitDir.x + rotMatrix[0][1]*positronInitDir.y +
        rotMatrix[0][2]*positronInitDir.z;
      v = rotMatrix[1][0]*positronInitDir.x + rotMatrix[1][1]*positronInitDir.y +
        rotMatrix[1][2]*positronInitDir.z;
      w = rotMatrix[2][0]*positronInitDir.x + rotMatrix[2][1]*positronInitDir.y +
        rotMatrix[2][2]*positronInitDir.z; 

      mag = sqrt(initialEnergyPositron*initialEnergyPositron + 2.0*initialEnergyPositron*511.0);
      positronInitDir.x = u*mag;
      positronInitDir.y = v*mag;
      positronInitDir.z = w*mag;
 
      neutronInitDir.x = -positronInitDir.x; 
      neutronInitDir.y = -positronInitDir.y;    
      neutronInitDir.z = energyIBD - positronInitDir.z;     

      mag = sqrt(neutronInitDir.x*neutronInitDir.x + neutronInitDir.y*neutronInitDir.y + 
        neutronInitDir.z*neutronInitDir.z);
      neutronInitDir.x = neutronInitDir.x/mag; 
      neutronInitDir.y = neutronInitDir.y/mag;    
      neutronInitDir.z = neutronInitDir.z/mag; 

      //Output the Euler angles for the positron
      initString(string0, stringLen_g);
      sprintf(string0, "  %lf %lf %lf ", alpha*(180.0/PI), acos(cBeta)*(180.0/PI),
        gamma*(180.0/PI));
      fputs(string0, event);

      //Set the Euler angles for the neutron      
      cBeta = neutronInitDir.z;
      sBeta = sqrt(1.0 - cBeta*cBeta);  
          
      cAlpha = neutronInitDir.x/sBeta;
      sAlpha = sqrt(1.0 - cAlpha*cAlpha);
      
      gamma = 2.0*PI*rng(genFlag, rGen);
      cGamma = cos(gamma);
      sGamma = sin(gamma);

      //Set the gamma angle for the neutron and the start time
      initString(string0, stringLen_g);
      sprintf(string0, "%lf\n", gamma*(180.0/PI));
      fputs(string0, event);
       
      initString(string0, stringLen_g);
      sprintf(string0, "  %lf\n", startTimePositron);
      fputs(string0, event);  
      
      //Define the rotation matrix      
      rotMatrix[0][0] = cGamma*cBeta*cAlpha - sGamma*sAlpha;
      rotMatrix[0][1] = cGamma*cBeta*sAlpha + sGamma*cAlpha;
      rotMatrix[0][2] = -cGamma*sBeta;
      rotMatrix[1][0] = -sGamma*cBeta*cAlpha - cGamma*sAlpha;
      rotMatrix[1][1] = -sGamma*cBeta*sAlpha + cGamma*cAlpha;
      rotMatrix[1][2] = sGamma*sBeta;
      rotMatrix[2][0] = sBeta*cAlpha;
      rotMatrix[2][1] = sBeta*sAlpha;
      rotMatrix[2][2] = cBeta;     
      
      //Apply the transformation and add to the detector
      nDepOutside=0;
      nInAfterOut=0;
      depCnt=0;
      for (int j=0; j<nNeutronDep; j++)
      {
        u = rotMatrix[0][0]*depNeutron[j].x + rotMatrix[0][1]*depNeutron[j].y + 
          rotMatrix[0][2]*depNeutron[j].z;
        v = rotMatrix[1][0]*depNeutron[j].x + rotMatrix[1][1]*depNeutron[j].y + 
          rotMatrix[1][2]*depNeutron[j].z;
        w = rotMatrix[2][0]*depNeutron[j].x + rotMatrix[2][1]*depNeutron[j].y + 
          rotMatrix[2][2]*depNeutron[j].z;             

        u = u + delta->x;  v = v + delta->y;  w = w + delta->z;
        depNeutron[j].x=u;  depNeutron[j].y=v;  depNeutron[j].z=w;
                
        check = addDeposit(depNeutron[j], neutronEnergy, neutronTimeFirst, neutronTimeCent);
        if (check==0) { nDepOutside++; }
        else
        {
          if (depNeutron[j].energy>E_MIN)
          {
           if (depCnt==0)
            {
              fprintf(neutronScript, "echo \"#Event %d\" >> outputs/lightTransNeutron.out\n",
                neutronEventNum);
            }
            
            printDeposit(depNeutron[j], depCnt, neutronScript, 4);
            depCnt++;
            if (nDepOutside!=0) nInAfterOut++;
          } //End check if the deposit is above threshold
        } //End check if in the detector 
      } //End loop over neutron deposits
      
      if (VERBOSE_LEVEL==1 && nDepOutside!=0 && nInAfterOut!=0)
      {
        printf("Neutron Event %d:\n", neutronEventNum);
        printf("  nDepOutside=%d\n", nDepOutside);
        printf("  nInAfterOut=%d\n", nInAfterOut);
      }

      if (depCnt!=0)
      {
        fprintf(neutronScript, "NUM_LINES=$(cat outputs/trans-pmt-neutron.dat | wc -l )\n");
        fprintf(neutronScript, "echo $NUM_LINES >> outputs/lightTransNeutron.out\n");
        fprintf(neutronScript, "cat outputs/trans-pmt-neutron.dat >> outputs/lightTransNeutron.out\n");
        
        fprintf(neutronScript, "NUM_LINES=$(cat outputs/trans-time-neutron.dat | wc -l )\n");
        fprintf(neutronScript, "echo $NUM_LINES >> outputs/lightTransNeutron.out\n");          
        fprintf(neutronScript, "cat outputs/trans-time-neutron.dat >> outputs/lightTransNeutron.out\n");  
      }
      else
      {
        fprintf(neutronScript, "echo \"#Event %d\" >> outputs/lightTransNeutron.out\n",
          neutronEventNum);
        fprintf(neutronScript, "echo \"0\" >> outputs/lightTransNeutron.out\n");    
        fprintf(neutronScript, "echo \"0\" >> outputs/lightTransNeutron.out\n");
      }
 
      if (RECON_PROGRAM==0) //If using Bruce's reconstruction then call here
      {
        fprintf(positronScript, "./chargeReconstruction %d %d %d %d\n", NCELL, LATTICE_FLAG,
          MIRROR_FLAG, neutronEventNum);
      } 
      
      //Output the cell hits to the event file. ****Note: keeping the positron cell hits
      //and the neutron cell hits seperate for now

      //Get the time centroids for the cells and the number of hit cells
      nLines=0;
      nPositronCells=0;
      nNeutronCells=0;
      for (int ii=0; ii<NCELL; ii++) {
        for (int j=0; j<NCELL; j++) {
          for (int k=0; k<NCELL; k++) {
            if (positronEnergy[ii][j][k]>E_MIN)
            {
              positronTimeCent[ii][j][k] = positronTimeCent[ii][j][k]/positronEnergy[ii][j][k];
              nPositronCells++;
            }
            
            if (neutronEnergy[ii][j][k]>E_MIN)
            {
              neutronTimeCent[ii][j][k] = neutronTimeCent[ii][j][k]/neutronEnergy[ii][j][k];
              nNeutronCells++;
            }
      } } }       
     
      nLines = nPositronCells + nNeutronCells;
     
      fputs("#Cells_deposits\n", event);
      initString(string0, stringLen_g);
      sprintf(string0, "%d %d %d\n", nLines, nPositronCells, nNeutronCells);
      fputs(string0, event);
      
      //The positron cells
      for (int ii=0; ii<NCELL; ii++) {
        for (int j=0; j<NCELL; j++) {
          for (int k=0; k<NCELL; k++) {
            if (positronEnergy[ii][j][k]>E_MIN)
            {
              initString(string0, stringLen_g);
              sprintf(string0, "  %d %d %d %lf %lf %lf\n", ii+1, j+1, k+1, 
                positronTimeFirst[ii][j][k], positronTimeCent[ii][j][k], positronEnergy[ii][j][k]);
              fputs(string0, event);
            }
      } } } 

      //The neutron cells
      for (int ii=0; ii<NCELL; ii++) {
        for (int j=0; j<NCELL; j++) {
          for (int k=0; k<NCELL; k++) {
            if (neutronEnergy[ii][j][k]>E_MIN)
            {
              initString(string0, stringLen_g);
              sprintf(string0, "  %d %d %d %lf %lf %lf\n", ii+1, j+1, k+1, 
                neutronTimeFirst[ii][j][k], neutronTimeCent[ii][j][k], neutronEnergy[ii][j][k]);
              fputs(string0, event);
            }
      } } }       
      
      //Initialize for the next event
      delete [] depNeutron;
      delete [] depPositron;
      deleteEvent(positronEnergy, positronTimeFirst, positronTimeCent);
      deleteEvent(neutronEnergy, neutronTimeFirst, neutronTimeCent);      
      eventCnt++;      

      initString(string0, stringLen_g);
    }

    //Close files and reinitialize
    if (VERBOSE_LEVEL==0) printf("%d events.\n", eventCnt);
    eventCnt=0;
    fclose(positron);
    fclose(neutron);
    fclose(positronScript);
    fclose(neutronScript);
    fclose(event);
  } //End loop over the positron and neutron files

  //clean-up
  delete delta;
  fclose(list);
  fclose(out);
  
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Set the Light Transport Input &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  transIn = fopen("inputs/trans.dat", "w");
  if (transIn==NULL)
  {
    printf("Could not set the light transport input file!\n");
    return 0;
  }
  
  printDetectorParameters(transIn);

  fclose(transIn);
  
  //&&&&&&&
  //& End &
  //&&&&&&&

  return 0;
}

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&

//*************************************************************************************************

void determineEventStartPoint(int distribFlag, vector_t * delta, int genFlag, gsl_rng * rGen)
{
  //Relocate in the detector
  if (distribFlag==0)       //In the center of the detector. 
  {
  
    delta->x = CELL_DIM*(0.5 - rng(genFlag, rGen));
    delta->y = CELL_DIM*(0.5 - rng(genFlag, rGen));
    delta->z = CELL_DIM*(0.5 - rng(genFlag, rGen));
  }
  else                      //Uniform in the detector
  {
    delta->x = NCELL*CELL_DIM*(0.5 - rng(genFlag, rGen));
    delta->y = NCELL*CELL_DIM*(0.5 - rng(genFlag, rGen));
    delta->z = NCELL*CELL_DIM*(0.5 - rng(genFlag, rGen));
  }
}

//*************************************************************************************************

void printDeposit(deposit_t dep, int depNumber, FILE * f, int particleID)
{
  //Declare needed variables
  char configFileName[stringLen_g];
  double x, y, z;

  //Get the output file names
  initString(configFileName, stringLen_g);

  if (particleID==3)      //Positron
  {
    sprintf(configFileName, "inputs/transConfigPositron.txt");
  }
  else if (particleID==4) //Neutron
  {
    sprintf(configFileName, "inputs/transConfigNeutron.txt");
  }
  else                    //Default behavior
  {
    printf("Wrong call to printDeposit!\n");
    sprintf(configFileName, "inputs/transConfig.txt");
  }
  
  //Get the position in units of cell dimension
  x = 1.0 + NCELL/2.0 + dep.x/CELL_DIM;
  y = 1.0 + NCELL/2.0 + dep.y/CELL_DIM;
  z = 1.0 + NCELL/2.0 + dep.z/CELL_DIM;
  
  //Print the call to the light transport
  if (depNumber>0)
  {
    fprintf(f, "./trans %d %lf %lf %lf %lf %lf %s add\n", dep.id, dep.energy, dep.t*1000.0, x, y, z,
      configFileName);
  }
  else
  {
    fprintf(f, "./trans %d %lf %lf %lf %lf %lf %s\n", dep.id, dep.energy, dep.t*1000.0, x, y, z,
      configFileName);
  }
}

//*************************************************************************************************

void deleteEvent(double energy[][NCELL][NCELL], double timeFirst[][NCELL][NCELL],
  double timeCent[][NCELL][NCELL])
{
  for (int i=0; i<NCELL; i++) {
    for (int j=0; j<NCELL; j++) {
      for (int k=0; k<NCELL; k++) {
        energy[i][j][k] = 0.0;
        timeFirst[i][j][k] = INFINITY;
        timeCent[i][j][k] = 0.0;
  } } }
}

//*************************************************************************************************
  
int addDeposit(deposit_t dep, double energy[][NCELL][NCELL], double timeFirst[][NCELL][NCELL], 
  double timeCent[][NCELL][NCELL])
{
  //Declare needed variables
  int nx, ny, nz;

  //Place in the detector
  nx = int(1.0 + NCELL/2.0 + dep.x/CELL_DIM);
  ny = int(1.0 + NCELL/2.0 + dep.y/CELL_DIM);
  nz = int(1.0 + NCELL/2.0 + dep.z/CELL_DIM);

  //Add the deposit to the detector
  if ( !((nx<1 || nx>NCELL) || (ny<1 || ny>NCELL) || (nz<1 || nz>NCELL)) )
  {
    if (dep.energy>E_MIN)
    {
      energy[nx-1][ny-1][nz-1] += dep.energy;
      if (dep.t<timeFirst[nx-1][ny-1][nz-1]) timeFirst[nx-1][ny-1][nz-1] = dep.t;
      timeCent[nx-1][ny-1][nz-1] += dep.energy*dep.t;
    }
    
    return 1;
  }
  else { return 0; }
}

//*************************************************************************************************
