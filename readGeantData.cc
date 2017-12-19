//&&&&&&&&&&&&&
//& C Headers &
//&&&&&&&&&&&&&
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
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
#define VERBOSE_LEVEL 0         //
#define SCRIPT_VERBOSE_LEVEL 1  //
#define NEW_DEPOSIT_FORMAT 0    //Normally this should be set to 1, unless you are running data in
                                //data/positron_reactor_spectrum_correct directory. In that case set
                                //this macro to 0

//&&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Declarations &
//&&&&&&&&&&&&&&&&&&&&&&&&&

int determineEventStartPoint(int distribFlag, int fiducialVolume, vector_t * delta, int genFlag,
  gsl_rng * rGen);
void printDeposit(deposit_t dep, int depNumber, FILE * f, int fileNum);
void deleteEvent(double energy[][NCELL][NCELL], double timeFirst[][NCELL][NCELL],
  double timeCent[][NCELL][NCELL]);
int addDeposit(deposit_t dep, int index, double energy[][NCELL][NCELL],
  double timeFirst[][NCELL][NCELL], double timeCent[][NCELL][NCELL]);

//*************************************************************************************************

//&&&&&&&&
//& Main &
//&&&&&&&&
int main ()
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Declare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  //Varibles for program configuration
  char basisFileName[stringLen_g];        //TO DO: Go through and comment these
  char dataFileName[stringLen_g];         //
  char scriptDirName[stringLen_g];        //
  char scriptName[stringLen_g];           //
  char eventDirName[stringLen_g];         //
  char eventFileName[stringLen_g];        //
  char string0[stringLen_g];              //
  char string1[stringLen_g];              //
  int genFlag, seed, seedFlag;            //
  int nCells, latticeFlag, mirrorFlag;    //
  int guideFlag, pmtCoupleFlag;           //
  int distribFlag, fiducialVolume;        //
  int quenchingFlag;                      //
  int scinTimeResponse, configCheck;      //

  //Variables for reading the geant data
  int check, nLines, eventNum, nDep;      //
  int nEntries, eventCnt=0, geantSeed;    //
  int startParticleID;                    //
  double initialEnergy, startTheta;       //
  double startPhi;                        //
  double startTime, particleEnergy;       //
  vector_t startPos;                      //
  deposit_t * dep, * tempDep;             //
  
  //Variables for generating events
  int coordinateFlag=0;                       //
  double cellEnergy[NCELL][NCELL][NCELL];     //
  double cellTimeFirst[NCELL][NCELL][NCELL];  //
  double cellTimeCent[NCELL][NCELL][NCELL];   //
  double rotMatrix[3][3];                     //
  double u, v, w;                             //
  int nu, nv, nw;                             //
  vector_t * delta = new vector_t;            //
  double alpha, gamma;                        //
  double cAlpha, sAlpha;                      //
  double cBeta, sBeta;                        //
  double cGamma, sGamma;                      //
  int firstAlpha;                             //
  int inDetector, nDepOutside, nInAfterOut;   //
  int hitCnt=0;                               //
  gsl_rng *  rGen;                            //
  const gsl_rng_type * T;                     //
  
  //Files
  DIR * dir;                                  //
  FILE * config, * list, * data, * script;    //
  FILE * out, * event, * transIn;             //
   

  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&

  //Clear the cell list
  deleteEvent(cellEnergy, cellTimeFirst, cellTimeCent);

  //Read the program control variables
  config = fopen("inputs/readConfig.txt", "r");
  if (config==NULL)
  {
    printf("Could not open the file readConfig.txt! Exiting...\n");
    exit(1);
  }

  if (getNumEntries(config)!=11)
  {
    printf("The file inputs/readConfig.txt is improperly formatted! Exiting...\n");
    exit(1);
  }

  for (int i=0; i<11; i++)
  {
    initString(string0, stringLen_g);
    initString(string1, stringLen_g);
    getEntry(config, string0, stringLen_g);  
    
    if      (i==0)  { sscanf(string0, "%s %d", string1, &genFlag);          }
    else if (i==1)  { sscanf(string0, "%s %d", string1, &seedFlag);         }
    else if (i==2)  { sscanf(string0, "%s %d", string1, &nCells);           }
    else if (i==3)  { sscanf(string0, "%s %d", string1, &latticeFlag);      }
    else if (i==4)  { sscanf(string0, "%s %d", string1, &mirrorFlag);       }
    else if (i==5)  { sscanf(string0, "%s %d", string1, &guideFlag);        }
    else if (i==6)  { sscanf(string0, "%s %d", string1, &pmtCoupleFlag);    }
    else if (i==7)  { sscanf(string0, "%s %d", string1, &distribFlag);      }
    else if (i==8)  { sscanf(string0, "%s %d", string1, &fiducialVolume);   }
    else if (i==9)  { sscanf(string0, "%s %d", string1, &quenchingFlag);    }
    else            { sscanf(string0, "%s %d", string1, &scinTimeResponse); }
  }
  
  fclose(config);
  
  if (fiducialVolume<0) { fiducialVolume = 0; }
  else if (fiducialVolume>(NCELL/2 - ((NCELL-1)%2)))
  {
    fiducialVolume = 0;
    distribFlag = 0;
  }

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
    printf("Could not open the file! Exiting...\n");
    exit(1);
  }


  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Confirm the Configuration &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  if (NCELL!=nCells || LATTICE_FLAG!=latticeFlag || MIRROR_FLAG!=mirrorFlag ||
    GUIDE_FLAG!=guideFlag || PMT_COUPLING_FLAG!=pmtCoupleFlag)
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
  if (distribFlag==-1)
  {
    printf("  The events will be generated in the EXACT CENTER of the detector.\n");    
  }
  else if (distribFlag==0)
  {
    printf("  The events will be generated in the CENTER of the detector.\n");  
  }
  else if (distribFlag==1)
  {
    if (fiducialVolume==0) { printf("  The events will be generated UNIFORM in the detector.\n"); }
    else
    {
      printf("  The events will be generated UNIFORM in the inner %dx%dx%d cells.\n",
        NCELL-2*fiducialVolume, NCELL-2*fiducialVolume, NCELL-2*fiducialVolume);
    }
  }
  else if (distribFlag==2)
  {
    printf("  The events will be generated UNIFORM on the SIDES of the detector.\n");    
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
  }
  else
  {
    printf("  Light guides are OFF.\n");
    guideFlag=0;  //To ensure a proper call to the light transport function
  }

  //Pmt optical coupling
  if (pmtCoupleFlag==1) { printf("  The PMTs are optically coupled.\n"); }
  else { printf("  The PMTs are NOT optically coupled.\n"); }

  //Check if there is basis data to do the reconstruction with 
  
  //TO DO: Implement the coupling/non-coupling option into the event reconstruction
  
  initString(basisFileName, stringLen_g);
  sprintf(basisFileName, "reconstructionData/reconBasisData_%d_%d%d%d_%d_%d.dat", NCELL,
  LATTICE_FLAG, MIRROR_FLAG, GUIDE_FLAG, TIME_BIN_SIZE, N_TIME_BIN_PRINT);
  
  struct stat fileStatus;
  if(stat(basisFileName, &fileStatus)==-1)
  {
    printf("\n  WARNING! The required basis data for reconstruction does not exist! Proceeding anyway...\n");
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

  list = fopen("./inputs/geantDataList.txt", "r");
  if (list==NULL)
  {
    printf("Cannot open the file! Exiting...\n");
    exit(1); 
  }
  nEntries = getNumEntries(list);
   

  //&&&&&&&&&&&&&&&&&
  //& Read the Data &
  //&&&&&&&&&&&&&&&&&

  for (int i=0; i<nEntries; i++)
  {
    //Get the file name and open
    getEntry(list, dataFileName, stringLen_g);

    data = fopen(dataFileName, "r");
    if (data==NULL)
    {
      printf("Could not open %s! Skipping to the next file in the list.\n", dataFileName);
      continue;
    }
    if (VERBOSE_LEVEL==0) printf("%s has ", dataFileName);

    //Check script directory, name it similar the data directory
    copyString(string0, dataFileName, stringLen_g);
    removePath(string0, 1);
    getDirName(string1, string0, 0);
    initString(scriptDirName, stringLen_g);
    sprintf(scriptDirName, "simulationScripts/");
    strcat(scriptDirName, string1);

    dir = opendir(scriptDirName);
    if (dir) { closedir(dir); } //Directory exists -> close and move on
    else if (errno==ENOENT)     //Directory does not exist -> create
    {
      check = mkdir(scriptDirName, 0777);
      if (check)
      {
        printf("Error in opening script directory! Skipping to the next file in the list.\n");
        fclose(data);
        continue;
      }
    }
    else                        //opendir() failed for some other reason
    {
      printf("Error in opening script directory! Skipping to the next file in the list.\n");
      fclose(data);
      continue;
    }
    strcat(scriptDirName, "/");

    //Open the bash script and name the script similar to the data file name
    copyString(scriptName, dataFileName, stringLen_g);
    removePath(scriptName, 0);
    removeFileExtension(scriptName, 4);
    strcat(scriptName, "_Run.sh");
    initString(string0, stringLen_g);
    copyString(string0, scriptDirName, stringLen_g);
    strcat(string0, scriptName);
    copyString(scriptName, string0, stringLen_g);

    fprintf(out, "%s\n", scriptName);
    script = fopen(scriptName, "w");
    if (script==NULL)
    {
      printf("Could not open %s! Skipping to the next file in the list\n", scriptName);
      fclose(data);
      continue;
    }

    chmod(scriptName, 0777); //Make the script executable

    //Add the light transport configuration file
    int scriptIndex = i-int(i/N_SPLIT)*N_SPLIT; 
    char configFileName[stringLen_g];
    
    initString(configFileName, stringLen_g);   
    if (N_SPLIT==1) sprintf(configFileName, "inputs/transConfig.txt");
    else sprintf(configFileName, "inputs/transConfig%d.txt", scriptIndex);
    
    fprintf(script, "rm %s\n", configFileName);
    fprintf(script, "echo \"%d\" > %s\n", NCELL, configFileName);
    fprintf(script, "echo \"%d\" >> %s\n", LATTICE_FLAG, configFileName);
    fprintf(script, "echo \"%d\" >> %s\n", MIRROR_FLAG, configFileName);
    fprintf(script, "echo \"%d\" >> %s\n", GUIDE_FLAG, configFileName);
    fprintf(script, "echo \"%d\" >> %s\n", PMT_COUPLING_FLAG, configFileName);
    
    if (N_SPLIT==1)
    {
      fprintf(script, "echo \"'outputs/trans-pmt.dat'\" >> %s\n", configFileName);
      fprintf(script, "echo \"'outputs/trans-time.dat'\" >> %s\n", configFileName);
    }           //Note: you need the single quotes here because of Fortan's poor string reading
    else
    {
      fprintf(script, "echo \"'outputs/trans-pmt%d.dat'\" >> %s\n", scriptIndex, configFileName);
      fprintf(script, "echo \"'outputs/trans-time%d.dat'\" >> %s\n", scriptIndex, configFileName);
    }
    
    //Add light transport output manipulation stuff
    if (N_SPLIT==1) fprintf(script, "rm outputs/lightTrans.out\n");
    else fprintf(script, "rm outputs/lightTrans%d.out\n", scriptIndex);
    
    if (RECON_PROGRAM==0) fprintf(script, "rm outputs/reconstruction.out\n"); //If using Bruce's
                                                                              //reconstrution
    
    //Open the event file directory
    initString(eventDirName, stringLen_g);
    sprintf(eventDirName, "eventFiles/"); 
    copyString(string0, dataFileName, stringLen_g);
    removePath(string0, 1);
    initString(string1, stringLen_g);
    getDirName(string1, string0, 0);
    strcat(eventDirName, string1);

    dir = opendir(eventDirName);
    if (dir) { closedir(dir); } //Directory exists -> close and move on
    else if (errno==ENOENT)     //Directory does not exist -> create
    {
      check = mkdir(eventDirName, 0777);
      if (check)
      {
        printf("Error in opening eventFile directory! Skipping to the next file in the list.\n");
        fclose(data);
        fclose(script);
        continue;
      }
    }
    else                        //opendir() failed for some other reason
    {
      printf("Error in opening eventFile directory! Skipping to the next file in the list.\n");
      fclose(data);
      fclose(script);
      continue;
    }
    strcat(eventDirName, "/");

    //Open the event file
    copyString(eventFileName, dataFileName, stringLen_g);
    removePath(eventFileName, 0);
    removeFileExtension(eventFileName, 4);
    strcat(eventFileName, "_Event.dat");
    copyString(string0, eventDirName, stringLen_g);
    strcat(string0, eventFileName);
    copyString(eventFileName, string0, stringLen_g);

    event = fopen(eventFileName, "w");
    if (event==NULL)
    {
      printf("Could not open %s! Skipping to the next file in the list\n", string0);
      fclose(data);
      fclose(script);
      continue;
    }
    
    //Check if this is a neutron data file and if so it is for boron or lithium loading
    int boronOrLithium=0;
    if (strstr(dataFileName, "boron")) { boronOrLithium = 1; }
    else if (strstr(dataFileName, "lithium")) { boronOrLithium = 2; }
    
    //Loop over the data and generate the script
    initString(string0, stringLen_g);
    while (fgets(string0, stringLen_g, data)!=NULL)
    {
      fputs(string0, event);

      if (strstr(string0, "#Event")!=NULL) //Check if at the event header
      {
        //Read the event information
        fgets(string0, stringLen_g, data);
        fputs(string0, event);
        sscanf(string0, "%d", &nLines);
        for (int j=0; j<nLines; j++)
        {
          fgets(string0, stringLen_g, data);
          fputs(string0, event);
          if (j==0) { sscanf(string0, "%d", &eventNum); }
          else if (j==1) { sscanf(string0, "%d", &geantSeed); }
        }

        //Read the event properties
        initString(string0, stringLen_g);
        fgets(string0, stringLen_g, data);
        fputs(string0, event);
        if (strstr(string0, "#Geant4_properties")==NULL) { continue; }
        
        fgets(string0, stringLen_g, data);
        fputs(string0, event);
        sscanf(string0, "%d", &nLines);
        for (int j=0; j<nLines; j++)
        {
          fgets(string0, stringLen_g, data);
          fputs(string0, event);
          if (j==0) { sscanf(string0, "%d", &startParticleID); }
          else if (j==1) { sscanf(string0, "%lf", &initialEnergy); }
          else if (j==2) { sscanf(string0, "%lf %lf %lf", &startPos.x, &startPos.y, &startPos.z); }
          else if (j==3) { sscanf(string0, "%lf %lf", &startTheta, &startPhi); }
          else if (j==4) { sscanf(string0, "%lf", &startTime); }
        }

        //Read the deposits
        initString(string0, stringLen_g);
        fgets(string0, stringLen_g, data);
        fputs(string0, event);
        if (strstr(string0, "#Energy_deposits")==NULL) { continue; }

        fgets(string0, stringLen_g, data);
        fputs(string0, event);
        sscanf(string0, "%d", &nLines);
        
        dep = new deposit_t [nLines];
        firstAlpha=1;
        for (int j=0; j<nLines; j++)
        {
          initString(string0, stringLen_g);
          fgets(string0, stringLen_g, data);
          if (NEW_DEPOSIT_FORMAT)
          {
            sscanf(string0, "%d %lf %lf %lf %lf %lf %lf", 
              &dep[j].id, &u, &v, &w, &dep[j].t, &particleEnergy, &dep[j].energy);
          }
          else
          {
            sscanf(string0, "%d %lf %lf %lf %lf %lf", 
              &dep[j].id, &u, &v, &w, &dep[j].t, &dep[j].energy);
            particleEnergy = 0.0;
          }

          //Adjust energy to reflect quenching in the scintillator
          if (quenchingFlag==0 && boronOrLithium!=0)
          {
            if (dep[j].id>3 || dep[j].id==0) //Zero energy deposits not an electron, positron,
            {                                //or gamma
              dep[j].energy = 0.0;                                    
            }
            
            if (dep[j].id==8 && firstAlpha==1) //If the particle is an alpha set it to the
            {                                  //electron equivalent for B-10 or Li-6 capture
              if (boronOrLithium==1) { dep[j].energy = 0.07; }
              else                   { dep[j].energy = 0.40; }
              firstAlpha=0;
            }
          }

          initString(string0, stringLen_g);
          sprintf(string0, "  %d %lf %lf %lf %lf %lf %lf\n", dep[j].id, u, v, w, dep[j].t,
            particleEnergy, dep[j].energy);
          fputs(string0, event);

          if (coordinateFlag==0 || coordinateFlag==5 || coordinateFlag==6) //xyz -> xyz
          {
            dep[j].x=u/10.0;  dep[j].y=v/10.0; dep[j].z=w/10.0;  //Convert to cm
          }
        } //End loop over the Geant4 hits

        //Set the deposit id's so that the light transport generates the scintillation light with
        //the correct time dependence
        if (scinTimeResponse==0)  //Pure exponential with a 2.2 ns decay time
        {
          for (int depNum=0; depNum<nLines; depNum++) { dep[depNum].id = 0; }
        }
        else                      //Uses the gamma, fast neutron, or thermal neutron time response
        {                         //from "Comparative analysis of pulse shape discrimination..."
          for (int depNum=0; depNum<nLines; depNum++)
          {
            if      (dep[depNum].id<=3) { dep[depNum].id = 1; }     //Electrons/gammas
            else if (dep[depNum].id==4) { dep[depNum].id = 2; }     //Proton scatters
            else if (dep[depNum].id==8) { dep[depNum].id = 3; }     //Neutron captures
            else                        { dep[depNum].energy=0.0; } //Axe deposits from other than
          }                                                         //e-, e+, gamma, proton, or alpha
        }

        //Remove deposits that have no energy    
        nDep=0;
        for (int depNum=0; depNum<nLines; depNum++)
        {
          if (dep[depNum].energy>0.0) nDep++;
        }
        
        if (nDep<nLines)
        {
          tempDep = new deposit_t [nDep];
          nDep=0;
          for (int depNum=0; depNum<nLines; depNum++)
          {
            if (dep[depNum].energy>0.0)
            {
              tempDep[nDep].id = dep[depNum].id;
              tempDep[nDep].x = dep[depNum].x;
              tempDep[nDep].y = dep[depNum].y;
              tempDep[nDep].z = dep[depNum].z;
              tempDep[nDep].t = dep[depNum].t;
              tempDep[nDep].energy = dep[depNum].energy;
              nDep++;
            }
          }
          
          delete [] dep;
          dep = new deposit_t [nDep];
          for (int depNum=0; depNum<nDep; depNum++)
          {
            dep[depNum].id = tempDep[depNum].id;
            dep[depNum].x = tempDep[depNum].x;
            dep[depNum].y = tempDep[depNum].y;
            dep[depNum].z = tempDep[depNum].z;
            dep[depNum].t = tempDep[depNum].t;
            dep[depNum].energy = tempDep[depNum].energy;
          }   
          delete [] tempDep;     
          nLines = nDep;
        }
        
        //Adjust the timing relative to the first event
        sortDeposits(dep, nLines, 1); //Sort by time with the earliest deposits first
        
        for (int depNum=1; depNum<nLines; depNum++)
        {
          dep[depNum].t = dep[depNum].t - dep[0].t;
        }
        dep[0].t=0.0;


        //Transform the deposits
        if (distribFlag<3)
        {
          //I think that requiring that the deposit first in time be in the detector is safe and
          //unquestionably physical, but killing a track after a deposit outside the detector may
          //not work out so well because that deposit could have been from a secondary and the
          //primary has future deposits in the detector.

          double r = sqrt(dep[0].x*dep[0].x + dep[0].y*dep[0].y + dep[0].z*dep[0].z);
          double rmax;
          if (distribFlag>0) rmax = sqrt(3.0)*(NCELL - 2*fiducialVolume)*CELL_DIM;
          else rmax = (double(NCELL)/2.0)*CELL_DIM/0.95;

          if (r<0.95*rmax)        //First deposit CAN be in the dectector, 0.95 was a nice fudge
          {                       //factor to include             
            inDetector=0;
            while (!inDetector) //Require that the first deposit is in the detector
            {
              //Set the offset
              delta->x=0.0; delta->y=0.0; delta->z=0.0;
              coordinateFlag = determineEventStartPoint(distribFlag, fiducialVolume, delta, genFlag, rGen);
             
              //Define the random Euler angles
              alpha = 2.0*PI*rng(genFlag, rGen);
              if (coordinateFlag==0) cBeta = 2.0*rng(genFlag, rGen) - 1.0;
              else cBeta = rng(genFlag, rGen);
              gamma = 2.0*PI*rng(genFlag, rGen);
              
              cAlpha = cos(alpha);
              sAlpha = sin(alpha);
              sBeta = sqrt(1.0 - cBeta*cBeta);
              cGamma = cos(gamma);
              sGamma = sin(gamma);
              
              if (distribFlag==-1)  //If testing use the following Euler angles
              {
                cAlpha = 1.0; //Make the rotation matrix the idenity
                sAlpha = 0.0;
                cBeta = 1.0;
                sBeta = 0.0;
                cGamma = 1.0;
                sGamma = 0.0;            
              }

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

              if (coordinateFlag%2==0 && coordinateFlag!=0)  //On a NCELL+1 side
              {
                for (int rowIndex=0; rowIndex<3; rowIndex++)
                {
                  rotMatrix[rowIndex][2] = -1.0*rotMatrix[rowIndex][2];
                }        
              }
              
              //Apply the transformation and check the first deposit
              u = rotMatrix[0][0]*dep[0].x + rotMatrix[0][1]*dep[0].y + rotMatrix[0][2]*dep[0].z;
              v = rotMatrix[1][0]*dep[0].x + rotMatrix[1][1]*dep[0].y + rotMatrix[1][2]*dep[0].z;
              w = rotMatrix[2][0]*dep[0].x + rotMatrix[2][1]*dep[0].y + rotMatrix[2][2]*dep[0].z;

              u = u + delta->x;  v = v + delta->y;  w = w + delta->z;
              
              nu = 1 + NCELL/2.0 + u/CELL_DIM;
              nv = 1 + NCELL/2.0 + v/CELL_DIM;
              nw = 1 + NCELL/2.0 + w/CELL_DIM;

              if ( !((nu<1 || nu>NCELL) || (nv<1 || nv>NCELL) || (nw<1 || nw>NCELL)) )
              {
                inDetector=1;
              }
            } //End get transformation
            
            //Set the light transport properties
            fputs("#Light_transport_properties\n3\n", event);

            initString(string0, stringLen_g);
            sprintf(string0, "  %lf %lf %lf\n", delta->x+startPos.x, delta->y+startPos.y, 
              delta->z+startPos.z);
            fputs(string0, event);

            initString(string0, stringLen_g);
            sprintf(string0, "  %lf %lf %lf\n", alpha*(180.0/PI), acos(cBeta)*(180.0/PI),
              gamma*(180.0/PI));
            fputs(string0, event);
            
            initString(string0, stringLen_g);
            sprintf(string0, "  %lf\n", startTime);
            fputs(string0, event);
            
            //Apply the transformation and add deposits to the detector
            nDepOutside=0;
            nInAfterOut=0;
            hitCnt=0;
            for (int j=0; j<nLines; j++)
            {
              u = rotMatrix[0][0]*dep[j].x + rotMatrix[0][1]*dep[j].y + rotMatrix[0][2]*dep[j].z;
              v = rotMatrix[1][0]*dep[j].x + rotMatrix[1][1]*dep[j].y + rotMatrix[1][2]*dep[j].z;
              w = rotMatrix[2][0]*dep[j].x + rotMatrix[2][1]*dep[j].y + rotMatrix[2][2]*dep[j].z;

              dep[j].x = u + delta->x;  dep[j].y = v + delta->y;  dep[j].z = w + delta->z;
              check = addDeposit(dep[j], j, cellEnergy, cellTimeFirst, cellTimeCent);
              
              if (check==0) { nDepOutside++; }
              else
              {
                if (hitCnt==0)
                {
                  if (N_SPLIT==1)
                  {
                    fprintf(script, "echo \"#Event %d\" >> outputs/lightTrans.out\n", eventNum);
                  }
                  else
                  {
                    fprintf(script, "echo \"#Event %d\" >> outputs/lightTrans%d.out\n", eventNum,
                      i-int(i/N_SPLIT)*N_SPLIT);
                  }
                                   
                  if (SCRIPT_VERBOSE_LEVEL==2)
                  {
                    if ((eventNum+1)%100==0) fprintf(script, "echo \"%d events processed\"\n",
                      eventNum+1);
                  }
                } //End check if the first deposit

                printDeposit(dep[j], j, script, i);
                hitCnt++;
                if (nDepOutside!=0) nInAfterOut++;
              } //End check if in the detector            
            } //End loop over deposits
                  
            if (VERBOSE_LEVEL==1 && nDepOutside!=0 && nInAfterOut!=0)
            {
              printf("Event %d:\n", eventNum);
              printf("  nDepOutside=%d\n", nDepOutside);
              printf("  nInAfterOut=%d\n", nInAfterOut);
            }
            
            //Print save the light transport output by getting the number of lines in the output files
            //and then concatenate the full output file with the number of lines and the contents.            
            
            if (hitCnt>0)
            {
              if (N_SPLIT==1)
              {
                //First to the PMT hit and their charges
                fprintf(script, "NUM_LINES=$(cat outputs/trans-pmt.dat | wc -l )\n");
                fprintf(script, "echo $NUM_LINES >> outputs/lightTrans.out\n");
                fprintf(script, "cat outputs/trans-pmt.dat >> outputs/lightTrans.out\n");
                
                //Now do the time bins of the hit PMTs
                fprintf(script, "NUM_LINES=$(cat outputs/trans-time.dat | wc -l )\n");
                fprintf(script, "echo $NUM_LINES >> outputs/lightTrans.out\n");          
                fprintf(script, "cat outputs/trans-time.dat >> outputs/lightTrans.out\n");
                
                //If using Bruce's reconstruction then call here
                if (RECON_PROGRAM==0)
                {
                  fprintf(script, "./chargeReconstruction %d %d %d %d\n", NCELL, LATTICE_FLAG,
                    MIRROR_FLAG, eventNum);
                }
              }
              else
              {
                //First to the PMT hit and their charges
                fprintf(script, "NUM_LINES=$(cat outputs/trans-pmt%d.dat | wc -l )\n", scriptIndex);
                fprintf(script, "echo $NUM_LINES >> outputs/lightTrans%d.out\n", scriptIndex);
                fprintf(script, "cat outputs/trans-pmt%d.dat >> outputs/lightTrans%d.out\n",
                  scriptIndex, scriptIndex);
                
                //Now do the time bins of the hit PMTs
                fprintf(script, "NUM_LINES=$(cat outputs/trans-time%d.dat | wc -l )\n", scriptIndex);
                fprintf(script, "echo $NUM_LINES >> outputs/lightTrans%d.out\n", scriptIndex);          
                fprintf(script, "cat outputs/trans-time%d.dat >> outputs/lightTrans%d.out\n",
                  scriptIndex, scriptIndex);
                  
                if (RECON_PROGRAM==0)
                {
                  printf("The option to use Bruce's reconstruction with split scripts is not supported.\n");
                }
              }
            }
            else
            {
              if (N_SPLIT==1)
              {
                fprintf(script, "echo \"#Event %d\" >> outputs/lightTrans.out\n", eventNum);
                fprintf(script, "echo \"0\" >> outputs/lightTrans.out\n");    
                fprintf(script, "echo \"0\" >> outputs/lightTrans.out\n");                        
              }
              else
              {
                fprintf(script, "echo \"#Event %d\" >> outputs/lightTrans%d.out\n", eventNum,
                  scriptIndex);
                fprintf(script, "echo \"0\" >> outputs/lightTrans%d.out\n", scriptIndex);    
                fprintf(script, "echo \"0\" >> outputs/lightTrans%d.out\n", scriptIndex);               
              }  
            } //End check if there were cell hits
            
            //Output the cell hits to the event file 
     
            //Get the time centroids for the cells and the number of hit cells
            nLines=0;
            for (int ii=0; ii<NCELL; ii++) {
              for (int j=0; j<NCELL; j++) {
                for (int k=0; k<NCELL; k++) {
                  if (cellEnergy[ii][j][k]>E_MIN)
                  {
                    cellTimeCent[ii][j][k] = cellTimeCent[ii][j][k]/cellEnergy[ii][j][k];
                    nLines++;
                  }
            } } } 

            fputs("#Cells_deposits\n", event);
            initString(string0, stringLen_g);
            sprintf(string0, "%d\n", nLines);
            fputs(string0, event);
            for (int ii=0; ii<NCELL; ii++) {
              for (int j=0; j<NCELL; j++) {
                for (int k=0; k<NCELL; k++) {
                  if (cellEnergy[ii][j][k]>E_MIN)
                  {
                    initString(string0, stringLen_g);
                    sprintf(string0, "  %d %d %d %lf %lf %lf\n", ii+1, j+1, k+1,
                      cellTimeFirst[ii][j][k], cellTimeCent[ii][j][k], cellEnergy[ii][j][k]);
                    fputs(string0, event);
                  }
            } } } 

            deleteEvent(cellEnergy, cellTimeFirst, cellTimeCent);
            eventCnt++;
          }
          else  //The first deposit is not in the detector
          {  
            fputs("#Light_transport_properties\n3\n", event);

            initString(string0, stringLen_g);
            sprintf(string0, "  %lf %lf %lf\n", startPos.x, startPos.y, startPos.z);
            fputs(string0, event);

            initString(string0, stringLen_g);
            sprintf(string0, "  %lf %lf\n", startTheta, startPhi);
            fputs(string0, event);
            
            initString(string0, stringLen_g);
            sprintf(string0, "  %lf\n", startTime);
            fputs(string0, event);
            
            fputs("#Cells_deposits\n0\n", event);
            
            if (N_SPLIT==1)
            {
              fprintf(script, "echo \"#Event %d\" >> outputs/lightTrans.out\n", eventNum);
              fprintf(script, "echo \"0\" >> outputs/lightTrans.out\n");    
              fprintf(script, "echo \"0\" >> outputs/lightTrans.out\n");                        
            }
            else
            {
              fprintf(script, "echo \"#Event %d\" >> outputs/lightTrans%d.out\n", eventNum,
                scriptIndex);
              fprintf(script, "echo \"0\" >> outputs/lightTrans%d.out\n", scriptIndex);    
              fprintf(script, "echo \"0\" >> outputs/lightTrans%d.out\n", scriptIndex);               
            }              

          } //End check if the first deposit can be in the detector
        }
        else  //Use the deposits directly from Geant
        {
          //Set the light transport properties
          fputs("#Light_transport_properties\n3\n", event);
          
          initString(string0, stringLen_g);
          sprintf(string0, "  0.0 0.0 0.0\n");
          fputs(string0, event);

          initString(string0, stringLen_g);
          sprintf(string0, "  0.0 1.0 0.0\n");
          fputs(string0, event);
          
          initString(string0, stringLen_g);
          sprintf(string0, "  0.0\n");
          fputs(string0, event);
          
          //Add deposits to the cells
          hitCnt=0;
          for (int j=0; j<nLines; j++)
          {
            check = addDeposit(dep[j], j, cellEnergy, cellTimeFirst, cellTimeCent);
            
            if (check==0) { printf("WARNING! Possible use of incompatible data!\n"); }
            else
            {
              if (dep[j].energy>E_MIN)
              {
                if (hitCnt==0)
                {
                  if (N_SPLIT==1)
                  {
                    fprintf(script, "echo \"#Event %d\" >> outputs/lightTrans.out\n", eventNum);
                  }
                  else
                  {
                    fprintf(script, "echo \"#Event %d\" >> outputs/lightTrans%d.out\n", eventNum,
                      i-int(i/N_SPLIT)*N_SPLIT);
                  }

                  if (SCRIPT_VERBOSE_LEVEL==2)
                  {
                    if ((eventNum+1)%100==0) fprintf(script, "echo \"%d events processed\"\n",
                      eventNum+1);
                  }
                } //End check if the first deposit

                printDeposit(dep[j], j, script, i);
                hitCnt++;
                                                
              } //End check if the deposit is above threshold
            } //End check if in the detector            
          } //End loop over deposits          
            
          //Print save the light transport output by getting the number of lines in the output files
          //and then concatenate the full output file with the number of lines and the contents.
          
          if (hitCnt>0)
          {
            if (N_SPLIT==1)
            {
              //First to the PMT hit and their charges
              fprintf(script, "NUM_LINES=$(cat outputs/trans-pmt.dat | wc -l )\n");
              fprintf(script, "echo $NUM_LINES >> outputs/lightTrans.out\n");
              fprintf(script, "cat outputs/trans-pmt.dat >> outputs/lightTrans.out\n");
              
              //Now do the time bins of the hit PMTs
              fprintf(script, "NUM_LINES=$(cat outputs/trans-time.dat | wc -l )\n");
              fprintf(script, "echo $NUM_LINES >> outputs/lightTrans.out\n");          
              fprintf(script, "cat outputs/trans-time.dat >> outputs/lightTrans.out\n");
              
              //If using Bruce's reconstruction then call here
              if (RECON_PROGRAM==0)
              {
                fprintf(script, "./chargeReconstruction %d %d %d %d\n", NCELL, LATTICE_FLAG,
                  MIRROR_FLAG, eventNum);
              }
            }
            else
            {
              fprintf(script, "NUM_LINES=$(cat outputs/trans-pmt%d.dat | wc -l )\n", scriptIndex);
              fprintf(script, "echo $NUM_LINES >> outputs/lightTrans%d.out\n", scriptIndex);
              fprintf(script, "cat outputs/trans-pmt%d.dat >> outputs/lightTrans%d.out\n",
                scriptIndex, scriptIndex);
              
              fprintf(script, "NUM_LINES=$(cat outputs/trans-time%d.dat | wc -l )\n", scriptIndex);
              fprintf(script, "echo $NUM_LINES >> outputs/lightTrans%d.out\n", scriptIndex);          
              fprintf(script, "cat outputs/trans-time%d.dat >> outputs/lightTrans%d.out\n",
                scriptIndex, scriptIndex);
                
              if (RECON_PROGRAM==0)
              {
                printf("The option to use Bruce's reconstruction with split scripts is not supported.\n");
              }
            }
          }
          else
          {
            if (N_SPLIT==1)
            {
              fprintf(script, "echo \"#Event %d\" >> outputs/lightTrans.out\n", eventNum);
              fprintf(script, "echo \"0\" >> outputs/lightTrans.out\n");    
              fprintf(script, "echo \"0\" >> outputs/lightTrans.out\n");                        
            }
            else
            {
              fprintf(script, "echo \"#Event %d\" >> outputs/lightTrans%d.out\n", eventNum,
                scriptIndex);
              fprintf(script, "echo \"0\" >> outputs/lightTrans%d.out\n", scriptIndex);    
              fprintf(script, "echo \"0\" >> outputs/lightTrans%d.out\n", scriptIndex);               
            }  
          } //End check if there were cell hits
          
          //Output the cell hits to the event file 
   
          //Get the time centroids for the cells and the number of hit cells
          nLines=0;
          for (int ii=0; ii<NCELL; ii++) {
            for (int j=0; j<NCELL; j++) {
              for (int k=0; k<NCELL; k++) {
                if (cellEnergy[ii][j][k]>E_MIN)
                {
                  cellTimeCent[ii][j][k] = cellTimeCent[ii][j][k]/cellEnergy[ii][j][k];
                  nLines++;
                }
          } } } 

          fputs("#Cells_deposits\n", event);
          initString(string0, stringLen_g);
          sprintf(string0, "%d\n", nLines);
          fputs(string0, event);
          for (int ii=0; ii<NCELL; ii++) {
            for (int j=0; j<NCELL; j++) {
              for (int k=0; k<NCELL; k++) {
                if (cellEnergy[ii][j][k]>E_MIN)
                {
                  initString(string0, stringLen_g);
                  sprintf(string0, "  %d %d %d %lf %lf %lf\n", ii+1, j+1, k+1,
                    cellTimeFirst[ii][j][k], cellTimeCent[ii][j][k], cellEnergy[ii][j][k]);
                  fputs(string0, event);
                }
          } } } 

          deleteEvent(cellEnergy, cellTimeFirst, cellTimeCent);
          eventCnt++;      

        } //End check if using the unmodified Geant deposits

        //printf("%d\n", eventCnt);            
        delete [] dep;
      } //End  check if at an event header

      initString(string0, stringLen_g);
    } //End the loop through the file

    //Close files and reinitialize
    if (VERBOSE_LEVEL==0) printf("%d events.\n", eventCnt);
    eventCnt=0;
    fclose(data);
    if (SCRIPT_VERBOSE_LEVEL>=1) fprintf(script, "echo \"%s processed!\"\n", scriptName);
    fclose(script);
    fclose(event);
  } //End loop over the files

  //clean-up
  if (genFlag!=0) gsl_rng_free(rGen);
  
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

int determineEventStartPoint(int distribFlag, int fiducialVolume, vector_t * delta, int genFlag,
  gsl_rng * rGen)
{
  int coordinateFlag=0;

  //Relocate in the detector
  if (distribFlag==0)       //In the center of the detector. 
  {
  
    delta->x = CELL_DIM*(0.5 - rng(genFlag, rGen));
    delta->y = CELL_DIM*(0.5 - rng(genFlag, rGen));
    delta->z = CELL_DIM*(0.5 - rng(genFlag, rGen));
  }
  else if (distribFlag==1)  //Uniform in the detector
  {
    delta->x = (NCELL-2*fiducialVolume)*CELL_DIM*(0.5 - rng(genFlag, rGen));
    delta->y = (NCELL-2*fiducialVolume)*CELL_DIM*(0.5 - rng(genFlag, rGen));
    delta->z = (NCELL-2*fiducialVolume)*CELL_DIM*(0.5 - rng(genFlag, rGen));
  }
  else if (distribFlag==2)
  {
  
    //TO DO: Implement code so that external gammas can be generated on all sides
  
    int side; //Use Bruce's light transport notation:
              //side==1 -> x=0 side
              //side==2 -> x=NCELL+1 side
              //side==3 -> y=0 side
              //side==4 -> y=NCELL+1 side
              //side==5 -> z=0 side
              //side==6 -> z=NCELL+1 side
              
    if (rng(genFlag, rGen)<0.5) side=5;
    else side=6;

    coordinateFlag = side;
    double u = NCELL*CELL_DIM*(0.5 - rng(genFlag, rGen));
    double v = NCELL*CELL_DIM*(0.5 - rng(genFlag, rGen));    

    delta->x = u;
    delta->y = v;
    delta->z = (side-5)*NCELL*CELL_DIM - (NCELL*CELL_DIM)/2.0;                                                
  }
  else if (distribFlag==-1)      //Place in the exact center of the detector for testing 
  {
    delta->x = 0.0;
    delta->y = 0.0;
    delta->z = 0.0;
  }
    
  return coordinateFlag;
}

//*************************************************************************************************

void printDeposit(deposit_t dep, int depNumber, FILE * f, int fileNum)
{
  //Declare needed variables
  double x, y, z;
  int scriptIndex = fileNum - int(fileNum/N_SPLIT)*N_SPLIT;
  
  //Get the position in units of cell dimension
  x = 1.0 + NCELL/2.0 + dep.x/CELL_DIM;
  y = 1.0 + NCELL/2.0 + dep.y/CELL_DIM;
  z = 1.0 + NCELL/2.0 + dep.z/CELL_DIM;
  
  //Print the call to the light transport
  if (depNumber>0)
  {
    if (N_SPLIT==1)
    {
      fprintf(f, "./trans %d %lf %lf %lf %lf %lf inputs/transConfig.txt add\n", dep.id, dep.energy,
        dep.t*1000.0, x, y, z);
    }
    else
    {
      fprintf(f, "./trans %d %lf %lf %lf %lf %lf inputs/transConfig%d.txt add\n", dep.id,
        dep.energy, dep.t*1000.0, x, y, z, scriptIndex);    
    }
  }
  else
  {
    if (N_SPLIT==1)
    {
      fprintf(f, "./trans %d %lf %lf %lf %lf %lf inputs/transConfig.txt\n", dep.id, dep.energy,
        dep.t*1000.0, x, y, z);
    }
    else
    {
      fprintf(f, "./trans %d %lf %lf %lf %lf %lf inputs/transConfig%d.txt\n", dep.id, dep.energy,
        dep.t*1000.0, x, y, z, scriptIndex);   
    }
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
  
int addDeposit(deposit_t dep, int index, double energy[][NCELL][NCELL],
  double timeFirst[][NCELL][NCELL], double timeCent[][NCELL][NCELL])
{
  //Declare needed variables
  int nx, ny, nz;
  int keepDeposit;
  
  if (index==0 || dep.energy>E_MIN) { keepDeposit = 1; }
  else { keepDeposit = 0; }

  //Place in the detector
  nx = 1 + NCELL/2.0 + dep.x/CELL_DIM;
  ny = 1 + NCELL/2.0 + dep.y/CELL_DIM;
  nz = 1 + NCELL/2.0 + dep.z/CELL_DIM;

  //Add the deposit to the detector
  if ( !((nx<1 || nx>NCELL) || (ny<1 || ny>NCELL) || (nz<1 || nz>NCELL)) )
  {
    if (keepDeposit==1)
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
