//&&&&&&&&&&&&&
//& C Headers &
//&&&&&&&&&&&&&
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//&&&&&&&&&&&&&&&&&&&&&&&&&&
//& GNU Scientific Library &
//&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//&&&&&&&&&&&&&&&&&&
//& Custom Headers &
//&&&&&&&&&&&&&&&&&&
#include "parsing.hh"
#include "randomNumberGenerator.hh"
#include "reconLib.hh"

//&&&&&&&&&&&&&
//& Constants &
//&&&&&&&&&&&&&
#define STREAM_FLAG 1             //0->write to /dev/null
#define RECON_VERBOSE_LEVEL 1     //
#define GET_CHARGE_ITERATIONS 0   //Get the charge reconstruction iterations

//&&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Declarations &
//&&&&&&&&&&&&&&&&&&&&&&&&&

void readHeaderBlock(FILE * f);
int getMaxEntryIndex(double * a, int size);
double getLeadingEdgeTime(int nx, int ny, int nz, component_t * pmtList, double ** hitPMT,
  int pmtCnt);
void bubbleSortAscendingInt(int * a, int len);

//*************************************************************************************************

//&&&&&&&&
//& Main &
//&&&&&&&&
int main (int argc, char * argv[])
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Declare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  //TO DO: Unless something absolutly has to be on the heap change it to be on the stack
  
  //Variables needed for reading the event file
  char eventFileName[stringLen_g];                //File name of data to read
  char basisFileName[stringLen_g];                //File name of the template waveforms
  char outputFileName[stringLen_g];               //File name of the output file
  char string0[stringLen_g];                      //Strings use to reading
  char string1[stringLen_g];                      //
  double ** primaryPMT;       //STACK-IT??         //Tempate waveforms for the primary PMT
  double ** sidePMT;          //STACK-IT??         //for the side PMT
  double primaryPMT_Q[NCELL];                     //The charge for primary PMTs. These will be set
                                                  //such that the ith element is the number of PEs
                                                  //needed to make 1 PE in the PMT. See below
  double sidePMT_Q[NCELL];                        //The charge for the side PMTs. See note above
  double primaryPMT_t[NCELL];                     //The start time of the primary PMT
  double primaryPMT_tPrime[NCELL];                //
  double binCount;                                //Used for reading the PMT waveforms
  int splitIndex, ibdFlag=-1;                     //Index to determine output files
  int reconType, nCells, latticeFlag, mirrorFlag; //Configuration variables
  int guideFlag;                                  //
  int binNum, nLines, eventNum;                   //Used for reading the waveforms and event file
  int continueFlag;                               //Flag used, if testing, to control the output of
                                                  //reconstruction infomation
  pmtTimeBin_t hit;                               //Used for reading the waveforms
  
  
  //Variables needed for the reconstruction
  component_t * pmtList, * cellList;  //The list of hit pmt's that stores the charge and start time
                                      //of the pulse. The list of cells hit that stores the energy                                 
                                      //and the time
  component_t * tempList;             //
  detector_t * initDet, * saveDet;    //Detector that is used to initialize the cell list
  double ** hitPMT;                   //Stores the hit PMT's waveformn
  double conversionSlope;             //The slope of the fit that converts PEs to energy
  double conversionIntercept;         //The intercept of the fit that converts PEs to energy
  double cellThres;                   //Threshold of a cell for a cell hit
  double x0, xN, y0, yN, z0, zN;      //Variables used to set the energy and the time of the hit cells
  double setTime0;                    //
  double chiSqrInit, chiSqrFinal;     //The initial and final chi-squared, used to determine if
                                      //reconstruction was successful
  double chiSqrOld, chiSqrNew;        //Used to check the minimization process
  int firstEvent=1;                   //Used for handleing memory allocation
  int pmtCnt, cellCnt, cellCntFinal;  //Number of hit pmts and cells
  int loopCnt, index, stop;           //Number of loops of reconstruction, variable for indexing
                                      //arrays
  int nReInit;                        //
  int eventCnt=0, pmtIndex;           //Event counter to track progress, index used to access array
                                      //elements
  int nx, ny, nz;                     //
  int genFlag, seedFlag, seed;        //
  gsl_rng *  rGen;                    //      
  const gsl_rng_type * T;             //

  //Files
  FILE * config, * event, * out;  //The configuration, event, and output files
  FILE * basis, * iter;          //The template waveform and recon interaction files  
  
  
  //&&&&&&&&&&&&&
  //& Initalize &
  //&&&&&&&&&&&&&
  
  //Get CL arguments
  if (argc<2)
  {
    printf("Incorrect call. Please give the event file to reconstruct.\n");
    printf("Exiting...\n");
    exit(1);
  }
  
  initString(eventFileName, stringLen_g);
  copyString(eventFileName, argv[1], int(strlen(argv[1])) );

  if (argc==3)
  {
    splitIndex = atoi(argv[2]);
  }
  else splitIndex = -1;
  
  if (strstr(eventFileName, "IBD")==NULL) { ibdFlag=0; }
  else { ibdFlag=1; } 


  //Initialize the random number generator (may not need)
  genFlag=1;
  seedFlag=0;
    
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


  //Open the event file
  event = fopen(eventFileName, "r");
  if (event==NULL)
  {
    printf("Could not open the file %s! Exiting...\n", eventFileName);
    exit(1);
  }


  //Open the configure file
  config = fopen("inputs/reconConfig.txt", "r");
  if (config==NULL)
  {
    printf("Could not open the configuration file! Exiting...\n");
    exit(1);
  }
  
  nLines = getNumEntries(config);
  if (nLines!=5)  //Note minimal error checking here
  {
    printf("Incorrect format for the configuration! Exiting...\n");
    exit(1);
  }
  
  //Get the configuration
  
  //TO DO: Include the dependence of the PMT coupling/non-coupling
  
  initString(string0, stringLen_g);
  getEntry(config, string0, stringLen_g);
  sscanf(string0, "%s %d", string1, &reconType);
  
  initString(string0, stringLen_g);
  getEntry(config, string0, stringLen_g);
  sscanf(string0, "%s %d", string1, &nCells);
    
  initString(string0, stringLen_g);
  getEntry(config, string0, stringLen_g);
  sscanf(string0, "%s %d", string1, &latticeFlag);  

  initString(string0, stringLen_g);
  getEntry(config, string0, stringLen_g);
  sscanf(string0, "%s %d", string1, &mirrorFlag);
  
  initString(string0, stringLen_g);
  getEntry(config, string0, stringLen_g);
  sscanf(string0, "%s %d", string1, &guideFlag);

  fclose(config);

  //Check the configuration  
  if (reconType<0 || reconType>2)
  {
    printf("There is not such reconstruction option %d. Preceeding with option 0.\n", reconType);
    reconType=0;
  }
  
  if (nCells!=NCELL || latticeFlag!=LATTICE_FLAG || mirrorFlag!=MIRROR_FLAG ||
    guideFlag!=GUIDE_FLAG)
  {
    printf("There is a mismatch between the macros defined in detectorParameters.hh and in reconConfig.txt!\n");
    printf("Adjust this before running. Exiting...\n");
    exit(1);
  }
  
  //Open the output file
  initString(outputFileName, stringLen_g);
  if (splitIndex==-1)
  {
    if (ibdFlag==0)
    {
      sprintf(outputFileName, "outputs/reconstruction.out");
    }
    else
    {
      if (strstr(eventFileName, "Positron"))
      {
        sprintf(outputFileName, "outputs/reconstructionPositron.out");      
      }
      else
      {
        sprintf(outputFileName, "outputs/reconstructionNeutron.out");            
      }
    }
  }
  else
  {
    sprintf(outputFileName, "outputs/reconstruction%d.out", splitIndex);
  }
  
  out = fopen(outputFileName, "w"); 
  if (out==NULL)
  {
    printf("Could not open the output file %s! Exiting...\n", outputFileName);
    exit(1);
  }

  
  //Get the overall factor to convert photoelectrons into energy
  getEnergyScaleFactors(&conversionSlope, &conversionIntercept);
  if (TEST_FLAG_RECON==1)
  {
    printf("The energy conversion function is %lfxPEs + %lf\n", conversionSlope,
      conversionIntercept);
  }  
  

  //Get the hit cell threshold
  if (CELL_THRESHOLD==0) { cellThres = 0.0; }
  else if (CELL_THRESHOLD==1)
  {
    cellThres = conversionSlope + conversionIntercept;
    if (cellThres<0.0)
    {
      printf("The cell hit threshold was negative! Setting to %lf!\n", E_MIN);
      cellThres = E_MIN;
    }
  }
  else { cellThres = E_MIN; }  
  
  
  //If printing the iterations initialize that
  int getIterationsFlag=0;
  if (GET_CHARGE_ITERATIONS)
  {
    //TO DO: Add the information on this functionality to the README
    iter = fopen("outputs/reconIterations.out", "w");
    if (iter==NULL)
    {
      printf("Could not open the file outputs/reconIterations.out!\n");
      printf("No reconstruction iterations will be printed!\n");
    }
    else { getIterationsFlag=1; }
  }
  

  //Allocate memory for the basis data
  primaryPMT = new double* [NCELL];
  sidePMT = new double* [NCELL];
  for (int i=0; i<NCELL; i++)
  {
    primaryPMT[i] = new double [N_TIME_BIN_PRINT];
    sidePMT[i] = new double [N_TIME_BIN_PRINT];
  }
  
  //Open the basis data file
  initString(basisFileName, stringLen_g);
  sprintf(basisFileName, "reconstructionData/reconBasisData_%d_%d%d%d_%d_%d.dat", NCELL,
    LATTICE_FLAG, MIRROR_FLAG, GUIDE_FLAG, TIME_BIN_SIZE, N_TIME_BIN_PRINT);
    
  basis = fopen(basisFileName, "r");
  if (basis==NULL)
  {
    printf("Could not open the tempate waveforms! Exiting...\n");
    exit(1);
  }
  
  //Get the basis data
  for (int i=0; i<NCELL; i++)
  {
    primaryPMT_Q[i]=0.0;
    sidePMT_Q[i]=0.0;
    
    primaryPMT_t[i]=-1.0;
    primaryPMT_tPrime[i]=0.0;    
            
    for (int j=0; j<N_TIME_BIN_PRINT; j++)
    {
      sidePMT[i][j]=0.0;
      primaryPMT[i][j]=0.0;
    }
  }
  
  //Check the preamble
  int checkFlag = checkReconData(basis);
  if (checkFlag==0)
  {
    printf("There is a mismatch betweent he currently defined macros and the reconsruction data!\n");
    printf("Would you like to proceed? Enter 0 to exit 1 to continue\n");
    printf(" > ");
    scanf("%d", &continueFlag);
    
    if (continueFlag==0)
    {
      printf("To set the reconstruction data see the README. Exiting...\n");    
      exit(1);
    }
  }
  
  for (int j=0; j<NCELL; j++)  //Read the waveforms
  {
    //Read the first two lines
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, basis);

    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, basis);
    
    //Read the primary waveform
    for (int k=0; k<N_TIME_BIN_PRINT; k++)
    {
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, basis);

      sscanf(string0, "%d %lf", &binNum, &binCount);
      primaryPMT[j][binNum-1] = binCount;
    }  
    
    //Read the next line
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, basis);
    
    //Read the side waveform
    for (int k=0; k<N_TIME_BIN_PRINT; k++)
    {
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, basis);

      sscanf(string0, "%d %lf", &binNum, &binCount);
      sidePMT[j][binNum-1] = binCount;
    }          
  } //End reading the wave forms
  fclose(basis);  
   
  //Get the basis chargess  
  for (int i=0; i<NCELL; i++)
  {
    for (int j=0; j<N_TIME_BIN_PRINT; j++)
    {
      primaryPMT_Q[i] += primaryPMT[i][j];
      sidePMT_Q[i] += sidePMT[i][j];
      
      if (primaryPMT_t[i]<0.0 && primaryPMT[i][j]>0.0) primaryPMT_t[i] = j;
    }
  }

  for (int i=0; i<NCELL; i++) //Set the side charges to be relative to the charge in the primary
  {
    sidePMT_Q[i] = sidePMT_Q[i]/primaryPMT_Q[i];  
  }

  for (int i=1; i<NCELL; i++) //Set the charges needed in the ith cell needed to get the same amount
  {                           //of light as 1 cell away. Effectivly this is getting the amount of 
                              //enery in the cell required for a unit charge in the PMT
    primaryPMT_Q[i] = primaryPMT_Q[0]/primaryPMT_Q[i];
  }
  primaryPMT_Q[0] = 1.0;
  
  for (int i=0; i<NCELL; i++)
  {
    primaryPMT_tPrime[i] = getMaxEntryIndex(primaryPMT[i], N_TIME_BIN_PRINT) - primaryPMT_t[i] + 1;
    primaryPMT_tPrime[i] = primaryPMT_tPrime[i]/2.0 + primaryPMT_t[i];
  }


  //Open the default output stream
  if (STREAM_FLAG==0)
  {
    if (freopen("/dev/null", "w", stdout)==NULL) exit(1);
  }
  

  //Output the tempate waveforms if desired
  if (TEST_FLAG_RECON==1)
  {
    printf("\nTo output the primary waveforms enter 1, otherwise enter 0\n");
    printf(" > ");
    scanf("%d", &continueFlag);

    if (continueFlag)
    {
      printf("Primary PMT Waveforms:\n");    
      for (int i=0; i<NCELL; i++)
      {
        printf("  %d Cell(s) Away\n", i+1);
        for (int j=0; j<N_TIME_BIN_PRINT; j++)
        {
          printf("    %d %lf\n", j, primaryPMT[i][j]);
        }
        
        //Check if continuing
        printf("\nTo continue with the next waveform enter 1, otherwise enter 0\n");
        printf(" > ");
        scanf("%d", &continueFlag);
        
        if (continueFlag) printf("\n");
        else break;
      }
      printf("\n");
    }
    
    printf("\nTo output the side waveforms enter 1, otherwise enter 0\n");
    printf(" > ");
    scanf("%d", &continueFlag);    
    
    if (continueFlag)
    {
      printf("Side PMT Waveforms:\n");    
      for (int i=0; i<NCELL; i++)
      {
        printf("  %d Cell(s) Away\n", i+1);
        for (int j=0; j<N_TIME_BIN_PRINT; j++)
        {
          printf("    %d %lf\n", j, sidePMT[i][j]);
        }
        
        //Check if continuing
        printf("\nTo continue with the next waveform enter 1, otherwise enter 0\n");
        printf(" > ");
        scanf("%d", &continueFlag);
        
        if (continueFlag) printf("\n");
        else break;
      }
      printf("\n");
    }
    
    printf("\nTo output the charge data enter 1, otherwise enter 0\n");
    printf(" > ");
    scanf("%d", &continueFlag);    
    
    if (continueFlag)
    {
      printf("Primary PMT Charge:\n");    
      for (int i=0; i<NCELL; i++)
      {
        printf("  %d %lf\n", i+1, primaryPMT_Q[i]);
      }
      
      printf("\nSide PMT Charge:\n");    
      for (int i=0; i<NCELL; i++)
      {
        printf("  %d %lf\n", i+1, sidePMT_Q[i]);
      }
      printf("\n");
    }
  }


  
  //&&&&&&&&&&&&&&&&&&&&&&&&
  //& Read and Reconstruct &
  //&&&&&&&&&&&&&&&&&&&&&&&&    
   
  initString(string0, stringLen_g);
  while (fgets(string0, stringLen_g, event)!=NULL)
  {
    //&&&&&&&&&&&&&&&&&&&
    //& Read Event File &
    //&&&&&&&&&&&&&&&&&&&  
    
    //TO DO: Make reading the unimportant information a function so that main() is a little cleaner
  
    //Check the event header
    if (strstr(string0, "#Event")==NULL)
    {
      printf("Error at event header!\n");
      printf("  %s", string0);
      continue;
    }
            
    //Read the event header
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, event);
    sscanf(string0, "%d", &nLines);
    
    //Read the block
    for (int j=0; j<nLines; j++)
    {
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event);
      
      if (j==0) sscanf(string0, "%d", &eventNum);
    }
    
    eventCnt++;
    if (RECON_VERBOSE_LEVEL>1)
    {
      int mod = pow(10.0, RECON_VERBOSE_LEVEL);
      if ((eventCnt+1)%mod==0 && RECON_VERBOSE_LEVEL>=2)
      {
        printf("%d events reconstructed\n", eventCnt+1);
      }
    }  
    
    if (GET_CHARGE_ITERATIONS==1 && getIterationsFlag==1) { fprintf(iter, "#Event %d\n", eventNum); }

    //Check the Geant4 properties header
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, event);
    if (strstr(string0, "#Geant4_properties")==NULL)
    {
      printf("Error at Geant4 properties header!\n");
      printf("  %s", string0);
      continue;
    } 

    //Read the Geant4 properties header
    readHeaderBlock(event);

    //Check energy deposits header
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, event);
    if (strstr(string0, "#Energy_deposits")==NULL)
    {
      printf("Error at energy properties header!\n");
      printf("  %s", string0);
      continue;
    }  

    //Read the energy deposits
    readHeaderBlock(event);
    
    //Check the light transport properties header
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, event);
    if (strstr(string0, "#Light_transport_properties")==NULL)
    {
      printf("Error at light transport properties header!\n");
      printf("  %s", string0);
      continue;
    }    
    
    //Read the light transport properties
    readHeaderBlock(event);
               
    //Check the cell deposits header
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, event);
    if (strstr(string0, "#Cells_deposits")==NULL)
    {
      printf("Error at cell deposits header!\n");
      printf("  %s", string0);
      continue;
    }       

    //Read the cell deposits properties
    readHeaderBlock(event);

    //Check the PMT charges header
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, event);
    if (strstr(string0, "#PMT_charges")==NULL)
    {
      printf("Error at PMT charges header!\n");
      printf("  %s", string0);
      continue;
    }       
    
    //Read the hit PMTs and their charges
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, event);
    sscanf(string0, "%d", &nLines);
    
    if (nLines==0)
    {
      printf("No light transport data to reconstruct for event %d!\n", eventNum);
      fprintf(out, "Event %d\n0\n", eventNum); 
      
      initString(string0, stringLen_g);   //Read the next two lines to get to the next event header
      fgets(string0, stringLen_g, event);
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event);      
      continue;     
    }
    
    pmtCnt = nLines;
    pmtList = new component_t [pmtCnt];
    clearComponents(pmtList, pmtCnt);
    for (int j=0; j<pmtCnt; j++)
    {
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event);
      
      sscanf(string0, "%d %d %d %lf", &pmtList[j].nx, &pmtList[j].ny, &pmtList[j].nz, 
        &pmtList[j].e);
      pmtList[j].t = INFINITY;
    } //End reading the hit pmts
    
    //Check the PMT waveforms header
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, event);
    if (strstr(string0, "#PMT_waveforms")==NULL)
    {
      printf("Error at PMT waveforms header!\n");
      printf("  %s", string0);
      continue;
    }       

    //Allocate the PMT waveforms
    hitPMT = new double* [pmtCnt];
    for (int j=0; j<pmtCnt; j++)
    {
      hitPMT[j] = new double [N_TIME_BIN_PRINT];
    } 

    for (int j=0; j<pmtCnt; j++)
    {
      for (int k=0; k<N_TIME_BIN_PRINT; k++)
      {
        hitPMT[j][k]=0.0;
      }
    }   
    
    //Read the PMT waveforms and their charges
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, event);
    sscanf(string0, "%d", &nLines);
    
    if (nLines==0)
    {
      printf("Error in event file formatting! Exiting...\n");
      exit(1);
    }
    
    for (int j=0; j<nLines; j++)
    {
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event);
      
      sscanf(string0, "%d %d %d %d %lf", &hit.nu, &hit.nv, &hit.nw, &hit.bin, &hit.count);
      
      for (int k=0; k<pmtCnt; k++)
      {
        if (hit.nu==pmtList[k].nx && hit.nv==pmtList[k].ny && hit.nw==pmtList[k].nz)
        {
          hitPMT[k][hit.bin-1] = hit.count; //Need the minus one here to convert from Fortran's
                                            //1-index arrays to C's 0-indexed arrays
                                            
          if (pmtList[k].t>(hit.bin-1)) pmtList[k].t = hit.bin-1;
          break;
        }
      }
    } //End reading the hits    


    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //& Initialize the Hit Cells &
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     
    //Add energy to the cells:
    //  Loop over all the cells and see if there are three orthogonal PMTs with light in them. If
    //  this is the case then add photoelectrons to the cell. The number of photoelectrons is
    //  computed by assuming that all of the charge in the hit PMTs comes only from the cell. The
    //  minimum of this set is taken as the initial value. 
      
    if (firstEvent==1)
    {
      initDet = new detector_t;
      saveDet = new detector_t;
      firstEvent=0;
    }
    clearDetector(initDet);
        
    for (int i=0; i<pmtCnt; i++)
    {
      initDet->component[pmtList[i].nx][pmtList[i].ny][pmtList[i].nz] = pmtList[i].e;
    }
    
    //Loop over the cells
    for (int i=1; i<=NCELL; i++) {
      for (int j=1; j<=NCELL; j++) {
        for (int k=1; k<=NCELL; k++) {

          //Check if three orthogonal PMTs have light move to next cell if this is not the case
          index=0;
          x0 = initDet->component[0][j][k];
          y0 = initDet->component[i][0][k];
          z0 = initDet->component[i][j][0];
          
          if (MIRROR_FLAG==0) //Use ALL of the PMTs
          {
            xN = initDet->component[NCELL+1][j][k];          
            yN = initDet->component[i][NCELL+1][k];          
            zN = initDet->component[i][j][NCELL+1];  
                    
            if ( (x0>0.0 || xN>0.0) && (y0>0.0 || yN>0.0) && (z0>0.0 || zN>0.0) )
            {
              index=1;
            }
          }
          else                //Use only the "0-side" PMTs
          {
            if ( x0>0.0 && y0>0.0 && z0>0.0 )
            {
              index=1;
            }        
          }
          if (index==0) { continue; } //Move to the next cell
          
        
          //Determine the charges from each of the PMTs
          
          //  Reconstructing events in a detector without light guides? READ THIS
          //
          //  If not using light guides should I adjust the initialization to take into account the
          //  non-uniform light collection from the cells in the outer layer?
          //    
          //  For now (2015-06-19) I will not change the initialization since this won't effect
          //  things much.
          
          x0 = initDet->component[0][j][k]*primaryPMT_Q[i-1];
          y0 = initDet->component[i][0][k]*primaryPMT_Q[j-1];
          z0 = initDet->component[i][j][0]*primaryPMT_Q[k-1];          
          
          if (MIRROR_FLAG==0)
          {
            xN = initDet->component[NCELL+1][j][k]*primaryPMT_Q[NCELL-i];
            yN = initDet->component[i][NCELL+1][k]*primaryPMT_Q[NCELL-j];
            zN = initDet->component[i][j][NCELL+1]*primaryPMT_Q[NCELL-k];
            
            if (x0>0.0 && xN>0.0) x0 = fmin(x0, xN);
            else if (x0<=0.0 && xN>0.0) x0 = xN;

            if (y0>0.0 && yN>0.0) y0 = fmin(y0, yN);
            else if (y0<=0.0 && yN>0.0) y0 = yN;

            if (z0>0.0 && zN>0.0) z0 = fmin(z0, zN);
            else if (z0<=0.0 && zN>0.0) z0 = zN;
          }
          
          initDet->component[i][j][k] = fmin(x0, fmin(y0, z0));
    } } } //End loop over cells

    //Allocate and assign the cell list
    cellCnt=0;
    for (int i=1; i<=NCELL; i++) {
      for (int j=1; j<=NCELL; j++) {
        for (int k=1; k<=NCELL; k++) {
          if (initDet->component[i][j][k]>0.0) cellCnt++;
    } } }
    
    cellList = new component_t [cellCnt];
    clearComponents(cellList, cellCnt);  
 
    cellCnt=0;
    for (int i=1; i<=NCELL; i++) {
      for (int j=1; j<=NCELL; j++) {
        for (int k=1; k<=NCELL; k++) {
          if (initDet->component[i][j][k]>0.0)
          {
            cellList[cellCnt].nx = i;
            cellList[cellCnt].ny = j;
            cellList[cellCnt].nz = k;
            cellList[cellCnt].e = initDet->component[i][j][k];
            //printf("%d %d %d %lf\n", i, j, k, getEnergy(cellList[cellCnt].e, conversionSlope, conversionIntercept) );
            cellCnt++;
          }
    } } } 
    
    //printf("\n");
    
    //Add time to the cells
    //  NOTE: Times are in units of the bin size and only on outputing are they converted to ns. 
    //
    //  Loop over all the cells and if there is energy in them calculate the first non-zero time bin
    //  of the waveforms of the PMTs that view the cell. Next calculate the start time for the cell
    //  in question by taking the PMT's start time and subtracting the start time of the basis
    //  waveform. For the start time of the cell take the average of the times.

    //TO DO: Think of an implement some different initialization procedures and decide which is best

    //Add the pmt start times to the detector 
    copyDetector(initDet, saveDet);
    clearDetector(initDet);
    
    stop = (1-MIRROR_FLAG)*(NCELL+1);
    for (int i=0; i<=stop; i+=NCELL+1) {
      for (int j=1; j<=NCELL; j++) {
        for (int k=1; k<=NCELL; k++) {
          initDet->component[i][j][k] = N_TIME_BIN_PRINT; //x-side PMTs
          initDet->component[j][i][k] = N_TIME_BIN_PRINT; //y-side PMTs
          initDet->component[j][k][i] = N_TIME_BIN_PRINT; //z-side PMTs
    } } }     
    
    for (int i=0; i<pmtCnt; i++)
    {
      initDet->component[pmtList[i].nx][pmtList[i].ny][pmtList[i].nz] = pmtList[i].t;
    }

    //Loop over the cells and determine their start times
    for (int i=1; i<=NCELL; i++) {
      for (int j=1; j<=NCELL; j++) {
        for (int k=1; k<=NCELL; k++) {
          if (saveDet->component[i][j][k]>0.0)
          {
            //Get the number of PMTs above threshold
            int numberPmtsAboveThreshold = 0;
            int proceedWithTimeInitialization = 0;
            if (getComponentIndex(pmtList, pmtCnt, 0, j, k)>=0) { numberPmtsAboveThreshold++; }
            if (getComponentIndex(pmtList, pmtCnt, i, 0, k)>=0) { numberPmtsAboveThreshold++; }
            if (getComponentIndex(pmtList, pmtCnt, i, j, 0)>=0) { numberPmtsAboveThreshold++; }
            
            if (MIRROR_FLAG==0)
            {
              if (getComponentIndex(pmtList, pmtCnt, NCELL+1, j, k)>=0) { numberPmtsAboveThreshold++; }
              if (getComponentIndex(pmtList, pmtCnt, i, NCELL+1, k)>=0) { numberPmtsAboveThreshold++; }
              if (getComponentIndex(pmtList, pmtCnt, i, j, NCELL+1)>=0) { numberPmtsAboveThreshold++; }            
            }
            
            if (MIRROR_FLAG==0 && numberPmtsAboveThreshold==6)
            {
              proceedWithTimeInitialization = 1;
            }
            else if (MIRROR_FLAG==0 && numberPmtsAboveThreshold==3)
            {
              proceedWithTimeInitialization = 1;
            }
            
            if (proceedWithTimeInitialization)
            {
              x0 = initDet->component[0][j][k] - primaryPMT_t[i-1];
              y0 = initDet->component[i][0][k] - primaryPMT_t[j-1];
              z0 = initDet->component[i][j][0] - primaryPMT_t[k-1];
              
              if (MIRROR_FLAG==0)
              {
                xN = initDet->component[NCELL+1][j][k] - primaryPMT_t[NCELL-i];
                yN = initDet->component[i][NCELL+1][k] - primaryPMT_t[NCELL-j];
                zN = initDet->component[i][j][NCELL+1] - primaryPMT_t[NCELL-k];

                x0 = (x0 + xN)/2.0;
                y0 = (y0 + yN)/2.0;
                z0 = (z0 + zN)/2.0;
              }
                          
              initDet->component[i][j][k] = (x0 + y0 + z0)/3.0;
            } //End check if all the PMTs that view the cell have light in them      
          } //End check if the cell has energy in it. 
    } } }     

    if (TIME_RECON)
    {
      setTime0 = INFINITY;  //If the earliest time is negative then shift all start times so that
                            //earliest time is at t=0.
      for (int i=1; i<=NCELL; i++) {
        for (int j=1; j<=NCELL; j++) {
          for (int k=1; k<=NCELL; k++) {
            if (saveDet->component[i][j][k]>0.0)
            {
              if (initDet->component[i][j][k]<setTime0) setTime0 = initDet->component[i][j][k];
            }
      } } }
    }
    else setTime0 = 0.0;

    cellCnt=0;
    for (int i=1; i<=NCELL; i++) {
      for (int j=1; j<=NCELL; j++) {
        for (int k=1; k<=NCELL; k++) {
          if (saveDet->component[i][j][k]>0.0)
          {
            if ( cellList[cellCnt].nx==i && cellList[cellCnt].ny==j && cellList[cellCnt].nz==k )
            {
              if (setTime0<0.0) cellList[cellCnt].t = initDet->component[i][j][k] - setTime0;
              else cellList[cellCnt].t = initDet->component[i][j][k];
              
              cellCnt++;
            }
            else
            {
              printf("Error in cell time assignment!\n");
            }
          }
    } } }

    
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //& Re-Initialize the PMT List &
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    //  NOTE: When reconstructing events in a detector that does not have light guides the cells  
    //  in the outter most layer have collection efficiencies that are dependent on the exact
    //  location of the deposit; this is due to the fact the the photocathode only covers some
    //  percentage of the cell's face (for 2.5" cubes and a 46mm diameter photocathode the coverage
    //  is only 41%).
    //
    //  To account for this effect PMTs that directly view (ie. are one cell away) hit cells are
    //  given a photoelectron charge of -1. Then, in the chi-squared function, this negative charge 
    //  serves as a flag that indicates the PMT should not be used in the calculation of the
    //  chi-squared. The determination of which PMTs to remove is done below
    
    if (GUIDE_FLAG==0)
    {
      //Get the indices of the PMTs to remove from the pmt list and set their charges to -1
      int removeIndex;
      for (int cellIndex=0; cellIndex<cellCnt; cellIndex++)
      {
        //x-sides
        if (cellList[cellIndex].nx==1)
        {
          removeIndex = getComponentIndex(pmtList, pmtCnt, 0, cellList[cellIndex].ny,
            cellList[cellIndex].nz);
          if (removeIndex!=-1) pmtList[removeIndex].e = -1.0;
        }
        else if (cellList[cellIndex].nx==NCELL)
        {
          removeIndex = getComponentIndex(pmtList, pmtCnt, NCELL+1, cellList[cellIndex].ny,
            cellList[cellIndex].nz);
          if (removeIndex!=-1) pmtList[removeIndex].e = -1.0;
        }
        
        //y-sides
        if (cellList[cellIndex].ny==1)
        {
          removeIndex = getComponentIndex(pmtList, pmtCnt, cellList[cellIndex].nx, 0,
            cellList[cellIndex].nz);
          if (removeIndex!=-1) pmtList[removeIndex].e = -1.0;
        }
        else if (cellList[cellIndex].ny==NCELL)
        {
          removeIndex = getComponentIndex(pmtList, pmtCnt, cellList[cellIndex].nx, NCELL+1,
            cellList[cellIndex].nz);
          if (removeIndex!=-1) pmtList[removeIndex].e = -1.0;
        }
        
        //z-sides
        if (cellList[cellIndex].nz==1)
        {
          removeIndex = getComponentIndex(pmtList, pmtCnt, cellList[cellIndex].nx,
            cellList[cellIndex].ny, 0);
          if (removeIndex!=-1) pmtList[removeIndex].e = -1.0;
        }
        else if (cellList[cellIndex].nz==NCELL)
        {
          removeIndex = getComponentIndex(pmtList, pmtCnt, cellList[cellIndex].nx, 
            cellList[cellIndex].ny, NCELL+1);
          if (removeIndex!=-1) pmtList[removeIndex].e = -1.0;
        }
      } //End loop over cells
      
    } //End re-initialization of PMT list if not using light guides
    
    //Output the pmt list if testing   
    if (TEST_FLAG_RECON==1)
    {
      printf("PMTs:\n");
      for (int i=0; i<pmtCnt; i++)
      {
        printf("%d %d %d %lf %lf\n", pmtList[i].nx, pmtList[i].ny, pmtList[i].nz,
          pmtList[i].t, pmtList[i].e);
      }
      printf("\n");

      printf("\nTo print the PMT waveforms enter 1, otherwise enter 0\n");
      printf(" > ");
      scanf("%d", &continueFlag);
      
      if (continueFlag)
      {
        for (int i=0; i<pmtCnt; i++)
        {
          for (int j=0; j<N_TIME_BIN_PRINT; j++)
          {
            if (hitPMT[i][j]>0.0)
            {
              printf("%d %d %d %d %lf\n", pmtList[i].nx, pmtList[i].ny, pmtList[i].nz, j,
                hitPMT[i][j]);
            }
          }
          
          //Check if continuing
          printf("\nTo continue with the next waveform enter 1, otherwise enter 0\n");
          printf(" > ");
          scanf("%d", &continueFlag);
        
          if (continueFlag) printf("\n");
          else break;           
        }
      }
    }
    
 
    //&&&&&&&&&&&&&&&
    //& Reconstruct &
    //&&&&&&&&&&&&&&&
    
    if (getTotalCharge(pmtList, pmtCnt)==-1.0*pmtCnt)  //If you removed all PMTs then there is nothing to do!
    {
      printf("No PMTs left to reconstruct for event %d!\n", eventNum);
      fprintf(out, "Event %d\n0\n", eventNum);
      continue;
    }
    
    if (reconType%2==0)
    {
      chiSqrInit = getChiSqrCharge(primaryPMT_Q, sidePMT_Q, pmtList, cellList, pmtCnt, cellCnt);
    }
    else
    {
      chiSqrInit = getChiSqrWaveform(primaryPMT, sidePMT, pmtList, hitPMT, cellList, pmtCnt,
        cellCnt);
    }

    //Initialize the reconstuction
    chiSqrOld = chiSqrInit;
    chiSqrNew = 0.0;
    loopCnt=0;
    
    //Reconstruct
    while (chiSqrNew<chiSqrOld)
    {
      if (TEST_FLAG_RECON==1)
      {
        printf("loopCnt %d\n", loopCnt);
        printf("Cells:\n");
        for (int i=0; i<cellCnt; i++)
        {
          printf("%d %d %d %lf %lf\n", cellList[i].nx, cellList[i].ny, cellList[i].nz,
            cellList[i].t, getEnergy(cellList[i].e, conversionSlope, conversionIntercept) );
        }
        printf("\n");
      }  
      
      if (getIterationsFlag)
      {
        fprintf(iter, "%d\n", cellCnt);
        for (int i=0; i<cellCnt; i++)
        {
          fprintf(iter, "%d %d %d %lf %lf\n", cellList[i].nx, cellList[i].ny, cellList[i].nz,
            cellList[i].t, getEnergy(cellList[i].e, conversionSlope, conversionIntercept) );
        }
      }
    
      //Switch the chi-squareds
      if (loopCnt>0)
      {
        chiSqrOld = chiSqrNew;
        chiSqrNew = 0.0;
      }
      
      //Call the gradient search. Always start with a E minimization
      gradientSearchCharge(primaryPMT_Q, sidePMT_Q, pmtList, cellList, pmtCnt, cellCnt);
      
      if (reconType%2==0)     //Charge only -> get the new chi-squared
      {
        chiSqrNew = getChiSqrCharge(primaryPMT_Q, sidePMT_Q, pmtList, cellList, pmtCnt, cellCnt);
      }
      else if (reconType==1)  //Charge and time -> step in t and get the new chi-squared
      {
        gradientSearchWaveform(primaryPMT, sidePMT, pmtList, hitPMT, cellList, pmtCnt, cellCnt);        
        chiSqrNew = getChiSqrWaveform(primaryPMT, sidePMT, pmtList, hitPMT, cellList, pmtCnt,
          cellCnt);
      }
              
      //Check the chi-squared
      if (TEST_FLAG_RECON==1) printf("new=%lf old=%lf\n", chiSqrNew, chiSqrOld);
          
      loopCnt++;    
    } //End reconstruction loop


    //Re-initialize the cell list with only cells above threshold
    cellCntFinal=0;
    for (int i=0; i<cellCnt; i++)
    {
      if (getEnergy(cellList[i].e, conversionSlope, conversionIntercept)>cellThres) cellCntFinal++;
    }

    if (cellCnt!=cellCntFinal)
    {
      tempList = new component_t [cellCntFinal];
      clearComponents(tempList, cellCntFinal);
      
      cellCntFinal=0;
      for (int i=0; i<cellCnt; i++)
      {
        if (getEnergy(cellList[i].e, conversionSlope, conversionIntercept)>cellThres)
        {
          tempList[cellCntFinal].nx = cellList[i].nx;
          tempList[cellCntFinal].ny = cellList[i].ny;
          tempList[cellCntFinal].nz = cellList[i].nz;
          tempList[cellCntFinal].t = cellList[i].t;
          tempList[cellCntFinal].e = cellList[i].e;
                  
          cellCntFinal++;
        }
      }     
    
      delete [] cellList;
      cellCnt = cellCntFinal;
      cellList = new component_t [cellCnt];
      clearComponents(cellList, cellCnt);

      for (int i=0; i<cellCnt; i++)
      {
        cellList[i].nx = tempList[i].nx;
        cellList[i].ny = tempList[i].ny;
        cellList[i].nz = tempList[i].nz;
        cellList[i].t = tempList[i].t;
        cellList[i].e = tempList[i].e;
      }

      delete [] tempList;
    }

    //Do the time minimization if needed
    if (reconType==2)
    {
      //Re-estimate the start times of the hit cells
      if (TIME_INIT)
      {
        for (int cellIndex=0; cellIndex<cellCnt; cellIndex++)
        {
          nx = cellList[cellIndex].nx;
          ny = cellList[cellIndex].ny;
          nz = cellList[cellIndex].nz;

          //Use pulse peak for guessing the time
          if (TIME_INIT==2)
          {
            //Use the leading edge for guessing the time
            x0 = getLeadingEdgeTime(0, ny, nz, pmtList, hitPMT, pmtCnt) - primaryPMT_tPrime[nx-1];
            y0 = getLeadingEdgeTime(nx, 0, nz, pmtList, hitPMT, pmtCnt) - primaryPMT_tPrime[ny-1];
            z0 = getLeadingEdgeTime(nx, ny, 0, pmtList, hitPMT, pmtCnt) - primaryPMT_tPrime[nz-1];
            
            if (MIRROR_FLAG==0)
            {
              xN = getLeadingEdgeTime(NCELL+1, ny, nz, pmtList, hitPMT, pmtCnt) -
                primaryPMT_tPrime[NCELL-nx];
              yN = getLeadingEdgeTime(nx, NCELL+1, nz, pmtList, hitPMT, pmtCnt) -
                primaryPMT_tPrime[NCELL-ny];
              zN = getLeadingEdgeTime(nx, ny, NCELL+1, pmtList, hitPMT, pmtCnt) -
                primaryPMT_tPrime[NCELL-nz]; 
              
              x0 = (x0 + xN)/2.0;       
              y0 = (y0 + yN)/2.0;       
              z0 = (z0 + zN)/2.0;       
            }
          } //End use 50% leading edge time
          else
          {
            int pmtIndex0 = getComponentIndex(pmtList, pmtCnt, 0, ny, nz);
            x0 = getMaxEntryIndex(hitPMT[pmtIndex0], N_TIME_BIN_PRINT) -
              getMaxEntryIndex(primaryPMT[nx-1], N_TIME_BIN_PRINT);

            pmtIndex0 = getComponentIndex(pmtList, pmtCnt, nx, 0, nz);
            y0 = getMaxEntryIndex(hitPMT[pmtIndex0], N_TIME_BIN_PRINT)-
              getMaxEntryIndex(primaryPMT[ny-1], N_TIME_BIN_PRINT);

            pmtIndex0 = getComponentIndex(pmtList, pmtCnt, nx, ny, 0);
            z0 = getMaxEntryIndex(hitPMT[pmtIndex0], N_TIME_BIN_PRINT) -
              getMaxEntryIndex(primaryPMT[nz-1], N_TIME_BIN_PRINT);

            if (MIRROR_FLAG==0)
            {
              pmtIndex0 = getComponentIndex(pmtList, pmtCnt, NCELL+1, ny, nz);
              xN = getMaxEntryIndex(hitPMT[pmtIndex0], N_TIME_BIN_PRINT) -
                getMaxEntryIndex(primaryPMT[NCELL-nx], N_TIME_BIN_PRINT);

              pmtIndex0 = getComponentIndex(pmtList, pmtCnt, nx, NCELL+1, nz);
              yN = getMaxEntryIndex(hitPMT[pmtIndex0], N_TIME_BIN_PRINT)-
                getMaxEntryIndex(primaryPMT[NCELL-ny], N_TIME_BIN_PRINT);

              pmtIndex0 = getComponentIndex(pmtList, pmtCnt, nx, ny, NCELL+1);
              zN = getMaxEntryIndex(hitPMT[pmtIndex0], N_TIME_BIN_PRINT) -
                getMaxEntryIndex(primaryPMT[NCELL-nz], N_TIME_BIN_PRINT);
        
              x0 = (x0 + xN)/2.0;       
              y0 = (y0 + yN)/2.0;       
              z0 = (z0 + zN)/2.0;        
            }
          } //End use max bin

          cellList[cellIndex].t = (x0 + y0 + z0)/3.0; 
          //printf("%d %d %d %lf %lf \n", nx, ny, nz, cellList[cellIndex].t, cellList[cellIndex].t*0.37);
        } //End cell hit time re-initialization loop 
      } //End check if reinitializing 
      
      
    
      //Get required info for testing statements
      if (TEST_FLAG_RECON==2)
      {
        pmtIndex = getHighEnergyCellsPmtIndex(pmtList, cellList, pmtCnt, cellCnt);
      } //End testing statements
    
      //Initialize the reconstuction
      chiSqrOld = getChiSqrWaveform(primaryPMT, sidePMT, pmtList, hitPMT, cellList, pmtCnt, cellCnt);
      chiSqrNew = 0.0;
      loopCnt=0;
      nReInit=0;
      
      //Reconstruct
      while (chiSqrNew<chiSqrOld)
      {
        //Test statements
        if (TEST_FLAG_RECON==1) printf("loopCnt %d\n", loopCnt);
 
        if (TEST_FLAG_RECON==2 && loopCnt==0) //No need in continuing to print this out
        {
          //Print the highest cell waveform
          printf("#Event %d\n", eventNum);
          printf("PMT Waveform:\n");
          for (int bin=0; bin<N_TIME_BIN_PRINT; bin++)
          {
            printf("%d %lf\n", bin, hitPMT[pmtIndex][bin]);
          }
          printf("\n");
                    
          //Note that the waveform for the reconstructed PMT will be printed from the chi-squared
          //function so that it does not have to be calculated seperatly here.           
        } //End test statement block
      
       
        //Switch the chi-squareds
        if (loopCnt>0)
        {
          chiSqrOld = chiSqrNew;
          chiSqrNew = 0.0;
        }
        
        //Call the gradient search.
        gradientSearchWaveform(primaryPMT, sidePMT, pmtList, hitPMT, cellList, pmtCnt, cellCnt);        
        chiSqrNew = getChiSqrWaveform(primaryPMT, sidePMT, pmtList, hitPMT, cellList, pmtCnt,
          cellCnt);
                
        //Check the chi-squared
        if (TEST_FLAG_RECON==1)
        {
          printf("new=%lf old=%lf\n", chiSqrNew, chiSqrOld);

          printf("Cells:\n");
          for (int i=0; i<cellCnt; i++)
          {
            printf("%d %d %d %lf %lf\n", cellList[i].nx, cellList[i].ny, cellList[i].nz,
              cellList[i].t, cellList[i].e);
          }
          printf("\n");          
        }
        
        //If the chi-squared didn't change on the first pass re-initialize the cell times
        if (TIME_RECON==1 && chiSqrOld==chiSqrNew && loopCnt==0 && nReInit<5) //Note that the 5 here 
        {                                                                  //was choosen to at random
          for (int i=0; i<cellCnt; i++)
          {
            cellList[i].t += gsl_ran_gaussian(rGen, 0.2); //The best time resolution is probably on
          }                                               //the order of 1/10th of the sampling
                                                          //period; hence the choice here.
          chiSqrOld = getChiSqrWaveform(primaryPMT, sidePMT, pmtList, hitPMT, cellList, pmtCnt,
            cellCnt);
          chiSqrNew=0.0;                                                
          nReInit++;
          loopCnt--;
        }
            
        loopCnt++;    
      } //End reconstruction loop
    } //End doing the time reconstruction only
    
    
    //Check the minizmation
    if (reconType==0)
    {
      chiSqrFinal = getChiSqrCharge(primaryPMT_Q, sidePMT_Q, pmtList, cellList, pmtCnt, cellCnt);
    }
    else
    {
      chiSqrFinal = getChiSqrWaveform(primaryPMT, sidePMT, pmtList, hitPMT, cellList, pmtCnt,
        cellCnt);
    }
    
    if (TEST_FLAG_RECON==1) printf("chiSqrInit = %lf, chiSqrFinal = %lf\n\n", chiSqrInit, chiSqrFinal);
    
    if (chiSqrInit<=chiSqrFinal && TEST_FLAG_RECON==1)
    {
      printf("The minimumization failed!\n\n");
    }


    //&&&&&&&&&&&
    //& Outputs &
    //&&&&&&&&&&&
    
    fprintf(out, "Event %d\n%d\n", eventNum, cellCnt);
    for (int i=0; i<cellCnt; i++)
    {
      fprintf(out, "%d %d %d %lf %lf\n", cellList[i].nx, cellList[i].ny, cellList[i].nz,
        (cellList[i].t*TIME_BIN_SIZE)/1000.0,
          getEnergy(cellList[i].e, conversionSlope, conversionIntercept) );
    }


    //&&&&&&&&&&&&
    //& Clean-up &
    //&&&&&&&&&&&&

    delete [] pmtList;
    for (int i=0; i<pmtCnt; i++) { delete [] hitPMT[i]; }
    delete [] hitPMT;
    delete [] cellList;
    
    if (TEST_FLAG_RECON)
    {
      printf("Enter 1 to continue, 0 to exit\n");
      printf(" > ");
      scanf("%d", &continueFlag);
      
      if (!continueFlag) break;
    }

  } //End reading the event file
  
  if (RECON_VERBOSE_LEVEL) printf("%s reconstructed!\n", eventFileName);

  //&&&&&&&
  //& End &
  //&&&&&&&

  //De-allocate memory
  if (genFlag!=0) gsl_rng_free(rGen);
  
  for (int i=0; i<NCELL; i++)
  {
    delete [] primaryPMT[i];
    delete [] sidePMT[i];
  }  
  delete [] primaryPMT;
  delete [] sidePMT;
  
  delete initDet;
  delete saveDet;
      
  fclose(event);
  fclose(out);
  if (getIterationsFlag) { fclose(iter); }
    
  return 0;
}

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&

//*************************************************************************************************

void readHeaderBlock(FILE * f)
{
  //Declare needed varialbes
  char s0[stringLen_g];
  int nLines;
  
  //Initialize
  initString(s0, stringLen_g);
  fgets(s0, stringLen_g, f);
  sscanf(s0, "%d", &nLines);
  
  //Read the block
  for (int j=0; j<nLines; j++)
  {
    initString(s0, stringLen_g);
    fgets(s0, stringLen_g, f);
  }  
}

//*************************************************************************************************

int getMaxEntryIndex(double * a, int size)
{
  int index=-1;
  double check=0.0;
  
  for (int i=0; i<size; i++)
  {
    if (a[i]>check)
    {
      check = a[i];
      index = i;
    }
  }
  
  return index;
}

//*************************************************************************************************

double getLeadingEdgeTime(int nx, int ny, int nz, component_t * pmtList, double ** hitPMT,
  int pmtCnt)
{
  //Get the pmt index and ensure that it is in the pmt list
  int pmtIndex = getComponentIndex(pmtList, pmtCnt, nx, ny, nz);
  if (pmtIndex==-1) return -1.0;
      
  //Get the first non-zero bin
  int n1;
  for (int binNum=0; binNum<N_TIME_BIN_PRINT; binNum++)
  {
    if (hitPMT[pmtIndex][binNum]>0.0)
    {
      n1 = binNum;
      break;
    }
  }
        
  //Get the max bin
  int n2 = getMaxEntryIndex(hitPMT[pmtIndex], N_TIME_BIN_PRINT);
  
  //Get the "leading-edge line"
  double c1 = 0.0;
  double c2 = hitPMT[pmtIndex][n2];
  n2++;
  
  double m = (c2 - c1)/double(n2 - n1);
  double b = -m*n1;
  
  //Return the leading edge time
  return (c2/2.0 - b)/m;
}

//*************************************************************************************************

void bubbleSortAscendingInt(int * a, int len)
{
 //Declare needed variables
  int temp;
  int swapCnt=1;
  
  //Bubble sort the array
  while (swapCnt)
  {
    swapCnt=0;
    for (int i=0; i<len-1; i++)
    {
      if (a[i]>a[i+1])
      {
        temp = a[i];
        a[i] = a[i+1];
        a[i+1] = temp;

        swapCnt++;
      } //End ordering check
    } //End loop over elments in the array
  } //End while loop
} //End function

//*************************************************************************************************
