//&&&&&&&&&&&&&
//& C Headers &
//&&&&&&&&&&&&&
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//&&&&&&&&&&&&&&&&&&
//& Custom Headers &
//&&&&&&&&&&&&&&&&&&
#include "parsing.hh"
#include "randomNumberGenerator.hh"
#include "reconLib.hh"

//&&&&&&&&&&&&&
//& Constants &
//&&&&&&&&&&&&&
#define STREAM_FLAG 1           //0->write to /dev/null
#define RECON_VERBOSE_LEVEL 0   //
#define GET_CHARGE_ITERATIONS 0 //Get the charge reconstruction iterations
#define N_HIT_PMTS_REQUIRED 6   //
#define HIT_PMT_THRESHOLD 3     //

//&&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Declarations &
//&&&&&&&&&&&&&&&&&&&&&&&&&

void readHeaderBlock(FILE * f);

//*************************************************************************************************

//&&&&&&&&
//& Main &
//&&&&&&&&
int main (int argc, char * argv[])
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Declare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  //Program control and initialization variables 
  char eventFileName[stringLen_g];    //File name of data to read
  char outputFileName[stringLen_g];   //File name of the output file    
  char basisFileName[stringLen_g];    //File name of the template waveforms
  char string0[stringLen_g];          //Strings use to reading
  char string1[stringLen_g];          //
  double binCount;                    //Used for reading the PMT waveforms  
  int reconType, nCells, latticeFlag; //Configuration variables
  int mirrorFlag, guideFlag;          //
  int splitIndex;                     //Index to determine output files if running in paralell
  int getIterationsFlag;              //Flag set by the program to determine if the minimization
                                      //steps are outputted or not.
  int continueFlag;                   //
  int nSidePMTs;                      //
  int binNum, nLines, eventNum;       //Used for reading the waveforms and event file 
  int eventCnt;
  int firstEvent;                     //Used for handleing memory allocation
      
  //Variables needed for the reconstruction
  component_t * pmtList, * cellList;  //The list of hit pmt's that stores the charge and start time
                                      //of the pulse. The list of cells hit that stores the energy                                 
                                      //and the time
  component_t * tempList;             //
  detector_t * initDet;               //Detector that is used to initialize the cell list
  double conversionSlope;             //The slope of the fit that converts PEs to energy
  double conversionIntercept;         //The intercept of the fit that converts PEs to energy
  double cellThres;                   //Threshold of a cell for a cell hit
  double * primaryCharge;             //The charge for primary PMTs. These will be set such that the
                                      //ith element is the number of PEs needed to make 1 PE in the
                                      //PMT. See below
  double ** sideCharge;               //The charge for the side PMTs. These will be set relative to
                                      //their primary's PMT charge
  double x0, xN, y0, yN, z0, zN;      //Variables used to set the energy and the time of the hit cells
  double chiSqrInit, chiSqrFinal;     //The initial and final chi-squared, used to determine if
                                      //reconstruction was successful
  double chiSqrOld, chiSqrNew;        //Used to check the minimization process
  int pmtCnt, cellCnt, cellCntFinal;  //Number of hit pmts and cells
  int index;                          //Variables for indexing arrays and loops
  int loopCnt;                        //Number of loops of reconstruction
                                           
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
  

  //Open the event file
  event = fopen(eventFileName, "r");
  if (event==NULL)
  {
    printf("Could not open the file %s! Exiting...\n", eventFileName);
    exit(1);
  }


  //Open the configure file
  config = fopen("inputs/reconTeflonConfig.txt", "r");
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
  if (reconType!=0) {reconType = 0; }
  
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
    sprintf(outputFileName, "outputs/reconstruction.out");
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
  getIterationsFlag=0;
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
  

  //Allocate memory for the basis data and initialize
  nSidePMTs = NCELL - 1;
    
  primaryCharge = new double [NCELL];
  sideCharge = new double* [NCELL];
  for (int i=0; i<NCELL; i++) { sideCharge[i] = new double [nSidePMTs]; }
  
  for (int i=0; i<NCELL; i++)
  {
    primaryCharge[i] = 0.0;
    
    for (int j=0; j<nSidePMTs; j++) { sideCharge[i][j] = 0.0; }
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
  
  //Check the basis data preamble
  if (checkReconData(basis)==0)
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
  
  //Get the basis data
  for (int i=0; i<NCELL; i++)
  {
    //Read the first two lines
    for (int j=0; j<2; j++)
    {
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, basis);
    }
    
    //Read the primary waveform
    for (int j=0; j<N_TIME_BIN_PRINT; j++)
    {
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, basis);
      sscanf(string0, "%d %lf", &binNum, &binCount);

      primaryCharge[i] += binCount;
    }
    
    //Read the side waveforms
    for (int j=0; j<nSidePMTs; j++)
    {
      //Read the first line
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, basis);      
    
      for (int k=0; k<N_TIME_BIN_PRINT; k++)
      {
        initString(string0, stringLen_g);
        fgets(string0, stringLen_g, basis);
        sscanf(string0, "%d %lf", &binNum, &binCount);

        sideCharge[i][j] += binCount;
      }
    }      
  } //End reading and writing the wave forms
  fclose(basis);  
  
  //Get the non-uniformity
  for (int i=0; i<NCELL; i++)
  {
    for (int j=0; j<nSidePMTs; j++)
    {
      sideCharge[i][j] = sideCharge[i][j]/primaryCharge[i];
    }
  }
  
  for (int i=1; i<NCELL; i++)
  {
    primaryCharge[i] = primaryCharge[0]/primaryCharge[i];
  }
  primaryCharge[0] = 1.0;


  //Open the default output stream
  if (STREAM_FLAG==0)
  {
    if (freopen("/dev/null", "w", stdout)==NULL) exit(1);
  }
  

  //Output the basis data if desired
  if (TEST_FLAG_RECON==1)
  {
    printf("\nTo output the charge data enter 1, otherwise enter 0\n");
    printf(" > ");
    scanf("%d", &continueFlag);    
    
    if (continueFlag)
    {
      printf("Primary PMT Charge:\n");    
      for (int i=0; i<NCELL; i++)
      {
        printf("  %d %lf\n", i+1, primaryCharge[i]);
      }
      
      printf("\nSide PMT Charge:\n");    
      for (int i=0; i<NCELL; i++)
      {
        for (int j=0; j<nSidePMTs; j++)
        {
          printf("  %d %d %lf\n", i+1, j+1, sideCharge[i][j]);
        }
      }
      printf("\n");
    }
  } //End optional output if testing


  //Initialize other needed variables
  eventCnt=0;
  firstEvent=1;
  initString(string0, stringLen_g);
  
  
  
  //&&&&&&&&&&&&&&&&&&&&&&&&
  //& Read and Reconstruct &
  //&&&&&&&&&&&&&&&&&&&&&&&&    

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

    //Read the PMT waveforms
    readHeaderBlock(event);


    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //& Initialize the Hit Cells &
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     
    //Add energy to the cells:
    
    //  Loop over all the cells and see if there are AT LEAST three orthogonal PMTs with light in
    //  them. If this is the case then check the total number of hit PMTs is satisfactory. For the
    //  half instrumented case three is the maximum number, for the full instrumented case the
    //  required number are controled with the N_HIT_PMTS_REQUIRED. If the hit condition is met then
    //  add photoelectrons to the cell. The number of photoelectrons is computed by assuming that
    //  all of the charge in the hit PMTs comes only from the cell. The minimum of this set is taken
    //  as the initial value. 
      
    if (firstEvent==1)
    {
      initDet = new detector_t;
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

          //Check which PMTs have light in them
          index=0;
          int nHitPMTs=0;
          x0 = initDet->component[0][j][k];
          y0 = initDet->component[i][0][k];
          z0 = initDet->component[i][j][0];
          
          if (MIRROR_FLAG==0) //Use ALL of the PMTs
          {
            xN = initDet->component[NCELL+1][j][k];          
            yN = initDet->component[i][NCELL+1][k];          
            zN = initDet->component[i][j][NCELL+1];  
            
            if (x0>=HIT_PMT_THRESHOLD) nHitPMTs++;
            if (y0>=HIT_PMT_THRESHOLD) nHitPMTs++;
            if (z0>=HIT_PMT_THRESHOLD) nHitPMTs++;
            if (xN>=HIT_PMT_THRESHOLD) nHitPMTs++;
            if (yN>=HIT_PMT_THRESHOLD) nHitPMTs++;
            if (zN>=HIT_PMT_THRESHOLD) nHitPMTs++;
            
            if (nHitPMTs>=N_HIT_PMTS_REQUIRED) { index=1; }
          }
          else                //Use only the "0-side" PMTs
          {
            if ( x0>=HIT_PMT_THRESHOLD && y0>=HIT_PMT_THRESHOLD && z0>=HIT_PMT_THRESHOLD )
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
          
          x0 = initDet->component[0][j][k]*primaryCharge[i-1];
          y0 = initDet->component[i][0][k]*primaryCharge[j-1];
          z0 = initDet->component[i][j][0]*primaryCharge[k-1];          
          
          if (MIRROR_FLAG==0)
          {
            xN = initDet->component[NCELL+1][j][k]*primaryCharge[NCELL-i];
            yN = initDet->component[i][NCELL+1][k]*primaryCharge[NCELL-j];
            zN = initDet->component[i][j][NCELL+1]*primaryCharge[NCELL-k];
            
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

            cellCnt++;
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
    
    //Output the lists if testing   
    if (TEST_FLAG_RECON==1)
    {
      printf("\nTo print the PMT list enter 1, otherwise enter 0\n");
      printf(" > ");
      scanf("%d", &continueFlag);

      if (continueFlag)
      {
        printf("PMTs:\n");
        for (int i=0; i<pmtCnt; i++)
        {
          printf("%d %d %d %lf\n", pmtList[i].nx, pmtList[i].ny, pmtList[i].nz, pmtList[i].e);
        }
        printf("\n");
      }
      
      printf("\nTo print the cell list enter 1, otherwise enter 0\n");
      printf(" > ");
      scanf("%d", &continueFlag);      

      if (continueFlag)
      {
        printf("Cells:\n");
        for (int i=0; i<cellCnt; i++)
        {
          printf("%d %d %d %lf\n", cellList[i].nx, cellList[i].ny, cellList[i].nz, cellList[i].e);
        }
        printf("\n");        
      }
    } //End optional test statements
    
 
    //&&&&&&&&&&&&&&&
    //& Reconstruct &
    //&&&&&&&&&&&&&&&
    
    if (getTotalCharge(pmtList, pmtCnt)==-1.0*pmtCnt)  //If you removed all PMTs then there is nothing to do!
    {
      printf("No PMTs left to reconstruct for event %d!\n", eventNum);
      fprintf(out, "Event %d\n0\n", eventNum);
      continue;
    }

    chiSqrInit = getChiSqrTeflon(primaryCharge, sideCharge, pmtList, cellList, pmtCnt, cellCnt);

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

      //Call the gradient search and get the new chi-squared
      gradientSearchTeflon(primaryCharge, sideCharge, pmtList, cellList, pmtCnt, cellCnt);
      chiSqrNew = getChiSqrTeflon(primaryCharge, sideCharge, pmtList, cellList, pmtCnt, cellCnt);
              
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
    
    //Check the minizmation
    chiSqrFinal = getChiSqrTeflon(primaryCharge, sideCharge, pmtList, cellList, pmtCnt, cellCnt);
    
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
       cellList[i].t, getEnergy(cellList[i].e, conversionSlope, conversionIntercept) );
    }


    //&&&&&&&&&&&&
    //& Clean-up &
    //&&&&&&&&&&&&

    delete [] pmtList;
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
  delete [] primaryCharge;
  for (int i=0; i<NCELL; i++) { delete [] sideCharge[i]; }
  delete [] sideCharge;
  delete initDet;
  
  //Close files    
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
