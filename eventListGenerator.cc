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
#include <gsl/gsl_randist.h>

//&&&&&&&&&&&&&&&&&&
//& Custom Headers &
//&&&&&&&&&&&&&&&&&&
#include "analysis.hh"
#include "detectorParameters.hh"
#include "reconLib.hh"
#include "parsing.hh"
#include "randomNumberGenerator.hh"

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
  char string0[stringLen_g];              //TO DO: Go through and comment these
  char string1[stringLen_g];              //
  char tempString[stringLen_g];           //
  int genFlag, seed, seedFlag;            //
  int nCells, latticeFlag;                //
  int mirrorFlag, guideFlag;              //
  int configCheck;                        //
  double reactorPower, detectorBaseline;  //
  double backgroundRatekHz;               //  

  //Variables for the fake data
  char dataFileBaseName[stringLen_g];   //
  char ** dataFileName;                 //
  char ** ibdEventFileName;             //
  char ** backgroundEventFileName;      //
  char tempFileName[stringLen_g];       //
  int nEntries, nSignalFiles;           //
  int nBackgroundFiles, nSignalEvents;  //
  int nBackgroundEvents;                //
  int nNeededSignalFiles, eventsLeft;   //
  int * nEventsFromFile, eventIndex;    //
  int totalEventCount;                  //
  int numberOfFiles;                    //
  double detectorVolume, detectorMass;  //
  double neutrinoRateHz;                //
  double avgNumSignalEvents;            //
  long double runTime_s, runTime_ns;    //
  gsl_rng *  rGen;                      //
  const gsl_rng_type * T;               //
  
  //Event lists
  nuLatCorrelatedEvent_t  * correlatedEventList;      //
  nuLatUncorrelatedEvent_t  * uncorrelatedEventList;  //
  nuLatEvent_t  * eventList;                          //
  
  //Files
  FILE * config, * list;        //
  FILE  * event, * data, * out; //


  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&

  //Read the program control variables
  config = fopen("inputs/eventListGeneratorConfig.txt", "r");
  if (config==NULL)
  {
    printf("Could not open the file inputs/eventListGeneratorConfig.txt! Exiting...\n");
    exit(1);
  }

  if (getNumEntries(config)!=9)
  {
    printf("The file eventListGeneratorConfig.txt is improperly formatted! Exiting...\n");
    exit(1);
  }

  for (int i=0; i<9; i++)
  {
    initString(string0, stringLen_g);
    initString(string1, stringLen_g);
    getEntry(config, string0, stringLen_g);  
    
    if      (i==0)  { sscanf(string0, "%s %d", string1, &genFlag);            }
    else if (i==1)  { sscanf(string0, "%s %d", string1, &seedFlag);           }
    else if (i==2)  { sscanf(string0, "%s %d", string1, &nCells);             }
    else if (i==3)  { sscanf(string0, "%s %d", string1, &latticeFlag);        }
    else if (i==4)  { sscanf(string0, "%s %d", string1, &mirrorFlag);         }
    else if (i==5)  { sscanf(string0, "%s %d", string1, &guideFlag);          }
    else if (i==6)  { sscanf(string0, "%s %lf", string1, &reactorPower);      }
    else if (i==7)  { sscanf(string0, "%s %lf", string1, &detectorBaseline);  }
    else            { sscanf(string0, "%s %lf", string1, &backgroundRatekHz); }
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

  //Confirm the configuration 
  if (NCELL!=nCells || LATTICE_FLAG!=latticeFlag || MIRROR_FLAG!=mirrorFlag ||
    GUIDE_FLAG!=guideFlag)
  {
    printf("There is a mismatch between the macros defined in detectorParameters.hh and in eventListGeneratorConfig.txt!\n");
    printf("Check what macros were for the event files that your are using\n");
    printf("Adjust this before running. Exiting...\n");
    exit(1);
  }
  
  printf("Dear user this is the configuration that you inputted:\n");
  printf("  The seed is %d\n", int(seed));

  printf("  The detector size will be %dx%dx%d cells.\n", NCELL, NCELL, NCELL);

  if (mirrorFlag) { printf("  Mirrors are ON.\n"); }
  else { printf("  Mirrors are OFF.\n"); }

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

  printf("  The reactor power is %4.1lf MW\n", reactorPower);
  printf("  The detector baseline is %3.1lf m\n", detectorBaseline);

  printf("\nIf this configuration is correct enter an a non-zero integer otherwise enter 0\n");
  printf("> ");
  scanf("%d", &configCheck);
  
  if (!configCheck)
  {
    printf("Edit the configuration in inputs/eventListGeneratorConfig.txt. Exiting...\n");
    return 0;
  }
  
  //Get the names of the correlated and uncorrelated events files
  list = fopen("inputs/eventListGenerator.in", "r");
  if (list==NULL)
  {
    printf("Cannot open the file inputs/eventListGenerator.in! Exiting...\n");
    exit(1); 
  }
  nEntries = getNumEntries(list);
  
  nSignalFiles=0;
  nBackgroundFiles=0; 
  for (int i=0; i<nEntries; i++)
  {
    getEntry(list, tempFileName, stringLen_g);  
    if (strstr(tempFileName, "IBD")) nSignalFiles++;
    else nBackgroundFiles++;
  }
  
  ibdEventFileName = new char* [nSignalFiles];
  backgroundEventFileName = new char* [nBackgroundFiles];
  
  for (int i=0; i<nSignalFiles; i++)
  {
    ibdEventFileName[i] = new char [stringLen_g];
    initString(ibdEventFileName[i], stringLen_g);
  }
  
  for (int i=0; i<nBackgroundFiles; i++)
  {
    backgroundEventFileName[i] = new char [stringLen_g];
    initString(backgroundEventFileName[i], stringLen_g);
  }
  
  nSignalFiles=0;
  nBackgroundFiles=0; 
  rewind(list);
  for (int i=0; i<nEntries; i++)
  {
    getEntry(list, tempFileName, stringLen_g);  
    if (strstr(tempFileName, "IBD"))
    {
      copyString(ibdEventFileName[nSignalFiles], tempFileName, strlen(tempFileName));
      nSignalFiles++;
    }
    else
    {
      copyString(backgroundEventFileName[nBackgroundFiles], tempFileName, strlen(tempFileName));
      nBackgroundFiles++;
    }
  }
  fclose(list);
  
  //Get the runtime
  
  //TO DO: Actualy determine the number of events in each file instead of assuming that there are 1000
  
  nBackgroundEvents = nBackgroundFiles*1000;  //****Assuming that each file has 1000 events,
                                              //****should be good enough...
  runTime_s = nBackgroundEvents/(backgroundRatekHz*1000.0);
  runTime_ns = runTime_s*1E9;
  
  //Get the number of neutrino events  
  detectorVolume = pow(NCELL*CELL_DIM, 3.0);
  detectorMass = 1.021*detectorVolume*1E-6;
  neutrinoRateHz = (NEUTRINO_EVENT_RATE_STANDARD/(365.25*24.0*60.0*60.0))*detectorMass
    *(1000.0/detectorBaseline)*(1000.0/detectorBaseline)*(reactorPower/1000.0);
  
  avgNumSignalEvents = neutrinoRateHz*runTime_s;
  if (avgNumSignalEvents>30.0)
  {
    nSignalEvents = int(round(avgNumSignalEvents + gsl_ran_gaussian(rGen, sqrt(avgNumSignalEvents))));
  }
  else
  {
    nSignalEvents = gsl_ran_poisson(rGen, avgNumSignalEvents);
  }
  
  if (nSignalEvents>nSignalFiles*1000)
  {
    printf("\n");
    printf("  Caution! There are not enough neutrino events listed in eventListGenerator.in. Continuing...\n");
  }

  //Write the run info
  printf("  The runtime is %Lf s\n", runTime_s);
  printf("  with %d neutrino events\n", nSignalEvents);
  printf("  and %d background events\n", nBackgroundEvents);
  
  //Get the number of files to generate and their names
  time_t mytime = time(NULL);
  initString(tempString, stringLen_g);
  sprintf(tempString, "fakeData/%s", ctime(&mytime));
  
  initString(dataFileBaseName, stringLen_g);
  copyString(dataFileBaseName, tempString, strlen(tempString)-1);
  for (int i=0; i<int(strlen(dataFileBaseName)); i++)
  {
    if (dataFileBaseName[i]==' ') dataFileBaseName[i] = '_';
  }  
  
  totalEventCount = 2*nSignalEvents + nBackgroundEvents;
  numberOfFiles = int(round(totalEventCount/1000.0));
  
  dataFileName = new char* [numberOfFiles];
  for (int i=0; i<numberOfFiles; i++)
  {
    dataFileName[i] = new char [stringLen_g];
    initString(dataFileName[i], stringLen_g);
    copyString(dataFileName[i], dataFileBaseName, strlen(dataFileBaseName));
    
    initString(tempString, stringLen_g);
    sprintf(tempString, "_%d.dat", i);
    strcat(dataFileName[i], tempString);
    
    data = fopen(dataFileName[i], "w");
    if (data==NULL)
    {
      printf("Could not open the file %s! Exiting...\n", dataFileName[i]);
      exit(1);
    }

    printDetectorParameters(data);
    fprintf(data, "The reactor power is %4.1lf MW\n", reactorPower);
    fprintf(data, "The detector baseline is %3.1lf m\n", detectorBaseline);
    fprintf(data, "There are %d signal events\n", nSignalEvents);
    fprintf(data, "There are %d background events\n", nBackgroundEvents);
    fprintf(data, "The inputted signal-to-background is %lf.\n",
      double(nSignalEvents)/double(nBackgroundEvents) );
    fprintf(data, "The runtime is %Lf s.\n", runTime_ns*1E-9);      
    fclose(data);
  }

  //Generate input file for the trigger program    
  out = fopen("inputs/getTiggerEvents.in", "w");
  if (out==NULL)
  {
    printf("Could not open the file inputs/getTiggerEvents.dat! Exiting...\n");
    exit(1);
  }

  for (int i=0; i<numberOfFiles; i++)
  {
    fprintf(out, "%s\n", dataFileName[i]);  
  }
  fclose(out);
  

  //&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Generate the Data Key &
  //&&&&&&&&&&&&&&&&&&&&&&&&&
  
  //Get the correlated events
  correlatedEventList = new nuLatCorrelatedEvent_t [nSignalEvents];
  initializeEventList(correlatedEventList, nSignalEvents);
  for (int i=0; i<nSignalEvents; i++)
  {
    correlatedEventList[i].eventID = i;
    copyString(correlatedEventList[i].eventFile, ibdEventFileName[i/1000], stringLen_g);
    correlatedEventList[i].eventNumber = i%1000;
    correlatedEventList[i].promptEventTime = runTime_ns*rng(genFlag, rGen);
  }
  
  nNeededSignalFiles = nSignalEvents/1000 + 1;
  nEventsFromFile = new int [nNeededSignalFiles];
 
  eventsLeft = nSignalEvents;
  for (int i=0; i<nNeededSignalFiles; i++)
  {
    if (eventsLeft<=1000)
    {
      nEventsFromFile[i] = eventsLeft;
    }
    else
    {
      nEventsFromFile[i] = 1000;
      eventsLeft -= 1000;
    }
  } 
  
  eventIndex=0;    
  for (int i=0; i<nNeededSignalFiles; i++)
  {
    event = fopen(ibdEventFileName[i], "r");
    if (event==NULL)
    {
      printf("Could not open the file %s! Skipping to the next file...\n", ibdEventFileName[i]);
      continue;
    }
    
    for (int j=0; j<nEventsFromFile[i]; j++)
    {
      long double captureTime = getNeutronCaptureTime(event);
      
      if (captureTime>0.0)
      {
        correlatedEventList[eventIndex].delayedEventTime =  captureTime + 
          correlatedEventList[eventIndex].promptEventTime;
      }
      else 
      {
        printf("Problem with %s event %d\n", ibdEventFileName[i], j);
      }
          
      eventIndex++;
    }
    
    fclose(event);
  }
  
  delete [] nEventsFromFile;

  //Get the uncorrelated events
  uncorrelatedEventList = new nuLatUncorrelatedEvent_t [nBackgroundEvents];
  initializeEventList(uncorrelatedEventList, nBackgroundEvents);
  for (int i=0; i<nBackgroundEvents; i++)
  {
    uncorrelatedEventList[i].eventID = i;
    copyString(uncorrelatedEventList[i].eventFile, backgroundEventFileName[i/1000], stringLen_g);
    uncorrelatedEventList[i].eventNumber = i%1000;
    uncorrelatedEventList[i].eventTime = runTime_ns*rng(genFlag, rGen);
  }

  //Combine the all events
  eventList = new nuLatEvent_t [totalEventCount];
  
  for (int i=0; i<nSignalEvents; i++)
  {
    eventList[2*i].eventID = 2*i;
    copyString(eventList[2*i].eventFile, correlatedEventList[i].eventFile, stringLen_g);
    eventList[2*i].correlatedFlag = 1;
    eventList[2*i].eventNumber = correlatedEventList[i].eventNumber;
    eventList[2*i].eventTime = correlatedEventList[i].promptEventTime;
    
    eventList[2*i+1].eventID = 2*i+1;
    copyString(eventList[2*i+1].eventFile, correlatedEventList[i].eventFile, stringLen_g);
    eventList[2*i+1].correlatedFlag = 2;
    eventList[2*i+1].eventNumber = correlatedEventList[i].eventNumber;
    eventList[2*i+1].eventTime = correlatedEventList[i].delayedEventTime;
  }
  
  for (int i=2*nSignalEvents; i<totalEventCount; i++)
  {
    eventList[i].eventID = i;
    copyString(eventList[i].eventFile, uncorrelatedEventList[i-2*nSignalEvents].eventFile, stringLen_g);
    eventList[i].correlatedFlag = 0;
    eventList[i].eventNumber = uncorrelatedEventList[i-2*nSignalEvents].eventNumber;
    eventList[i].eventTime = uncorrelatedEventList[i-2*nSignalEvents].eventTime;  
  }
  
  
  //&&&&&&&&&&
  //& Output &
  //&&&&&&&&&&
  
  timeOrderEvents(eventList, totalEventCount);
  int eventNum=0;
  for (int i=0; i<numberOfFiles; i++)
  {
    data = fopen(dataFileName[i], "a");
    if (data==NULL)
    {
      printf("Could not open the file %s! Skipping to the next file...\n", dataFileName[i]);
      continue;
    }
    
    if (eventNum+1000<=totalEventCount)
    {
      for (int j=0; j<1000; j++)
      {
        fprintf(data, "%d\n", i*1000+j);
        fprintf(data, "  %s\n", eventList[i*1000+j].eventFile);
        fprintf(data, "  %d\n", eventList[i*1000+j].eventNumber);
        fprintf(data, "  %d\n", eventList[i*1000+j].correlatedFlag);
        fprintf(data, "  %Lf\n", eventList[i*1000+j].eventTime);
        
        eventNum++;     
      }
    }
    else
    {
      for (int j=0; j<totalEventCount-eventNum; j++)
      {
        fprintf(data, "%d\n", i*1000+j);
        fprintf(data, "  %s\n", eventList[i*1000+j].eventFile);
        fprintf(data, "  %d\n", eventList[i*1000+j].eventNumber);
        fprintf(data, "  %d\n", eventList[i*1000+j].correlatedFlag);
        fprintf(data, "  %Lf\n", eventList[i*1000+j].eventTime);
        
        eventNum++;  
      }    
    }
    
    fclose(data);
  }

  
  //&&&&&&&
  //& End &
  //&&&&&&&
  
  for (int i=0; i<nSignalFiles; i++)
  {
    delete [] ibdEventFileName[i];
  }
  
  for (int i=0; i<nBackgroundFiles; i++)
  {
    delete [] backgroundEventFileName[i];
  }
  
  for (int i=0; i<numberOfFiles; i++)
  {
    delete [] dataFileName[i];  
  }
  
  delete [] ibdEventFileName;
  delete [] backgroundEventFileName;
  delete [] dataFileName;
  delete [] correlatedEventList;
  delete [] uncorrelatedEventList;
  delete [] eventList;

  return 0;
}

//*************************************************************************************************
