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
#include "analysis.hh"
#include "detectorParameters.hh"
#include "parsing.hh"
#include "reconLib.hh"

//&&&&&&&&&&&&&
//& Constants &
//&&&&&&&&&&&&&
#define ENERGY_THRESHOLD_FOR_TRIGGER 0.20   //Based off of 400 keVee neutron tag with the worst
                                            //collection possible in a 15x15x15 NuLat with half
                                            //instrumentation
#define HIT_PMT_THRESHOLD 10.0               //Set roughly based off of the mis-channled light
#define GET_AND_PRINT_EVENT_NUMBER 1        //Get and print the event number of a trigger

//&&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Declarations &
//&&&&&&&&&&&&&&&&&&&&&&&&&

int thresholdTrigger(component_t * pmtList, int length);
int cellMultiplicityTrigger(component_t * pmtList, int length);

//*************************************************************************************************

//&&&&&&&&
//& Main &
//&&&&&&&&
int main ()
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Declare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  //Variables for determining trigger events
  char dataFileName[stringLen_g];     //
  char triggerFileName[stringLen_g];  //
  char string0[stringLen_g];          //
  int nEntries, correlatedFlag;       //
  int nLines, triggerEventNumber;     //
  int nPositronLines, nNeutronLines;  //
  int trigger;                        //
  nuLatEvent_t currentEvent;          //
  component_t * pmtList;              //

  //Files
  FILE * list, * in;    //   
  FILE * event, * trig; //
  FILE * out;           // 


  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&

  list = fopen("inputs/getTiggerEvents.in", "r");
  if (list==NULL)
  {
    printf("Could not open the file inputs/getTiggerEvents.dat! Exiting...\n");
    exit(1);
  }
  
  out = fopen("inputs/recordTriggerEvents.in", "w");
  if (out==NULL)
  {
    printf("Could not open the file inputs/recordTriggerEvents.in! Exiting...\n");
    exit(1);
  }

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Determine Trigger Events &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  nEntries = getNumEntries(list);
  for (int i=0; i<nEntries; i++)
  {
    //Open the data listing file
    initString(dataFileName, stringLen_g);
    getEntry(list, dataFileName, stringLen_g);
    in = fopen(dataFileName, "r");
    if (in==NULL)
    {
      printf("Could not open the file %s! Continuing to the next file...\n", dataFileName);
      continue;
    }
    
    //Get the output file name and open
    initString(triggerFileName, stringLen_g);
    copyString(triggerFileName, dataFileName, strlen(dataFileName)-4);
    strcat(triggerFileName, "_trigger.dat");
    
    fprintf(out, "%s\n", triggerFileName);
    
    trig = fopen(triggerFileName, "w");
    if (trig==NULL)
    {
      printf("Could not open the file %s! Skipping to the next file...\n", triggerFileName);
      fclose(in);
      continue;
    }
    
    //Read the preamble
    for (int j=0; j<36; j++)
    {
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, in);
    }
    
    //Read the events and check if they are a valid trigger
    for (int eventNum=0; eventNum<1000; eventNum++)
    {
      if ((eventNum+1)%100==0) printf("%s %d events completed!\n", dataFileName, eventNum+1);
    
      currentEvent = readNuLatEvent(in);
      
      event = fopen(currentEvent.eventFile, "r");
      if (event==NULL)
      {
        printf("Could not open the file %s! Skipping to the next event...\n", currentEvent.eventFile);
        continue;
      }
      
      //Get the PMT charges
      goToEventNumber(event, currentEvent.eventNumber, currentEvent.correlatedFlag);
      
      if (currentEvent.correlatedFlag>0) correlatedFlag = 1;
      else correlatedFlag = 0;
      
      //Check the event header
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event);
      if (strstr(string0, "#Event")==NULL)
      {
        printf("Error at the event header!\n");
      }
      if (GET_AND_PRINT_EVENT_NUMBER)
      {
        initString(string0, stringLen_g);
        fgets(string0, stringLen_g, event);
      
        sscanf(string0, "%d", &nLines);
        for (int j=0; j<nLines; j++)
        {
          initString(string0, stringLen_g);
          fgets(string0, stringLen_g, event);
          if (j==0) sscanf(string0, "%d", &triggerEventNumber);
        }
      }
      else 
      {
        readBlock(event, 0);
      }
      
      //Check the Geant4 properties header
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event);
      if (strstr(string0, "#Geant4_properties")==NULL)
      {
        printf("Error at the Geant4 properties header!\n");
      }
      readBlock(event, 0);    

      //Check the energy deposits header
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event);
      if (strstr(string0, "#Energy_deposits")==NULL)
      {
        printf("Error at the energy deposits header!\n");
      }
      readBlock(event, correlatedFlag);
      
      //Check the light transport properties header
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event);
      if (strstr(string0, "#Light_transport_properties")==NULL)
      {
        printf("Error at the light transport header!\n");
      }
      readBlock(event, 0);  

      //Check the cell deposits header
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event);
      if (strstr(string0, "#Cells_deposits")==NULL)
      {
        printf("Error at the cell deposits header!\n");
      }
      readBlock(event, correlatedFlag);  
          
      //Check the PMT charges header
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event);
      if (strstr(string0, "#PMT_charges")==NULL)
      {
        printf("Error at the PMT charges header!\n");
      }

      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event);
      
      if (correlatedFlag==0)
      {
        sscanf(string0, "%d", &nLines);
        pmtList = new component_t [nLines];
        clearComponents(pmtList, nLines);
        
        for (int j=0; j<nLines; j++)
        {
          initString(string0, stringLen_g);
          fgets(string0, stringLen_g, event);
          
          sscanf(string0, "%d %d %d %lf", &pmtList[j].nx, &pmtList[j].ny, &pmtList[j].nz,
            &pmtList[j].e);
        }
      }
      else
      {
        sscanf(string0, "%d %d %d", &nLines, &nPositronLines, &nNeutronLines);
        
        if (currentEvent.correlatedFlag==1) { nLines = nPositronLines; }
        else { nLines = nNeutronLines; }

        pmtList = new component_t [nLines];
        clearComponents(pmtList, nLines);
        
        for (int j=0; j<nPositronLines; j++)
        {
          initString(string0, stringLen_g);
          fgets(string0, stringLen_g, event);
          
          if (currentEvent.correlatedFlag==1)
          {
            sscanf(string0, "%d %d %d %lf", &pmtList[j].nx, &pmtList[j].ny, &pmtList[j].nz,
              &pmtList[j].e);
          }        
        }

        for (int j=0; j<nNeutronLines; j++)
        {
          initString(string0, stringLen_g);
          fgets(string0, stringLen_g, event);
          
          if (currentEvent.correlatedFlag==2)
          {
            sscanf(string0, "%d %d %d %lf", &pmtList[j].nx, &pmtList[j].ny, &pmtList[j].nz,
              &pmtList[j].e);
          }        
        }  
      }
      
      fclose(event);

      //Determine the if the event is a valid trigger
      trigger = thresholdTrigger(pmtList, nLines);
      if (trigger==0) { continue; }
      
      trigger = cellMultiplicityTrigger(pmtList, nLines);
      if (trigger==0) { continue; }
            
      fprintf(trig, "%d\n", currentEvent.eventID);
      fprintf(trig, "  %s\n", currentEvent.eventFile);
      fprintf(trig, "  %d\n", currentEvent.eventNumber);
      fprintf(trig, "  %d\n", currentEvent.correlatedFlag);
      fprintf(trig, "  %Lf\n", currentEvent.eventTime);
      
      if (GET_AND_PRINT_EVENT_NUMBER)
      {
        printf("%d\n", triggerEventNumber);
      }
            
    } //End reading events
    
    fclose(in);
    fclose(trig);
  } //End loop over data list entries 

  
  //&&&&&&&
  //& End &
  //&&&&&&&
  
  fclose(list);
  fclose(out);
  
  return 0;
}

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&

int thresholdTrigger(component_t * pmtList, int length)
{
  //Declare needed variables
  int trigger=0;              //
  double peSum=0.0;           //
  double conversionSlope;     //
  double conversionIntercept; //
  double eventEnergy;         //
  
  //Initialize
  getEnergyScaleFactors(&conversionSlope, &conversionIntercept);
  
  //Get PE sum  
  for (int i=0; i<length; i++)
  {
    if (pmtList[i].nx==0) peSum += pmtList[i].e;
  }

  //Convert to energy and return if triggered
  eventEnergy = conversionSlope*peSum + conversionIntercept;
  
  if (eventEnergy>=ENERGY_THRESHOLD_FOR_TRIGGER) trigger=1;
  
  return trigger;
}

//*************************************************************************************************

int cellMultiplicityTrigger(component_t * pmtList, int length)
{
  //Declare needed variables
  component_t * xCellList;
  component_t * yCellList;
  component_t * hitCellList;
  int trigger=0;
  int numberX=0, numberY=0;
  int numberOfHitCells=0;

  //Get the cells viewed by the hit PMTs
  for (int i=0; i<length; i++)
  {
    if (pmtList[i].nx==0 && pmtList[i].e>HIT_PMT_THRESHOLD) numberX++;
    else if (pmtList[i].ny==0 && pmtList[i].e>HIT_PMT_THRESHOLD) numberY++;
  }
  
  numberX = numberX*NCELL;
  numberY = numberY*NCELL;
  xCellList = new component_t [numberX];
  yCellList = new component_t [numberY];
  
  clearComponents(xCellList, numberX);
  clearComponents(yCellList, numberY);
  
  numberX=0;
  numberY=0;
  for (int i=0; i<length; i++)
  {
    if (pmtList[i].nx==0 && pmtList[i].e>HIT_PMT_THRESHOLD)
    {
      for (int j=0; j<NCELL; j++)
      {
        xCellList[numberX].nx = j+1;
        xCellList[numberX].ny = pmtList[i].ny;
        xCellList[numberX].nz = pmtList[i].nz;
        numberX++;
      }
    }
    else if (pmtList[i].ny==0 && pmtList[i].e>HIT_PMT_THRESHOLD)
    {
      for (int j=0; j<NCELL; j++)
      {
        yCellList[numberY].nx = pmtList[i].nx;
        yCellList[numberY].ny = j+1;
        yCellList[numberY].nz = pmtList[i].nz;
        numberY++;
      }    
    }
  }
  
  //Get the hit cells and return
  for (int i=0; i<numberX; i++)
  {
    for (int j=0; j<numberY; j++)
    {
      if (xCellList[i].nx==yCellList[j].nx && xCellList[i].ny==yCellList[j].ny &&
        xCellList[i].nz==yCellList[j].nz)
      {
        numberOfHitCells++;
      }
    }
  }
  
  hitCellList = new component_t [numberOfHitCells];
  clearComponents(hitCellList, numberOfHitCells);
  numberOfHitCells=0;
  for (int i=0; i<numberX; i++)
  {
    for (int j=0; j<numberY; j++)
    {
      if (xCellList[i].nx==yCellList[j].nx && xCellList[i].ny==yCellList[j].ny &&
        xCellList[i].nz==yCellList[j].nz)
      {
        hitCellList[numberOfHitCells].nx = xCellList[i].nx;
        hitCellList[numberOfHitCells].ny = xCellList[i].ny;
        hitCellList[numberOfHitCells].nz = xCellList[i].nz;
        
        numberOfHitCells++;
      }
    }
  }  
  
  //Clean-up and return
  delete [] xCellList;
  delete [] yCellList;
  delete [] hitCellList;

  if (numberOfHitCells==1) trigger=1;
  return trigger;
}

//*************************************************************************************************
