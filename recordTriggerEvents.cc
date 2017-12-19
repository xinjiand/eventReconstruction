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
#define COINCIDENCE_WINDOW 75000.0 //The coincidence window in ns

//&&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Declarations &
//&&&&&&&&&&&&&&&&&&&&&&&&&


//*************************************************************************************************

//&&&&&&&&
//& Main &
//&&&&&&&&
int main ()
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Declare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  //Variables for recording the trigger events
  char string0[stringLen_g];
  char ** triggerFileName;
  char ** eventListFileName;
  char ** fakeDataFileName;
  int nEntries, numberOfLines;
  int numberOfTriggers;
  int numberOfEventsInList;
  nuLatTriggerEvent_t * coincidenceEvent;
  nuLatEvent_t * triggerEvent;
  nuLatEvent_t currentEvent;
  long double windowStartTime, windowStopTime;
  long double lastTimeOfEventInList;
  
  //Files
  FILE * list, * trig, * eventList;
  FILE * fakeData;
  

  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&

  list = fopen("inputs/recordTriggerEvents.in", "r");
  if (list==NULL)
  {
    printf("Could not open the file inputs/recordTriggerEvents.in! Exiting...\n");
    exit(1);
  }

  nEntries = getNumEntries(list);
  triggerFileName = new char* [nEntries];
  eventListFileName = new char* [nEntries];
  fakeDataFileName = new char* [nEntries];
  for (int i=0; i<nEntries; i++)
  {
    triggerFileName[i] = new char [stringLen_g];
    eventListFileName[i] = new char [stringLen_g];
    fakeDataFileName[i] = new char [stringLen_g];
    
    initString(triggerFileName[i], stringLen_g);
    initString(eventListFileName[i], stringLen_g);
    initString(fakeDataFileName[i], stringLen_g);
    
    getEntry(list, triggerFileName[i], stringLen_g);
    
    copyString(eventListFileName[i], triggerFileName[i], strlen(triggerFileName[i])-12);
    strcat(eventListFileName[i], ".dat");
    
    copyString(fakeDataFileName[i], triggerFileName[i], strlen(triggerFileName[i])-12);
    strcat(fakeDataFileName[i], "_Fake_Data.dat");
  }
  
  
  //&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Get Coincident Events &
  //&&&&&&&&&&&&&&&&&&&&&&&&&

  lastTimeOfEventInList = 0.0;
  for (int i=0; i<nEntries; i++)
  {
    //Open the needed files
    trig = fopen(triggerFileName[i], "r");
    if (trig==NULL)
    {
      printf("Could not open the file %s! Skipping to the next file...\n", triggerFileName[i]);
      continue;
    }
    
    eventList = fopen(eventListFileName[i], "r");
    if (eventList==NULL)
    {
      printf("Could not open the file %s! Skipping to the next file...\n", eventListFileName[i]);
      fclose(trig);
      continue;    
    }
    
    fakeData = fopen(fakeDataFileName[i], "w");
    if (fakeData==NULL)
    {
      printf("Could not open the file %s! Skipping to the next file...\n", fakeDataFileName[i]);
      fclose(eventList);
      fclose(trig);
      continue;
    }
    
    //Get the trigger events
    numberOfLines = 0;
    while (fgets(string0, stringLen_g, trig)!=NULL) { numberOfLines++; }
    rewind(trig);
    
    if (numberOfLines%5!=0) 
    {
      printf("There was a formatting error in the file %s! Skipping to the next trigger file...\n", 
        triggerFileName[i]);
      fclose(trig);
      fclose(eventList);
      fclose(fakeData);
      continue;
    }
    
    numberOfTriggers = numberOfLines/5;
    triggerEvent = new nuLatEvent_t [numberOfTriggers];
    initializeEventList(triggerEvent, numberOfTriggers);
    
    coincidenceEvent = new nuLatTriggerEvent_t [numberOfTriggers];
    
    for (int j=0; j<numberOfTriggers; j++)
    {
      for (int k=0; k<5; k++)
      {
        initString(string0, stringLen_g);
        fgets(string0, stringLen_g, trig);
        
        if      (k==0) { sscanf(string0, "%d", &triggerEvent[j].eventID);        }
        else if (k==1) { sscanf(string0, "%s", triggerEvent[j].eventFile);       }       
        else if (k==2) { sscanf(string0, "%d", &triggerEvent[j].eventNumber);    }
        else if (k==3) { sscanf(string0, "%d", &triggerEvent[j].correlatedFlag); }
        else           { sscanf(string0, "%Lf", &triggerEvent[j].eventTime);     }
      }
      
      coincidenceEvent[j].triggerID = j;
      coincidenceEvent[j].numberOfEvents = 0;
    }
    
    //For each trigger determine the number of events that fall within the window
    for (int j=0; j<numberOfTriggers; j++)
    {
      windowStartTime = triggerEvent[j].eventTime - COINCIDENCE_WINDOW;
      if (windowStartTime<0.0) windowStartTime = 0.0;
      windowStopTime = triggerEvent[j].eventTime;
      
      if (windowStartTime<lastTimeOfEventInList)
      {
        printf("There are missing events due to file splitting for trigger %d\n", j);
      }
      
      //Set-up the event list for reading
      rewind(eventList);
      numberOfEventsInList = getNumberOfEventsInList(eventList);
      if (numberOfEventsInList==-1)
      {
        printf("Error in formatting of %s! Skipping to the next trigger file...\n", eventListFileName[i]);
        break;
      }
      
      rewind(eventList);
      for (int k=0; k<36; k++)
      {
        initString(string0, stringLen_g);
        fgets(string0, stringLen_g, eventList);
      }
      
      //Read the events in the list to get the numer of events in coincidence with the trigger
      for (int k=0; k<numberOfEventsInList; k++)
      {
        currentEvent = readNuLatEvent(eventList);
        
        if (windowStartTime<=currentEvent.eventTime && currentEvent.eventTime<windowStopTime)
        {
          coincidenceEvent[j].numberOfEvents++;
        }
      }
      
      //Get the coincidence events
      rewind(eventList);
      for (int k=0; k<36; k++)
      {
        initString(string0, stringLen_g);
        fgets(string0, stringLen_g, eventList);
      }
      
      coincidenceEvent[j].event = new nuLatEvent_t [coincidenceEvent[j].numberOfEvents+1];
      int eventIndex=0;
      for (int k=0; k<numberOfEventsInList; k++)
      {
        currentEvent = readNuLatEvent(eventList);
        
        if (windowStartTime<=currentEvent.eventTime && currentEvent.eventTime<windowStopTime)
        {
          coincidenceEvent[j].event[eventIndex].eventID        = currentEvent.eventID;
          copyString(coincidenceEvent[j].event[eventIndex].eventFile, currentEvent.eventFile,
            strlen(currentEvent.eventFile));
          coincidenceEvent[j].event[eventIndex].eventNumber    = currentEvent.eventNumber;
          coincidenceEvent[j].event[eventIndex].correlatedFlag = currentEvent.correlatedFlag;
          coincidenceEvent[j].event[eventIndex].eventTime      = currentEvent.eventTime;
          
          eventIndex++;
        }
      }
      
      coincidenceEvent[j].event[eventIndex].eventID        = triggerEvent[j].eventID;
      copyString(coincidenceEvent[j].event[eventIndex].eventFile, triggerEvent[j].eventFile,
        strlen(triggerEvent[j].eventFile));
      coincidenceEvent[j].event[eventIndex].eventNumber    = triggerEvent[j].eventNumber;
      coincidenceEvent[j].event[eventIndex].correlatedFlag = triggerEvent[j].correlatedFlag;
      coincidenceEvent[j].event[eventIndex].eventTime      = triggerEvent[j].eventTime;   
    }
    
    //Output the coincidence events
    for (int j=0; j<numberOfTriggers; j++)
    {
      fprintf(fakeData, "%d\n", coincidenceEvent[j].triggerID);
      fprintf(fakeData, "  %d\n", coincidenceEvent[j].numberOfEvents);
      
      for (int k=0; k<coincidenceEvent[j].numberOfEvents+1; k++)
      {
        fprintf(fakeData, "    %d\n", coincidenceEvent[j].event[k].eventID);
        fprintf(fakeData, "      %s\n", coincidenceEvent[j].event[k].eventFile);
        fprintf(fakeData, "      %d\n", coincidenceEvent[j].event[k].eventNumber);
        fprintf(fakeData, "      %Lf\n", coincidenceEvent[j].event[k].eventTime);
      }
    }
    
    //Get the time of the last event in the file
    rewind(eventList);
    for (int k=0; k<36; k++)
    {
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, eventList);
    }
    
    for (int k=0; k<numberOfEventsInList; k++)
    {
      currentEvent = readNuLatEvent(eventList);
    }
    
    lastTimeOfEventInList = currentEvent.eventTime;
    
    //Clean-up
    fclose(trig);
    fclose(eventList);
    fclose(fakeData);    
  }
  
  
  //&&&&&&&
  //& End &
  //&&&&&&&
  
  fclose(list);
  
  return 0;
}

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&

//*************************************************************************************************
