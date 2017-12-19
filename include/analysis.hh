#ifndef ANALYSIS_H
#define ANALYSIS_H

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
#include "detectorParameters.hh"
#include "reconLib.hh"
#include "parsing.hh"

//&&&&&&&&&&&&&
//& Constants &
//&&&&&&&&&&&&&

// Event rate in plastic scintillator:
//  - 1 ton
//  - 1 km
//  - 1 GW
//  -> 180 neutrinos per year 

#define NEUTRINO_EVENT_RATE_STANDARD 180.0

//&&&&&&&&&&&&&&&&&&&
//& Data Structures &
//&&&&&&&&&&&&&&&&&&&
typedef struct nuLatCorrelatedEvent_t
{
  int eventID;
  char eventFile[stringLen_g];
  int eventNumber;
  long double promptEventTime;
  long double delayedEventTime;
} nuLatCorrelatedEvent_t;

typedef struct nuLatUncorrelatedEvent_t
{
  int eventID;
  char eventFile[stringLen_g];
  int eventNumber;
  long double eventTime;
} nuLatUncorrelatedEvent_t;

typedef struct nuLatEvent_t
{
  int eventID;
  char eventFile[stringLen_g];
  int eventNumber;
  int correlatedFlag;
  long double eventTime;
} nuLatEvent_t;

typedef struct nuLatTriggerEvent_t
{
  int triggerID;
  int numberOfEvents;
  nuLatEvent_t * event;
} nuLatTriggerEvent_t;


//&&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Declarations &
//&&&&&&&&&&&&&&&&&&&&&&&&&

void initializeEventList(nuLatCorrelatedEvent_t * list, int length);
void initializeEventList(nuLatUncorrelatedEvent_t * list, int length);
void initializeEventList(nuLatEvent_t * list, int length);
long double getNeutronCaptureTime(FILE * f);
void readBlock(FILE * f, int headerType);
void timeOrderEvents(nuLatEvent_t * list, int length);
nuLatEvent_t readNuLatEvent(FILE * f);
void goToEventNumber(FILE * f, int eventNum, int correlatedFlag);
void readEventFileEvent(FILE * f, int correlatedFlag);
int getNumberOfEventsInList(FILE * f);

//*************************************************************************************************

#endif
