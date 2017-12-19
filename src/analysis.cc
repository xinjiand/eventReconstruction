#include "analysis.hh"

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&

void initializeEventList(nuLatCorrelatedEvent_t * list, int length)
{
  for (int i=0; i<length; i++)
  {
    list[i].eventID = -1;
    initString(list[i].eventFile, stringLen_g);
    list[i].eventNumber = -1;
    list[i].promptEventTime = 0.0;
    list[i].delayedEventTime = 0.0;
  }
}

//*************************************************************************************************

void initializeEventList(nuLatUncorrelatedEvent_t * list, int length)
{
  for (int i=0; i<length; i++)
  {
    list[i].eventID = -1;
    initString(list[i].eventFile, stringLen_g);
    list[i].eventNumber = -1;
    list[i].eventTime = 0.0;
  }
}

//*************************************************************************************************

void initializeEventList(nuLatEvent_t * list, int length)
{
  for (int i=0; i<length; i++)
  {
    list[i].eventID = -1;
    initString(list[i].eventFile, stringLen_g);
    list[i].eventNumber = -1;
    list[i].correlatedFlag=-1;
    list[i].eventTime = 0.0;
  }
}

//*************************************************************************************************

long double getNeutronCaptureTime(FILE * f)
{
  //Declare needed variables
  char string0[stringLen_g];      //
  long double captureTime = -1.0; //
  double tempD;                   //
  int nLines;                     //
  int nPositronDep, nNeutronDep;  //
  deposit_t * neutronDeposits;    //
  
  //Read the event
  
  //Check the event header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#Event")==NULL)
  {
    printf("Error at the event header!\n");
  }
  readBlock(f, 0);
  
  //Check the Geant4 properties header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#Geant4_properties")==NULL)
  {
    printf("Error at the Geant4 properties header!\n");
  }
  readBlock(f, 0);    

  //Check the energy deposits header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#Energy_deposits")==NULL)
  {
    printf("Error at the energy deposits header!\n");
  }
  
  //Read the neutron deposits
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  sscanf(string0, " %d %d %d", &nLines, &nPositronDep, &nNeutronDep);
  
  neutronDeposits = new deposit_t [nNeutronDep];
  
  for (int i=0; i<nPositronDep; i++)
  {
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, f);    
  }
  
  for (int i=0; i<nNeutronDep; i++)
  {
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, f);
    
    sscanf(string0, "%d %lf %lf %lf %lf %lf %lf", &neutronDeposits[i].id, &neutronDeposits[i].x,
      &neutronDeposits[i].y, &neutronDeposits[i].z, &neutronDeposits[i].t, &tempD, 
        &neutronDeposits[i].energy);
  }    

  //Check the light transport properties header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#Light_transport_properties")==NULL)
  {
    printf("Error at the light transport header!\n");
  }
  readBlock(f, 0);  

  //Check the cell deposits header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#Cells_deposits")==NULL)
  {
    printf("Error at the cell deposits header!\n");
  }
  readBlock(f, 1);  
      
  //Check the PMT charges header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#PMT_charges")==NULL)
  {
    printf("Error at the PMT charges header!\n");
  }
  readBlock(f, 1);  
      
  //Check the PMT waveforms header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#PMT_waveforms")==NULL)
  {
    printf("Error at the PMT waveforms header!\n");
  }
  readBlock(f, 1);
      
  //Check the reconstruction header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#Reconstruction")==NULL)
  {
    printf("Error at the reconstruction header!\n");
  }
  readBlock(f, 1);

  //Get the neutron capture time
  sortDeposits(neutronDeposits, nNeutronDep, 1);
  int depositIndex=0;
  while (neutronDeposits[depositIndex].id!=2 && neutronDeposits[depositIndex].id!=8)
  {
    depositIndex++;
  }
  
  if (neutronDeposits[depositIndex].id==2 || neutronDeposits[depositIndex].id==8)
  {
    captureTime = neutronDeposits[depositIndex].t;
  }
  
  delete [] neutronDeposits;
  
  return captureTime;
}

//*************************************************************************************************

void readBlock(FILE * f, int headerType)
{
  //Declare needed variables
  char string0[stringLen_g];        //
  int nLines, nNeutron, nPositron;  //
  
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (headerType==0) sscanf(string0, "%d", &nLines);
  else sscanf(string0, "%d %d %d", &nLines, &nNeutron, &nPositron);
  
  for (int i=0; i<nLines; i++)
  {
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, f);  
  }
}

//*************************************************************************************************

void timeOrderEvents(nuLatEvent_t * list, int length)
{
  //Declare needed variables
  nuLatEvent_t temp;
  int swapCnt=1;
  
  //Bubble sort the array
  while (swapCnt)
  {
    swapCnt=0;
    for (int i=0; i<length-1; i++)
    {
      if (list[i].eventTime>list[i+1].eventTime)
      {
        temp.eventID = list[i].eventID;
        copyString(temp.eventFile, list[i].eventFile, stringLen_g);
        temp.eventNumber = list[i].eventNumber;
        temp.correlatedFlag = list[i].correlatedFlag;
        temp.eventTime = list[i].eventTime;
      
        list[i].eventID = list[i+1].eventID;
        copyString(list[i].eventFile, list[i+1].eventFile, stringLen_g);
        list[i].eventNumber = list[i+1].eventNumber;
        list[i].correlatedFlag = list[i+1].correlatedFlag;
        list[i].eventTime = list[i+1].eventTime;
        
        list[i+1].eventID = temp.eventID;
        copyString(list[i+1].eventFile, temp.eventFile, stringLen_g);
        list[i+1].eventNumber = temp.eventNumber;
        list[i+1].correlatedFlag = temp.correlatedFlag;
        list[i+1].eventTime = temp.eventTime;        
        
        swapCnt++;
      } //End ordering check      
    } //End loop over elments in the array
  } //End while loop
}

//*************************************************************************************************

nuLatEvent_t readNuLatEvent(FILE * f)
{
  char string0[stringLen_g];
  nuLatEvent_t event;
  
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  sscanf(string0, "%d", &event.eventID);
  
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  initString(event.eventFile, stringLen_g);
  int charIndex=0;
  for (int i=0; i<int(strlen(string0))-1; i++)
  {
    if (string0[i]!=' ') 
    {
      event.eventFile[charIndex] = string0[i];
      charIndex++;
    }
  }
  
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  sscanf(string0, "%d", &event.eventNumber);
  
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  sscanf(string0, "%d", &event.correlatedFlag);

  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  sscanf(string0, "%Lf", &event.eventTime);
  
  return event;
}

//*************************************************************************************************

void goToEventNumber(FILE * f, int eventNum, int correlatedFlag)
{
  if (correlatedFlag>0) correlatedFlag=1;
  for (int i=0; i<eventNum; i++)
  {
    readEventFileEvent(f, correlatedFlag);  
  }
}

//*************************************************************************************************

void readEventFileEvent(FILE * f, int correlatedFlag)
{
  //Read the event
  char string0[stringLen_g];
  
  //Check the event header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#Event")==NULL)
  {
    printf("Error at the event header!\n");
  }
  readBlock(f, 0);
  
  //Check the Geant4 properties header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#Geant4_properties")==NULL)
  {
    printf("Error at the Geant4 properties header!\n");
  }
  readBlock(f, 0);    

  //Check the energy deposits header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#Energy_deposits")==NULL)
  {
    printf("Error at the energy deposits header!\n");
  }
  readBlock(f, correlatedFlag);
  
  //Check the light transport properties header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#Light_transport_properties")==NULL)
  {
    printf("Error at the light transport header!\n");
  }
  readBlock(f, 0);  

  //Check the cell deposits header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#Cells_deposits")==NULL)
  {
    printf("Error at the cell deposits header!\n");
  }
  readBlock(f, correlatedFlag);  
      
  //Check the PMT charges header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#PMT_charges")==NULL)
  {
    printf("Error at the PMT charges header!\n");
  }
  readBlock(f, correlatedFlag);  
      
  //Check the PMT waveforms header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#PMT_waveforms")==NULL)
  {
    printf("Error at the PMT waveforms header!\n");
  }
  readBlock(f, correlatedFlag);
      
  //Check the reconstruction header
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, f);
  if (strstr(string0, "#Reconstruction")==NULL)
  {
    printf("Error at the reconstruction header!\n");
  }
  readBlock(f, correlatedFlag);
}

//*************************************************************************************************

int getNumberOfEventsInList(FILE * f)
{
  char s[stringLen_g];
  int numberOfLines=0;
  
  rewind(f);
  initString(s, stringLen_g);
  while (fgets(s, stringLen_g, f)!=NULL) { numberOfLines++; }
  
  numberOfLines -= 36;
  
  if (numberOfLines%5==0) return numberOfLines/5;
  else return -1;
}

//*************************************************************************************************
