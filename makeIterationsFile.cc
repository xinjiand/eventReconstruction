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

//&&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Declarations &
//&&&&&&&&&&&&&&&&&&&&&&&&&

int moveToEvent(FILE * event, int eventNumber);
void readEventFileBlock(FILE * event);

//*************************************************************************************************

//&&&&&&&&
//& Main &
//&&&&&&&&
int main(int argc, char * argv[])
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Delcare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  char eventFileName[stringLen_g];      //
  char outputFileName[stringLen_g];     //
  char string0[stringLen_g];            //
  char string1[stringLen_g];            //
  int eventNumber, checkSuccess;        //
  int foundEvent, currentEvent;         //
  int loopCnt, nLines;                  //
  FILE * event, * temp, * iter, * out;  //


  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&

  //Initialize strings 
  initString(eventFileName, stringLen_g); 
  initString(outputFileName, stringLen_g);
  initString(string0, stringLen_g);

  //Get the CL arguments
  if (argc<3)
  {
    printf("Incorrect call, wrong number of arguments!\n");
    printf("Call using: $ ./makeInterationsFile eventFileName eventNumber\n");
    exit(1);
  }
  
  copyString(eventFileName, argv[1], int(strlen(argv[1])) );
  
  eventNumber = atoi(argv[2]);
  
  
  //&&&&&&&&&&&&&&&&&
  //& Get the Event &
  //&&&&&&&&&&&&&&&&&
  
  event = fopen(eventFileName, "r");
  if (event==NULL)
  {
    printf("Could not open the file %s! Exiting...\n", eventFileName);
    exit(1);
  }
  
  temp = fopen(".iterationsTemp.dat", "w");
  if (temp==NULL)
  {
    printf("Could not open the file .iterationsTemp.dat! Exiting...\n");
    exit(1);
  }
  
  checkSuccess = moveToEvent(event, eventNumber);

  if (checkSuccess==0)
  {
    printf("There was an error in reading the event file!\n");
    fclose(event);
    fclose(temp);
    exit(1);
  }  
  
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, event);
  while (strstr(string0, "#Reconstruction")==NULL)
  {
    fputs(string0, temp);

    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, event);
  } //End event reading loop
  
  fclose(event);
  fclose(temp);
  
  
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Make the Iterations File &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  //Open the output file
  initString(outputFileName, stringLen_g);
  copyString(outputFileName, eventFileName, stringLen_g);
  removeFileExtension(outputFileName, 4);
  
  initString(string0, stringLen_g);
  sprintf(string0, "_%d.dat", eventNumber);
  strcat(outputFileName, string0);
  
  out = fopen(outputFileName, "w");
  if (out==NULL)
  {
    printf("Could not open the file %s! Exiting...\n", outputFileName);
    exit(1);
  }
  
  //Open the file with the iterations
  iter = fopen("outputs/reconIterations.out", "r");
  if (iter==NULL)
  {
    printf("Could not open the file outputs/reconIterations.out! Exiting...\n");
    fclose(out);
    exit(1);
  }
 
  //Find the desired event
  foundEvent = 0;
  while (foundEvent==0)
  {
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, iter);
    
    if (strstr(string0, "#Event")!=NULL)
    {
      sscanf(string0, "%s %d", string1, &currentEvent);
      
      if (currentEvent==eventNumber) { foundEvent = 1; }
    }
  } //End looking for the event
  
  //Open the temp file for reading
  temp = fopen(".iterationsTemp.dat", "r");
  if (temp==NULL)
  {
    printf("Could not open the file .iterationsTemp.dat! Exiting...\n");
    fclose(out);
    fclose(iter);
    exit(1);
  }
  
  //Make the iterations file
  printf("Copy the following into the event display input:\n");
  loopCnt = 0;
  initString(string0, stringLen_g);
  fgets(string0, stringLen_g, iter);
  while (strstr(string0, "#Event")==NULL)
  {
    //Printf the other headers
    initString(string1, stringLen_g);
    sprintf(string1, "#Event\n2\n  %d\n  0\n", loopCnt);
    fputs(string1, out);
    
    rewind(temp);
    initString(string1, stringLen_g);
    while (fgets(string1, stringLen_g, temp)!=NULL)
    {
      fputs(string1, out);
      initString(string1, stringLen_g);
    } //End the other heads loop
    
    //Write the "recon'ed" cells
    initString(string1, stringLen_g);
    sprintf(string1, "#Reconstruction\n");
    fputs(string1, out);
    fputs(string0, out);
        
    sscanf(string0, "%d", &nLines);
    for (int i=0; i<nLines; i++)
    {
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, iter);
      fputs(string0, out);
    }
    
    //Write event display input
    if (loopCnt==0) { printf("%s %d 1\n", eventFileName, eventNumber); }
    printf("%s %d 0\n", outputFileName, loopCnt);
    
    //Re-init for the next loop
    loopCnt++;
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, iter);
  }

  printf("%s %d 0\n", eventFileName, eventNumber); 
  //&&&&&&&
  //& End &
  //&&&&&&&
  
  fclose(out);
  fclose(temp);
  
  return 0;
}

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&

//*************************************************************************************************

int moveToEvent(FILE * event, int eventNumber)
{
  //Declare needed variables
  char string[stringLen_g];
  int successFlag=0;
  int nLines;
  int currentEventNumber = eventNumber + 1;
  
  //Start fresh
  rewind(event);
  
  //Loop until the correct event is found
  initString(string, stringLen_g);  
  while (currentEventNumber!=eventNumber && fgets(string, stringLen_g, event)!=NULL)
  {
    //Check the event header
    if (strstr(string, "#Event")==NULL)
    {
      return successFlag;
    }
    
    initString(string, stringLen_g);
    fgets(string, stringLen_g, event); 
    sscanf(string, "%d", &nLines);
    
    for (int i=0; i<nLines; i++)
    {
      initString(string, stringLen_g);
      fgets(string, stringLen_g, event);
      if (i==0) sscanf(string, "%d", &currentEventNumber);
    }  
    
    //Check if at destination
    if (currentEventNumber!=eventNumber)
    {
      for (int i=0; i<7; i++)
      {
        readEventFileBlock(event);
      }
    }
    else { successFlag=1; }
  
    initString(string, stringLen_g);  
  } //End search loop

  //End
  return successFlag;  
}

//*************************************************************************************************

void readEventFileBlock(FILE * event)
{
  char string[stringLen_g];
  int nLines;
  
  initString(string, stringLen_g);
  fgets(string, stringLen_g, event);
  
  initString(string, stringLen_g);    
  fgets(string, stringLen_g, event);
  sscanf(string, "%d", &nLines);

  for (int i=0; i<nLines; i++)
  {
    initString(string, stringLen_g);
    fgets(string, stringLen_g, event);
  }
}

//*************************************************************************************************
