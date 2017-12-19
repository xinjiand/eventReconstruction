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

void readAndWriteEventFileBlock(FILE * in, FILE * out1, FILE * out2, int skipHeaderName, int ibdFlag);

//*************************************************************************************************

//&&&&&&&&
//& Main &
//&&&&&&&&
int main(int argc, char * argv[])
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Delcare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  char eventFileName[stringLen_g];            //TO DO: Go through and comment these
  char positronEventFileName[stringLen_g];    //
  char neutronEventFileName[stringLen_g];     //
  char sEvent[stringLen_g];                   //
  int nLines, nPositronLines, nNeutronLines;  //
  FILE * event, * positron, * neutron;        //  
  
  
  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&

  //Initialize the file names
  initString(eventFileName, stringLen_g);
  initString(positronEventFileName, stringLen_g);
  initString(neutronEventFileName, stringLen_g);

  //Get the CL arguments
  if (argc<3)
  {
    printf("Incorrect call, wrong number of arguments!\n");
    exit(1);
  }
  
  copyString(eventFileName, argv[1], int(strlen(argv[1])) );
  copyString(positronEventFileName, argv[2], int(strlen(argv[2])) );
  copyString(neutronEventFileName, argv[3], int(strlen(argv[3])) );
  
  //Open the needed files 
  event = fopen(eventFileName, "r");
  if (event==NULL)
  {
    printf("Could not open the file %s! Exiting...\n", eventFileName);
    exit(1);  
  } 

  positron = fopen(positronEventFileName, "w");
  if (positron==NULL)
  {
    printf("Could not open the file %s! Exiting...\n", positronEventFileName);
    fclose(event);
    exit(1);  
  }
  
  neutron = fopen(neutronEventFileName, "w");
  if (neutron==NULL)
  {
    printf("Could not open the file %s! Exiting...\n", neutronEventFileName);
    fclose(event);
    fclose(positron);
    exit(1);  
  }
  

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Generate the New Event Files &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  initString(sEvent, stringLen_g);
  while (fgets(sEvent, stringLen_g, event)!=NULL)
  {
    //Check the event header
    if (strstr(sEvent, "#Event")==NULL) 
    {
      printf("Incorrect format! Exiting loop!\n");
      break;
    }

    //Write the event header
    fputs(sEvent, positron);
    fputs(sEvent, neutron);
    readAndWriteEventFileBlock(event, positron, neutron, 1, 0); 

    //Read and write the pre-light transport headers
    readAndWriteEventFileBlock(event, positron, neutron, 0, 0); //Geant4 properties
    readAndWriteEventFileBlock(event, positron, neutron, 0, 1); //Energy deposits
    readAndWriteEventFileBlock(event, positron, neutron, 0, 0); //Light transport properties
    readAndWriteEventFileBlock(event, positron, neutron, 0, 1); //Cell deposits
    
    //Split the light transport results
    
    //The PMT charges
    initString(sEvent, stringLen_g);
    fgets(sEvent, stringLen_g, event);
    fputs(sEvent, positron);
    fputs(sEvent, neutron);
    
    initString(sEvent, stringLen_g); 
    fgets(sEvent, stringLen_g, event);
    sscanf(sEvent, "%d %d %d", &nLines, &nPositronLines, &nNeutronLines);
    
    //The positron file
    initString(sEvent, stringLen_g);
    sprintf(sEvent, "%d\n", nPositronLines);
    fputs(sEvent, positron);
    
    for (int i=0; i<nPositronLines; i++)
    {
      initString(sEvent, stringLen_g);
      fgets(sEvent, stringLen_g, event);
      fputs(sEvent, positron);    
    }
    
    //The neutron file
    initString(sEvent, stringLen_g);
    sprintf(sEvent, "%d\n", nNeutronLines);
    fputs(sEvent, neutron);
    
    for (int i=0; i<nNeutronLines; i++)
    {
      initString(sEvent, stringLen_g);
      fgets(sEvent, stringLen_g, event);
      fputs(sEvent, neutron);    
    }
    
    //The PMT waveforms
    initString(sEvent, stringLen_g);
    fgets(sEvent, stringLen_g, event);
    fputs(sEvent, positron);
    fputs(sEvent, neutron);
    
    initString(sEvent, stringLen_g); 
    fgets(sEvent, stringLen_g, event);
    sscanf(sEvent, "%d %d %d", &nLines, &nPositronLines, &nNeutronLines);
    
    //The positron file
    initString(sEvent, stringLen_g);
    sprintf(sEvent, "%d\n", nPositronLines);
    fputs(sEvent, positron);
    
    for (int i=0; i<nPositronLines; i++)
    {
      initString(sEvent, stringLen_g);
      fgets(sEvent, stringLen_g, event);
      fputs(sEvent, positron);    
    }
    
    //The neutron file
    initString(sEvent, stringLen_g);
    sprintf(sEvent, "%d\n", nNeutronLines);
    fputs(sEvent, neutron);
    
    for (int i=0; i<nNeutronLines; i++)
    {
      initString(sEvent, stringLen_g);
      fgets(sEvent, stringLen_g, event);
      fputs(sEvent, neutron);    
    }
  }

  //&&&&&&&
  //& End &
  //&&&&&&&
  
  fclose(event);
  fclose(positron);
  fclose(neutron);
  
  return 0;
}

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&

//*************************************************************************************************

void readAndWriteEventFileBlock(FILE * in, FILE * out1, FILE * out2, int skipHeaderName, int ibdFlag)
{
  char string[stringLen_g];
  int nLines, nt1, nt2;

  if (skipHeaderName==0)
  {
    initString(string, stringLen_g);
    fgets(string, stringLen_g, in);
    fputs(string, out1);
    fputs(string, out2);
  }
    
  initString(string, stringLen_g);    
  fgets(string, stringLen_g, in);
  fputs(string, out1);
  fputs(string, out2);
  if (ibdFlag==0) sscanf(string, "%d", &nLines);
  else sscanf(string, "%d %d %d", &nLines, &nt1, &nt2);
  
  for (int i=0; i<nLines; i++)
  {
    initString(string, stringLen_g);
    fgets(string, stringLen_g, in);
    fputs(string, out1);
    fputs(string, out2);
  }
}

//*************************************************************************************************
