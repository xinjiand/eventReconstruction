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

void readAndAddBlock(FILE * positron, FILE * neutron, FILE * out);

//*************************************************************************************************

//&&&&&&&&
//& Main &
//&&&&&&&&
int main(int argc, char * argv[])
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Delcare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  char positronFileName[stringLen_g];     //TO DO: Comment these
  char neutronFileName[stringLen_g];      //
  char outputFileName[stringLen_g];       //
  char sPositron[stringLen_g];            //
  char sNeutron[stringLen_g];             //
  char string0[stringLen_g];              //
  int contentFlag;                        //
  int positronEventNum, neutronEventNum;  //
  FILE * positron, * neutron, * out;      //

  
  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&
  
  //Get command-line argument
  if (argc!=2)
  {
    printf("Incorrect call! Wrong number of arguments!\n");
    exit(1);
  }
  
  contentFlag = atoi(argv[1]);

  //Open the files
  initString(positronFileName, stringLen_g);
  initString(neutronFileName, stringLen_g);
  initString(outputFileName, stringLen_g);
    
  if (contentFlag==0)
  {
    sprintf(positronFileName, "outputs/lightTransPositron.out");
    sprintf(neutronFileName, "outputs/lightTransNeutron.out");
    sprintf(outputFileName, "outputs/lightTrans.out");
  }
  else
  {
    sprintf(positronFileName, "outputs/reconstructionPositron.out");
    sprintf(neutronFileName, "outputs/reconstructionNeutron.out");  
    sprintf(outputFileName, "outputs/reconstruction.out");
  }
  
  positron = fopen(positronFileName, "r");
  if (positron==NULL)
  {
    printf("Could not open the file %s! Exiting...\n", positronFileName);
    exit(1);
  }
  
  neutron = fopen(neutronFileName, "r");
  if (neutron==NULL)
  {
    printf("Could not open the file %s! Exiting...\n", neutronFileName);
    exit(1);
  }
  
  //Open the output file
  out = fopen(outputFileName, "w");
  if (out==NULL)
  {
    printf("Could not open the file %s! Exiting...\n", outputFileName);
    exit(1);
  }

  
  //&&&&&&&&&&&&&&&&&&&&&
  //& Combine the Files &
  //&&&&&&&&&&&&&&&&&&&&&

  initString(sPositron, stringLen_g);
  initString(sNeutron, stringLen_g);  
  while (fgets(sPositron, stringLen_g, positron)!=NULL &&
    fgets(sNeutron, stringLen_g, neutron)!=NULL)
  {
    //Check event header
    if (strstr(sPositron, "Event")==NULL || strstr(sNeutron, "Event")==NULL)
    {
      printf("Error at the event header! Continue reading...\n");
      continue;
    }
    
    sscanf(sPositron, "%s %d", string0, &positronEventNum);
    sscanf(sNeutron, "%s %d", string0, &neutronEventNum);
    
    if (positronEventNum!=neutronEventNum)
    {
      printf("Event number do not match! Will continue anyway.\n");
    }
    
    fputs(sPositron, out);
    
    if (contentFlag==0) //Light transport
    {
      readAndAddBlock(positron, neutron, out);  //Read and write the hit PMTs and their charges
      readAndAddBlock(positron, neutron, out);  //Read and write the hit time bins
    }
    else                //Reconstruction
    {
      readAndAddBlock(positron, neutron, out);  //Read and write the hit cells
    }  
  } //End reading
  
  
  //&&&&&&&
  //& End &
  //&&&&&&&
  
  fclose(positron);
  fclose(neutron);
  fclose(out);
  
  return 0;
}

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&

//*************************************************************************************************

void readAndAddBlock(FILE * positron, FILE * neutron, FILE * out)
{
  //Declare needed variables
  char sPositron[stringLen_g];    //
  char sNeutron[stringLen_g];     //
  char string0[stringLen_g];      //
  int nHitPositron, nHitNeutron;  //

  //Initialize
  initString(sPositron, stringLen_g);
  initString(sNeutron, stringLen_g);
  
  fgets(sPositron, stringLen_g, positron);
  fgets(sNeutron, stringLen_g, neutron);

  sscanf(sPositron, "%d", &nHitPositron);
  sscanf(sNeutron, "%d", &nHitNeutron);
  
  initString(string0, stringLen_g);
  sprintf(string0, "%d %d %d\n", nHitPositron + nHitNeutron, nHitPositron, nHitNeutron);
  fputs(string0, out);

  //Read and write the positron block  
  for (int i=0; i<nHitPositron; i++)
  {
    initString(sPositron, stringLen_g);
    fgets(sPositron, stringLen_g, positron);
    fputs(sPositron, out);
  }

  //Read and write the neutron block  
  for (int i=0; i<nHitNeutron; i++)
  {
    initString(sNeutron, stringLen_g);    
    fgets(sNeutron, stringLen_g, neutron);
    fputs(sNeutron, out);
  }
}

//*************************************************************************************************
