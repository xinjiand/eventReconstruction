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
#include "detectorParameters.hh"
#include "reconLib.hh"

//&&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Declarations &
//&&&&&&&&&&&&&&&&&&&&&&&&&

void readEventFileBlock(FILE * in, FILE * out);

//*************************************************************************************************

//&&&&&&&&
//& Main &
//&&&&&&&&
int main(int argc, char * argv[])
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Delcare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  //Events for reading and writing
  char eventFileName[stringLen_g];    //Name of event file that you want to remove content from
  char sEvent[stringLen_g];           //String for reading the event file
  char string0[stringLen_g];          //String for reading
  int headerFlag, wrongFormatFlag=0;  //Flag that id's the content to remove, error flag
  int nLines;                         //Used for counting in loops
  
  //Files
  FILE * event, * temp;               //The event and a temp file
  
  
  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&

  //Initialize strings 
  initString(eventFileName, stringLen_g); 
  initString(sEvent, stringLen_g);
  initString(string0, stringLen_g);  

  //Get the CL arguments
  if (argc!=3)
  {
    printf("Incorrect call, wrong number of arguments! Call using: headerFlag fileName\n");
    exit(1);
  }
  
  headerFlag = atoi(argv[1]);
  if (headerFlag<0 || headerFlag>4)
  {
    printf("Incorrect call. See the README for more information. Exiting...\n");
    exit(1);
  }
  
  copyString(eventFileName, argv[2], int(strlen(argv[2])) );
   
  //Open the event file
  event = fopen(eventFileName, "r");
  if (event==NULL)
  {
    printf("Could not open the event file! Exiting...\n");
    exit(1);  
  }

  //Open the temp file
  temp = fopen(".temp.txt", "w");
  if (temp==NULL)
  {
    printf("Could not open the temp file! Exiting...\n");
    exit(1);  
  }
  

  //&&&&&&&&&&&&&&&&&&&&&&&
  //& Copy the Event File &
  //&&&&&&&&&&&&&&&&&&&&&&&

  //TO DO: Use systems calls for this

  while (fgets(sEvent, stringLen_g, event)!=NULL)
  {
    fputs(sEvent, temp);
    initString(sEvent, stringLen_g);
  }

  fclose(temp);
  fclose(event);

  
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Add the Light Transport Output &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  //Re-open the event file
  event = fopen(eventFileName, "w");
  if (event==NULL)
  {
    printf("Could not open the event file! Exiting...\n");
    exit(1);  
  }

  //Re-open the temp file
  temp = fopen(".temp.txt", "r");
  if (temp==NULL)
  {
    printf("Could not open the temp file! Exiting...\n");
    exit(1);  
  }  

  initString(sEvent, stringLen_g);  
  while (fgets(sEvent, stringLen_g, temp)!=NULL)
  {
    //Event header
    if (strstr(sEvent, "#Event")==NULL)
    {
      wrongFormatFlag=1;
      printf("Error at the event header!\n");
      break;
    }
    fputs(sEvent, event);
    readEventFileBlock(temp, event);

    //Geant4 header
    initString(sEvent, stringLen_g);
    fgets(sEvent, stringLen_g, temp); 
    if (strstr(sEvent, "#Geant4_properties")==NULL)
    {
      wrongFormatFlag=1;
      printf("Error at the Geant4 header!\n");
      break;
    }
    fputs(sEvent, event);
    readEventFileBlock(temp, event);
    
    //Energy deposit header
    initString(sEvent, stringLen_g);
    fgets(sEvent, stringLen_g, temp);
    if (strstr(sEvent, "#Energy_deposits")==NULL)
    {
      wrongFormatFlag=1;
      printf("Error at the energy deposits header!\n");
      break;
    }
    fputs(sEvent, event);
    readEventFileBlock(temp, event);    
    
    //Repeat with only the headers that need to be kept
    for (int headerNum=4; headerNum>=0; headerNum--)
    {
      //Check the headers
      initString(sEvent, stringLen_g);
      fgets(sEvent, stringLen_g, temp);
      
      if (headerNum==4)       //Light transport header
      {
        if (strstr(sEvent, "#Light_transport_properties")==NULL)
        {
          wrongFormatFlag=1;
          printf("Error at the light transport header!\n");
          break;
        }      
      }
      else if (headerNum==3)  //Cell deposits header
      {
        if (strstr(sEvent, "#Cells_deposits")==NULL)
        {
          wrongFormatFlag=1;
          printf("Error at the cell deposits header!\n");
          break;
        }           
      }
      else if (headerNum==2)  //PMT charges header
      {
        if (strstr(sEvent, "#PMT_charges")==NULL)
        {
          wrongFormatFlag=1;
          printf("Error at the pmt charges header!\n");
          break;
        }       
      }
      else if (headerNum==1)  //PMT waveform header
      {
        if (strstr(sEvent, "#PMT_waveforms")==NULL)
        {
          wrongFormatFlag=1;
          printf("Error at the pmt waveform header!\n");
          break;
        }       
      }
      else                    //Reconstruction header
      {
        if (strstr(sEvent, "#Reconstruction")==NULL)
        {
          wrongFormatFlag=1;
          printf("Error at the reconstruction header!\n");
          break;
        }       
      }
      
      //Check if keeping or removing and do so
      if (headerFlag<headerNum) //Keep the block
      {
        fputs(sEvent, event);
        readEventFileBlock(temp, event);       
      }
      else                      //Remove the block
      {  
        initString(sEvent, stringLen_g);    
        fgets(sEvent, stringLen_g, temp);
        sscanf(sEvent, "%d", &nLines);
        for (int lineNum=0; lineNum<nLines; lineNum++)
        {
          initString(sEvent, stringLen_g);
          fgets(sEvent, stringLen_g, temp);
        }      
      } //End check keep/remove block      
    } //End loop over removable blocks

    if (wrongFormatFlag==1) break;
    
  } //End reading the temp file


  //&&&&&&&&&&&&&&&&&&&&&&&&
  //& Reset the Event File &
  //&&&&&&&&&&&&&&&&&&&&&&&&
  
  if (wrongFormatFlag==1)
  {
    //Re-open the event file
    fclose(event);
    event = fopen(eventFileName, "w");
    if (event==NULL)
    {
      printf("Could not open the event file! Exiting...\n");
      exit(1);  
    }

    //Re-open the temp file
    fclose(temp);
    temp = fopen(".temp.txt", "r");
    if (temp==NULL)
    {
      printf("Could not open the temp file! Exiting...\n");
      exit(1);  
    }
    
    //Echo the temp file back to the event file
    initString(string0, stringLen_g);
    while ( fgets(string0, stringLen_g, temp)!=NULL )
    {
      fputs(string0, event);
      initString(string0, stringLen_g);      
    }
  }


  //&&&&&&&
  //& End &
  //&&&&&&&
  
  fclose(event);
  fclose(temp);
  
  return 0;
}

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&

//*************************************************************************************************

void readEventFileBlock(FILE * in, FILE * out)
{
  char string[stringLen_g];
  int nLines;
  
  initString(string, stringLen_g);    
  fgets(string, stringLen_g, in);
  fputs(string, out);
  sscanf(string, "%d", &nLines);
  for (int i=0; i<nLines; i++)
  {
    initString(string, stringLen_g);
    fgets(string, stringLen_g, in);
    fputs(string, out);        
  }
}

//*************************************************************************************************
