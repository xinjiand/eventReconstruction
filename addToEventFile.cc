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

void readEventFileBlock(FILE * in, FILE * out, int ibdFlag);

//*************************************************************************************************

//&&&&&&&&
//& Main &
//&&&&&&&&
int main(int argc, char * argv[])
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Delcare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  char eventFileName[stringLen_g];      //TO DO: Go through and comment these
  char outputFileName[stringLen_g];     //
  char tempFileName[stringLen_g];       //
  char sEvent[stringLen_g];             //
  char sAdd[stringLen_g];               //
  char string0[stringLen_g];            //
  char * checkReturnChar;               //
  int contentFlag, wrongFormatFlag=0;   //
  int splitIndex, ibdFlag;              //
  int nLines, eventNum, checkEventNum;  //
  int nPositronLines, nNeutronLines;    //
  int nx, ny, nz;                       //
  FILE * add, * event, * temp;          //
  pmtTimeBin_t readLightTransLine;      //
  component_t readReconLine;            //
  
  
  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&

  //Initialize strings 
  initString(eventFileName, stringLen_g); 
  initString(sEvent, stringLen_g);
  initString(sAdd, stringLen_g);
  initString(string0, stringLen_g);

  //Get the CL arguments
  if (argc<3)
  {
    printf("Incorrect call, wrong number of arguments!\n");
    exit(1);
  }
  
  contentFlag = atoi(argv[1]);
  if (contentFlag!=0 && contentFlag!=1)
  {
    printf("Incorrect call. Use either 0 or 1. Exiting...\n");
    exit(1);
  }
  
  copyString(eventFileName, argv[2], int(strlen(argv[2])) );
  
  if (argc==4) splitIndex = atoi(argv[3]);
  else splitIndex=-1;
   
  //Open the needed files 
   
  //Open the correct output file
  initString(outputFileName, stringLen_g);
  
  if (contentFlag==0)
  {
    if (splitIndex==-1) sprintf(outputFileName, "outputs/lightTrans.out");
    else sprintf(outputFileName, "outputs/lightTrans%d.out", splitIndex);
  }
  else
  {
    if (splitIndex==-1) sprintf(outputFileName, "outputs/reconstruction.out");
    else sprintf(outputFileName, "outputs/reconstruction%d.out", splitIndex);  
  }
  
  add = fopen(outputFileName, "r");  
  if (add==NULL)
  {
    printf("Could not open the file %s! Exiting...\n", outputFileName);
    exit(1);
  }

  //Open the event file
  event = fopen(eventFileName, "r");
  if (event==NULL)
  {
    printf("Could not open the event file! Exiting...\n");
    exit(1);  
  }

  //Open the temp file
  initString(tempFileName, stringLen_g);
  if (splitIndex==-1) sprintf(tempFileName, ".temp.txt");
  else sprintf(tempFileName, ".temp%d.txt", splitIndex);
  
  temp = fopen(tempFileName, "w");
  if (temp==NULL)
  {
    printf("Could not open the temp file! Exiting...\n");
    exit(1);  
  }
  
  //Get other needed parameters
  ibdFlag=-1;
  if (strstr(eventFileName, "IBD")==NULL) { ibdFlag=0; }
  else { ibdFlag=1; }

  //&&&&&&&&&&&&&&&&&&&&&&&
  //& Copy the Event File &
  //&&&&&&&&&&&&&&&&&&&&&&&

  //TO DO: Use system calls to copy files

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
  temp = fopen(tempFileName, "r");
  if (temp==NULL)
  {
    printf("Could not open the temp file! Exiting...\n");
    exit(1);  
  }  

  initString(sEvent, stringLen_g);  
  while (fgets(sEvent, stringLen_g, temp)!=NULL)
  {
    if (strstr(sEvent, "#Event")!=NULL) //Check if at the event header
    {
      //Copy the first part of the event structure back to the event file
      fputs(sEvent, event);
      
      //Read the event number and the Geant4 seed
      initString(sEvent, stringLen_g);
      fgets(sEvent, stringLen_g, temp);
      fputs(sEvent, event);
      sscanf(sEvent, "%d", &nLines);
      for (int i=0; i<nLines; i++)
      {
        initString(sEvent, stringLen_g);
        fgets(sEvent, stringLen_g, temp);
        fputs(sEvent, event);
        if (i==0) sscanf(sEvent, "%d", &eventNum);
      }

      //Read the other blocks of the event structure
      readEventFileBlock(temp, event, 0);       //Geant4 properties
      readEventFileBlock(temp, event, ibdFlag); //Energy deposits
      readEventFileBlock(temp, event, 0);       //Light transport properties

      initString(sEvent, stringLen_g);  //Cells deposits
      fgets(sEvent, stringLen_g, temp);
      fputs(sEvent, event);
      
      initString(sEvent, stringLen_g);
      fgets(sEvent, stringLen_g, temp);
      fputs(sEvent, event);
      sscanf(sEvent, "%d", &nLines);
      
      for (int i=0; i<nLines; i++)
      {
        initString(sEvent, stringLen_g);
        fgets(sEvent, stringLen_g, temp);
        fputs(sEvent, event);
      }

      if (contentFlag==1)
      {
        readEventFileBlock(temp, event, ibdFlag); //PMT charges
        readEventFileBlock(temp, event, ibdFlag); //PMT waveforms    
      }

      //Check that I am adding the same event number in both add and event
      initString(sAdd, stringLen_g);
      fgets(sAdd, stringLen_g, add);
      
      if (strstr(sAdd, "Event")==NULL)
      {
        printf("Error in the formatting of %s!\n", outputFileName);

        wrongFormatFlag=1;
        break;
      }
      
      sscanf(sAdd, "%s %d", string0, &checkEventNum);

      if (eventNum!=checkEventNum)
      {
        printf("Error in event alignment in %s %d v. %d! Exiting...\n", outputFileName, eventNum,
          checkEventNum);
        wrongFormatFlag=1;
        break;
      }

      //Get and add the needed info            
      if (contentFlag==0) //Add the light transport 
      {
        //Get the number of hit PMTs
        initString(sAdd, stringLen_g);
        checkReturnChar = fgets(sAdd, stringLen_g, add);
        if (checkReturnChar) {}
        
        if (ibdFlag==0) sscanf(sAdd, "%d", &nLines);
        else sscanf(sAdd, "%d %d %d", &nLines, &nPositronLines, &nNeutronLines);
        
        //Add the header
        initString(string0, stringLen_g);
        if (ibdFlag==0) sprintf(string0, "#PMT_charges\n%d\n", nLines);
        else sprintf(string0, "#PMT_charges\n%d %d %d\n", nLines, nPositronLines, nNeutronLines);
        fputs(string0, event);
        
        //Read the data and add
        for (int i=0; i<nLines; i++)
        {
          //Read
          initString(sAdd, stringLen_g);
          checkReturnChar = fgets(sAdd, stringLen_g, add);
          
          //Transform from Bruce's coordinates to the normal coordinates
          readLightTransLine.bin=0;       
          sscanf(sAdd, "%d %d %d %lf", &readLightTransLine.nu, &readLightTransLine.nv, 
            &readLightTransLine.nw, &readLightTransLine.count);
          
          if (readLightTransLine.nw==1 || readLightTransLine.nw==2)       //An x-side
          {
            nx = (readLightTransLine.nw-1)*(NCELL+1);
            ny = readLightTransLine.nu;
            nz = readLightTransLine.nv;          
          }
          else if (readLightTransLine.nw==3 || readLightTransLine.nw==4)  //A y-side
          {
            nx = readLightTransLine.nu;
            ny = (readLightTransLine.nw-3)*(NCELL+1);
            nz = readLightTransLine.nv;          
          }
          else                                                                //A z-side
          {
            nx = readLightTransLine.nu;
            ny = readLightTransLine.nv;
            nz = (readLightTransLine.nw-5)*(NCELL+1);        
          }
          
          readLightTransLine.nu = nx;
          readLightTransLine.nv = ny;
          readLightTransLine.nw = nz;
          
          //Write
          initString(string0, stringLen_g);
          sprintf(string0, "  %d %d %d %lf\n", readLightTransLine.nu, readLightTransLine.nv, 
            readLightTransLine.nw, readLightTransLine.count);
          fputs(string0, event);
        } //End reading and adding the data


        //Get the hit pmt time bins
        initString(sAdd, stringLen_g);
        checkReturnChar = fgets(sAdd, stringLen_g, add);
        if (ibdFlag==0) sscanf(sAdd, "%d", &nLines);
        else sscanf(sAdd, "%d %d %d", &nLines, &nPositronLines, &nNeutronLines);
        
        //Add the header
        initString(string0, stringLen_g);
        initString(string0, stringLen_g);
        if (ibdFlag==0) sprintf(string0, "#PMT_waveforms\n%d\n", nLines);
        else sprintf(string0, "#PMT_waveforms\n%d %d %d\n", nLines, nPositronLines, nNeutronLines);        
        fputs(string0, event);
        
        //Read the data and add
        for (int i=0; i<nLines; i++)
        {
          //Read
          initString(sAdd, stringLen_g);
          checkReturnChar = fgets(sAdd, stringLen_g, add);
          
          //Transform from Bruce's coordinates to the normal coordinates
          sscanf(sAdd, "%d %d %d %d %lf", &readLightTransLine.nu, &readLightTransLine.nv, 
            &readLightTransLine.nw, &readLightTransLine.bin, &readLightTransLine.count);
          
          if (readLightTransLine.nw==1 || readLightTransLine.nw==2)       //An x-side
          {
            nx = (readLightTransLine.nw-1)*(NCELL+1);
            ny = readLightTransLine.nu;
            nz = readLightTransLine.nv;          
          }
          else if (readLightTransLine.nw==3 || readLightTransLine.nw==4)  //A y-side
          {
            nx = readLightTransLine.nu;
            ny = (readLightTransLine.nw-3)*(NCELL+1);
            nz = readLightTransLine.nv;          
          }
          else                                                            //A z-side
          {
            nx = readLightTransLine.nu;
            ny = readLightTransLine.nv;
            nz = (readLightTransLine.nw-5)*(NCELL+1);        
          }
          
          readLightTransLine.nu = nx;
          readLightTransLine.nv = ny;
          readLightTransLine.nw = nz;
          
          //Write
          initString(string0, stringLen_g);
          sprintf(string0, "  %d %d %d %d %lf\n", readLightTransLine.nu, readLightTransLine.nv, 
            readLightTransLine.nw, readLightTransLine.bin, readLightTransLine.count);
          fputs(string0, event);
        } //End reading and adding the data
      }
      else                //The reconstruction data
      {
        //Get the number of hit cells
        initString(sAdd, stringLen_g);
        checkReturnChar = fgets(sAdd, stringLen_g, add);
        if (ibdFlag==0) sscanf(sAdd, "%d", &nLines);
        else sscanf(sAdd, "%d %d %d", &nLines, &nPositronLines, &nNeutronLines);
        
        //Add the header
        initString(string0, stringLen_g);
        if (ibdFlag==0) sprintf(string0, "#Reconstruction\n%d\n", nLines );
        else sprintf(string0, "#Reconstruction\n%d %d %d\n", nLines, nPositronLines, nNeutronLines);
        fputs(string0, event);
        
        //Read the data and add
        for (int i=0; i<nLines; i++)
        {
          //Read
          initString(sAdd, stringLen_g);
          checkReturnChar = fgets(sAdd, stringLen_g, add);
          sscanf(sAdd, "%d %d %d %lf %lf", &readReconLine.nx, &readReconLine.ny, &readReconLine.nz,
            &readReconLine.t, &readReconLine.e);         
             
          //Write
          initString(string0, stringLen_g);
          sprintf(string0, "  %d %d %d %lf %lf\n", readReconLine.nx, readReconLine.ny,
            readReconLine.nz, readReconLine.t, readReconLine.e); 
          fputs(string0, event); 
        } //End reading and adding the data for the reconstruction data
      } //End reading and writing
    } //End check if at an event header in the event file
  } //End loop over the temp file 


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
    temp = fopen(tempFileName, "r");
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
  
  fclose(add);
  fclose(event);
  fclose(temp);
  
  return 0;
}

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&

//*************************************************************************************************

void readEventFileBlock(FILE * in, FILE * out, int ibdFlag)
{
  char string[stringLen_g];
  int nLines, nt1, nt2;
  
  initString(string, stringLen_g);
  fgets(string, stringLen_g, in);
  fputs(string, out);
  
  initString(string, stringLen_g);    
  fgets(string, stringLen_g, in);
  fputs(string, out);
  if (ibdFlag==0) sscanf(string, "%d", &nLines);
  else sscanf(string, "%d %d %d", &nLines, &nt1, &nt2);
  
  for (int i=0; i<nLines; i++)
  {
    initString(string, stringLen_g);
    fgets(string, stringLen_g, in);
    fputs(string, out);        
  }
}

//*************************************************************************************************
