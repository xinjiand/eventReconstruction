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
#include "detectorParameters.hh"
#include "parsing.hh"
#include "reconLib.hh"

//*************************************************************************************************

//&&&&&&&&
//& Main &
//&&&&&&&&
int main ()
{
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	//& Declare Needed Variables &
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	char string0[stringLen_g];				                    //TO DO: Go through and comment on these
	char scriptName[stringLen_g];				                  //
  char eventFileName[stringLen_g];                      //
  char reconProgramName[stringLen_g];                   //
  char ltSplitScriptFileName[N_SPLIT][stringLen_g];     //  
  char reconSplitScriptFileName[N_SPLIT][stringLen_g];  //
	int nEntries, scriptIndex;                            //
	FILE * list, * ltScript, * reconScript;               //
	FILE * ltSplitScript[N_SPLIT];                        //
	FILE * reconSplitScript[N_SPLIT];                     //


	//&&&&&&&&&&&&&&
	//& Initialize &
	//&&&&&&&&&&&&&&
	
	//Set the correct reconstruction program
	initString(reconProgramName, stringLen_g);
	if (LATTICE_FLAG==1) { sprintf(reconProgramName, "reconstruct"); }
	else { sprintf(reconProgramName, "reconstructTeflon"); }
	
  //Open the input
	list = fopen("./inputs/scriptGenerator.in", "r");
	if (list==NULL)
	{
		printf("Cannot open the file! Exiting...\n");
		exit(1); 
	}
	
	//Open the split scripts
  if (N_SPLIT>1)
  {
    //For the light transport
    for (int i=0; i<N_SPLIT; i++)
    {
      initString(ltSplitScriptFileName[i], stringLen_g);
      sprintf(ltSplitScriptFileName[i], "runLightTransport%d.sh", i);
      
      ltSplitScript[i] = fopen(ltSplitScriptFileName[i], "w");
      if (ltSplitScript[i]==NULL)
      {
        printf("Could not open the file %s! Exiting...\n", ltSplitScriptFileName[i]);
        exit(1);
      }
    }
    
    //For the reconstruction
    for (int i=0; i<N_SPLIT; i++)
    {
      initString(reconSplitScriptFileName[i], stringLen_g);
      sprintf(reconSplitScriptFileName[i], "runReconstruction%d.sh", i);
      
      reconSplitScript[i] = fopen(reconSplitScriptFileName[i], "w");
      if (reconSplitScript[i]==NULL)
      {
        printf("Could not open the file %s! Exiting...\n", reconSplitScriptFileName[i]);
        exit(1);
      }
    }
  }

  //Open the light transport script
	ltScript = fopen("runLightTransport.sh", "w");
	if (ltScript==NULL)
	{
		printf("Cannot open the light transport script! Exiting...\n");
		exit(1); 
	}

  fprintf(ltScript, "#This is a bash script for running the light transport code\n");
	fprintf(ltScript, "echo \"Start Time:\" ;date\n");
  if (N_SPLIT>1)
	{
	  for (int i=0; i<N_SPLIT; i++)
	  {
	    fprintf(ltScript, "./%s &\n", ltSplitScriptFileName[i]);
	  }
	  fprintf(ltScript,
	    "echo \"Wait until the light transport finishes before running the reconstruction.\"\n");
	  fprintf(ltScript, "exit\n");
	}
	
  //Open the reconstruction script
	reconScript = fopen("runReconstruction.sh", "w");
	if (reconScript==NULL)
	{
		printf("Cannot open the reconstruction script! Exiting...\n");
		exit(1); 
	}

	fprintf(reconScript, "#This is a bash script for running the reconstruction code\n");
	fprintf(reconScript, "echo \"Start Time:\" ;date\n");
  if (N_SPLIT>1)
	{
	  for (int i=0; i<N_SPLIT; i++)
	  {
	    fprintf(reconScript, "./%s &\n", reconSplitScriptFileName[i]);
	  }
	  fprintf(reconScript,
	    "echo \"Wait until the reconstruction finishes before running any analysis.\"\n");
	  fprintf(reconScript, "exit\n");
	}

	
	//&&&&&&&&&&&&&&&&&&
	//& Read the Input &
	//&&&&&&&&&&&&&&&&&&

	nEntries = getNumEntries(list);	
	for (int i=0; i<nEntries; i++)
	{
	  //Get the script and event file names
    getEntry(list, string0, stringLen_g);
    copyString(scriptName, string0, int(strlen(string0)));
    
    initString(eventFileName, stringLen_g);
    sprintf(eventFileName, "eventFiles/");
    
		removePath(string0, 1);
		strcat(eventFileName, string0);
		
		removeFileExtension(eventFileName, 6);
		strcat(eventFileName, "Event.dat");
		
		//Determine if an IBD event file
    int ibdFlag=-1;
    if (strstr(eventFileName, "IBD")==NULL) { ibdFlag=0; }
    else { ibdFlag=1; }  

	  //Generate the light transport script
	  if (N_SPLIT==1)
	  {
		  fprintf(ltScript, "./%s\n", scriptName);    
      fprintf(ltScript, "./addToEventFile 0 %s\n", eventFileName);
    }
	  else
	  {
	    scriptIndex = i - int(i/N_SPLIT)*N_SPLIT;
	    fprintf(ltSplitScript[scriptIndex], "./%s\n", scriptName);
	    fprintf(ltSplitScript[scriptIndex], "./addToEventFile 0 %s %d\n", eventFileName, scriptIndex);
	  }

    //Generate the reconstruction script
    if (RECON_PROGRAM==0) //Use the original Fortran code
    {
      fprintf(reconScript, "echo \"Nothing to do here, which is normal!\"\n");
      fprintf(ltScript, "./addToEventFile 1 %s\n", eventFileName);
    }
    else                  //Use the charge and timing reconstruction code
    { 
      if (N_SPLIT==1)               
      {
        if (ibdFlag==0)
        {
          fprintf(reconScript, "./%s %s\n", reconProgramName, eventFileName);
        }
        else
        {
          //Get the file names for the split event files
          char positronEventFileName[stringLen_g];
          char neutronEventFileName[stringLen_g];
          
          initString(positronEventFileName, stringLen_g);
          initString(neutronEventFileName, stringLen_g);
          
          sprintf(positronEventFileName, "%s", eventFileName);
          sprintf(neutronEventFileName, "%s", eventFileName);
          
          removeFileExtension(positronEventFileName, 4);
          removeFileExtension(neutronEventFileName, 4);
          
          strcat(positronEventFileName, "_Positron.dat");
          strcat(neutronEventFileName, "_Neutron.dat");
          
          //Print the needed commands
          fprintf(reconScript, "./splitIbdForReconstruction %s %s %s\n", eventFileName,
            positronEventFileName, neutronEventFileName);
          fprintf(reconScript, "./%s %s\n", reconProgramName, positronEventFileName);
          fprintf(reconScript, "./%s %s\n", reconProgramName, neutronEventFileName);  
          fprintf(reconScript, "./combineIBD_Events 1\n");
          //fprintf(reconScript, "rm %s\n", positronEventFileName);
          //fprintf(reconScript, "rm %s\n", neutronEventFileName);
        }
        
        fprintf(reconScript, "./addToEventFile 1 %s\n", eventFileName);
      }
      else
      {
	      scriptIndex = i - int(i/N_SPLIT)*N_SPLIT;
	      fprintf(reconSplitScript[scriptIndex], "./%s %s %d\n", reconProgramName, eventFileName, scriptIndex);
	      fprintf(reconSplitScript[scriptIndex], "./addToEventFile 1 %s %d\n", eventFileName, scriptIndex);
	      
	      //TO DO: Implement this script splitting for IBD events
      }
    }
    
	}	//End loop over entries
	
	if (N_SPLIT==1)
	{
	  fprintf(ltScript, "echo \"End Time:\" ;date\n");	
  	fprintf(reconScript, "echo \"End Time:\" ;date\n");
	}
	else
	{
    for (int i=0; i<N_SPLIT; i++)
    {
  	  fprintf(ltSplitScript[i], "echo \"%s End Time:\" ;date\n", ltSplitScriptFileName[i]);
  	  fprintf(reconSplitScript[i], "echo \"%s End Time:\" ;date\n", reconSplitScriptFileName[i]);
    }	
	}

	
	//&&&&&&&
	//& End &
	//&&&&&&&

  //Close files
	fclose(list);
	fclose(ltScript);
  fclose(reconScript);
	if (N_SPLIT>1)
	{
    for (int i=0; i<N_SPLIT; i++)
    {
      fclose(ltSplitScript[i]);
      fclose(reconSplitScript[i]);      
    }	
	}

  //Make	the scripts executable
	chmod("runLightTransport.sh", 0777);
  chmod("runReconstruction.sh", 0777);
	if (N_SPLIT>1)
	{
    for (int i=0; i<N_SPLIT; i++)
    {
      chmod(ltSplitScriptFileName[i], 0777);
      chmod(reconSplitScriptFileName[i], 0777);
    }	
	}

	return 0;
}

//*************************************************************************************************
