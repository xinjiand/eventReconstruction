//&&&&&&&&&&&&&
//& C Headers &
//&&&&&&&&&&&&&
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//&&&&&&&&&&&&&&&
//& C++ Headers &
//&&&&&&&&&&&&&&&
#include <vector>
using namespace std;

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
int main(int argc, char * argv[])
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Delcare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  char basisFileName[stringLen_g];      //TO DO: Go through and comment these
  double energy, startTime;             //
  int nxTrue, nyTrue, nzTrue;           //
  int check, nx, ny, nz;                // 
  int nSidePMTs;                        //
  pmtTimeBin_t readLightTransLine;      //
  vector<pmtTimeBin_t> lightTransEntry; //
  pmtTimeBin_t * primary;               //
  pmtTimeBin_t ** side;                 //
  FILE * basis, * light;                //
  
  
  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&

  //Open the basis data file  
  initString(basisFileName, stringLen_g);
  sprintf(basisFileName, "reconstructionData/reconBasisData_%d_%d%d%d_%d_%d.dat", NCELL,
    LATTICE_FLAG, MIRROR_FLAG, GUIDE_FLAG, TIME_BIN_SIZE, N_TIME_BIN_PRINT);
    
  basis = fopen(basisFileName, "a");
  if (basis==NULL)
  {
    printf("Could not open the basis file! Exiting...\n");
    exit(1);
  }
  
  //Open the light transport output file
  light = fopen("outputs/trans-time.dat", "r");
  if (light==NULL)
  {
    printf("Could not open the light transport output file! Exiting...\n");
    exit(1);
  }
  
  //Get the command line arguments
  if (argc!=6)
  {
    printf("Incorrect call! Exiting...\n");
    exit(1);
  }
  
  energy = atof(argv[1]);
  startTime = atof(argv[2]);
  nxTrue = atoi(argv[3]);
  nyTrue = atoi(argv[4]);
  nzTrue = atoi(argv[5]);
  
  if (startTime!=0.0)
  {
    printf("start times not equal to zero are not implemented yet.\n");
  }
 
  //Initialize the arrays to hold the waveforms
  primary = new pmtTimeBin_t [N_TIME_BIN_PRINT];
  if (LATTICE_FLAG==1) { nSidePMTs = 1; }
  else { nSidePMTs = NCELL-1; }
  
  side = new pmtTimeBin_t* [nSidePMTs];
  for (int i=0; i<nSidePMTs; i++)
  {
    side[i] = new pmtTimeBin_t [N_TIME_BIN_PRINT];
  }  

  //Initalize the tempate waveform
  for (int i=0; i<nSidePMTs; i++)
  {
    for (int j=0; j<N_TIME_BIN_PRINT; j++)
    {
      if (i==0)
      {
        primary[j].nu = nxTrue;
        primary[j].nv = nyTrue;
        primary[j].nw = 0;
        primary[j].bin = j+1;
        primary[j].count = 0.0;
      }
      
      side[i][j].nu = nxTrue;
      side[i][j].nv = nyTrue;
      side[i][j].nw = 0;
      side[i][j].bin = j+1;
      side[i][j].count = 0.0;
    }
  } 
 
  
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Read the Light Transport Output &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  check = fscanf(light, "%d %d %d %d %lf", &readLightTransLine.nu, &readLightTransLine.nv,
    &readLightTransLine.nw, &readLightTransLine.bin, &readLightTransLine.count);
  lightTransEntry.push_back(readLightTransLine);
  
  while (check!=EOF)
  {
    check = fscanf(light, "%d %d %d %d %lf", &readLightTransLine.nu, &readLightTransLine.nv,
      &readLightTransLine.nw, &readLightTransLine.bin, &readLightTransLine.count);
    lightTransEntry.push_back(readLightTransLine);  
  }
  fclose(light);
  
  
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Output the Template Waveforms &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  //Transform from Bruce's coordinates to the normal coordinates  
  for (int i=0; i<int(lightTransEntry.size()); i++)
  {
    if (lightTransEntry[i].nw==1 || lightTransEntry[i].nw==2)       //An x-side
    {
      nx = (lightTransEntry[i].nw-1)*(NCELL+1);
      ny = lightTransEntry[i].nu;
      nz = lightTransEntry[i].nv;          
    }
    else if (lightTransEntry[i].nw==3 || lightTransEntry[i].nw==4)  //A y-side
    {
      nx = lightTransEntry[i].nu;
      ny = (lightTransEntry[i].nw-3)*(NCELL+1);
      nz = lightTransEntry[i].nv;          
    }
    else                                                            //A z-side
    {
      nx = lightTransEntry[i].nu;
      ny = lightTransEntry[i].nv;
      nz = (lightTransEntry[i].nw-5)*(NCELL+1);        
    }
    
    lightTransEntry[i].nu = nx;
    lightTransEntry[i].nv = ny;
    lightTransEntry[i].nw = nz;
    lightTransEntry[i].count = lightTransEntry[i].count;
  }
  
  //Get the waveforms, take the z=0 side for the data
  for (int i=0; i<int(lightTransEntry.size()); i++) //Loop over all the outputted bins
  {
    for (int j=0; j<N_TIME_BIN_PRINT; j++)          //Loop over all the bins of the primary and 
    {                                               //side pmts
      if ( lightTransEntry[i].nu==primary[j].nu && lightTransEntry[i].nv==primary[j].nv //Primary pmt
        && lightTransEntry[i].nw==primary[j].nw && lightTransEntry[i].bin==primary[j].bin)
      {
        primary[j].count += lightTransEntry[i].count;
      }
      
      for (int k=0; k<nSidePMTs; k++)              //Loop over all possible side band pmts
      {
        if (abs(lightTransEntry[i].nu-side[k][j].nu)==k+1 && lightTransEntry[i].nv==side[k][j].nv
          && lightTransEntry[i].nw==side[k][j].nw && lightTransEntry[i].bin==side[k][j].bin)
        {
          side[k][j].count += lightTransEntry[i].count;
        }
        else if (lightTransEntry[i].nu==side[k][j].nu && abs(lightTransEntry[i].nv-side[k][j].nv)==k+1
          && lightTransEntry[i].nw==side[k][j].nw && lightTransEntry[i].bin==side[k][j].bin)
        {
          side[k][j].count += lightTransEntry[i].count;
        }      
      } //End loop over the side PMTs
    } //End the loop over time bins
  } //End loop over vector entries
  lightTransEntry.clear();
  
  //Normalize the waveforms
  for (int i=0; i<nSidePMTs; i++)
  {
    for (int j=0; j<N_TIME_BIN_PRINT; j++)
    {
      if (i==0)
      {
        primary[j].count = primary[j].count/energy;
      }
      
      side[i][j].count = side[i][j].count/(2.0*energy);
    }
  } 

  fprintf(basis, "%d Cell(s) Away:\n", nzTrue);
  fprintf(basis, "  Primary:\n");
  for (int i=0; i<N_TIME_BIN_PRINT; i++)
  {
    fprintf(basis, "    %d %lf\n", primary[i].bin, primary[i].count );
  }
  
  for (int i=0; i<nSidePMTs; i++)
  {
    if (nSidePMTs==1) { fprintf(basis, "  Side:\n"); }
    else { fprintf(basis, "  Side %d:\n", i+1); }

    for (int j=0; j<N_TIME_BIN_PRINT; j++)
    {
      fprintf(basis, "    %d %lf\n", side[i][j].bin, side[i][j].count );
    }    
  }
  fclose(basis);

  
  //&&&&&&&
  //& End &
  //&&&&&&&
  
  delete [] primary;
  for (int i=0; i<nSidePMTs; i++) { delete [] side[i]; }
  delete [] side;
  
  return 0;
}

//*************************************************************************************************
