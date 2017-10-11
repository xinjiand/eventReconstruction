//&&&&&&&&&&&&&&&&&&&&&&
//& Custom Headerfiles &
//&&&&&&&&&&&&&&&&&&&&&&
#include "3dvector.hh"
#include "crossSection.hh"
#include "deposit.hh"
#include "emProcesses.hh"
#include "gammaRay.hh"
#include "parsing.hh"
#include "randomNumberGenerator.hh"

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

//&&&&&&&&&&&&&&&&&&&
//& C++ Headerfiles &
//&&&&&&&&&&&&&&&&&&&
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

//&&&&&&&&&&&&&&&&&&&&&&&&&&
//& GNU Scientific Library &
//&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <gsl/gsl_rng.h>

//&&&&&&&&&&&&&
//& Constants &
//&&&&&&&&&&&&&
#define PI 3.14159265

#define RATIO_A2_TO_TAG_EVENTS 60             //For normal opertation this should be set to 60
#define PRINT_INTERACTION_PROBABILITIES 0     //
#define IN_LOADED_SCINTILLATOR_DENSITY 0.966  //In units of grams per cc
#define ALPHA_613_KEV_STATE 0.96              //The ratio of conversion electron to gamma emission
                                              //for the 613 keV state of Sn-115.
#define TRACKING_ENERGY_CUTOFF 6.4            //The energy below which we stop tracking the gamma.
                                              //It is set to the approximate value of the reconstruction
                                              //threshold for the mini-LENS reconstruction

//&&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Declarations &
//&&&&&&&&&&&&&&&&&&&&&&&&&

double pathLength(int genFlag, gsl_rng * rGen, double * sigmaTotal, gammaRay * g0);
int getInteraction(int genFlag, gsl_rng * rGen, double * probPelec, double * probComp,
  double * probPpro, gammaRay * g0);
void setDirection(int genFlag, gsl_rng * rGen, gammaRay * g0, double cTheta0);
int getGammaDeposits(gammaRay * g0, double * sigmaTotal, double * probPelec, double * probComp,
  double * probPpro, int genFlag, gsl_rng * rGen, vector<Deposit> * deposits);
  
//&&&&&&&&
//& Main &
//&&&&&&&&

int main ()
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Declare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  //Program control variables
  int genFlag, seedFlag;  //
  int seed;               //
  int numTagEvents;       //
  int numA2Events;        //
  int checkSuccess;       //
  double quenchFactor;    //
  gsl_rng *  rGen;        //
  const gsl_rng_type * T; //
  
  //Variables for reading data
  char string0[stringLen_g];  //
  char string1[stringLen_g];  //
  int nEntries;               //  

  //Cross-section and other data variables
  double * sigmaTotal;      //Tables for the total cross secion,
  double * photoElectric;   //the photoelectric cross section,
  double * compton;         //the Compton cross section, 
  double * pairProduction;  //and the pair production cross section.
  double * probPelec;       //Tables for the probability at each energy of the photoelectric effect, 
  double * probComp;        //Compton scattering, 
  double * probPpro;        //and of pair production
  double * quenchData;      //
  
  //Variables for tracking the gammas
  CVector r, v;                       //
  vector<Deposit> depList;            //
  double conversionElectronFraction;  //
  double cTheta, sTheta;              //
  double phi, cPhi, sPhi;             //
  
  //Outputing variables
  char tagDirName[stringLen_g];     //
  char a2DirName[stringLen_g];      //
  char outputFileName[stringLen_g]; //
  int nFiles;                       //
  
  //Files and directories
  FILE * in, * quench, * out; //
  DIR * dir;                  //


  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&

  //Get the input variables
  in = fopen("inputs/quench.in", "r");
  if (in==NULL)
  {
    printf("Could not open the input file! Exiting...\n");
    exit(1);
  }
  
  nEntries = getNumEntries(in);
  if (nEntries!=4)
  {
    printf("The input file is not properly formatted! Exiting...\n");
    exit(1); 
  }
  
  for (int i=0; i<4; i++)
  {
    getEntry(in, string0, stringLen_g);
    
    if      (i==0) { sscanf(string0, "%s %d", string1, &genFlag);       }
    else if (i==1) { sscanf(string0, "%s %d", string1, &seedFlag);      }
    else if (i==2) { sscanf(string0, "%s %lf", string1, &quenchFactor); }
    else           { sscanf(string0, "%s %d", string1, &numTagEvents);  }
  }
  fclose(in);
  
  if (genFlag<0 || genFlag>13) { genFlag = 1; }

  if (quenchFactor<=0.0) { quenchFactor = 1.0; }
  
  if (numTagEvents%1000!=0) { numTagEvents = 1000; }
  
  numA2Events = numTagEvents*RATIO_A2_TO_TAG_EVENTS;

  //Get the quenching data
  quenchData = new double [int(E_MAX)];
  for (int i=0; i<int(E_MAX); i++)
  {
    quenchData[i] = 1.0;
  }
  
  quench = fopen("data/quench.dat", "r");
  if (quench==NULL)
  {
    printf("Could not open the quenching data file! Exiting...\n");
    exit(1);
  }

  nEntries = getNumEntries(quench);
  for (int i=0; i<nEntries; i++)
  {
    double tempD;
    getEntry(quench, string0, stringLen_g);
    sscanf(string0, "%lf %lf", &tempD, &quenchData[i]);
    
    quenchData[i] = quenchData[i]*quenchFactor;
    if (quenchData[i]>1.0) { quenchData[i] = 1.0; }
  }
  fclose(quench);
  
  //Initialize the random number generator
  if (genFlag!=0)
  {
    if (genFlag==1)        T = gsl_rng_default;
    else if (genFlag==2)  T = gsl_rng_ranlxs0;
    else if (genFlag==3)  T = gsl_rng_ranlxs1;
    else if (genFlag==4)  T = gsl_rng_ranlxs2;
    else if (genFlag==5)  T = gsl_rng_ranlxd1;
    else if (genFlag==6)  T = gsl_rng_ranlxd2;
    else if (genFlag==7)  T = gsl_rng_ranlux;
    else if (genFlag==8)  T = gsl_rng_ranlux389;
    else if (genFlag==9)  T = gsl_rng_cmrg;
    else if (genFlag==10)  T = gsl_rng_mrg;
    else if (genFlag==11)  T = gsl_rng_taus;
    else if (genFlag==12)  T = gsl_rng_taus2;
    else if (genFlag==13)  T = gsl_rng_gfsr4;
    else 
    {
      printf("You entered an illegal value!\n");
      printf("Generator set to GSL default.\n");
      T = gsl_rng_default;
    }

    rGen = gsl_rng_alloc(T);
  }

  if (seedFlag==1) seed = time(NULL);
  else seed = 1013;                        //A prime number
  rngSeed(genFlag, seed, rGen);

  //Get the cross-section
  sigmaTotal     = new double [int(E_MAX)];
  photoElectric  = new double [int(E_MAX)];  
  compton        = new double [int(E_MAX)];
  pairProduction = new double [int(E_MAX)];

  probPelec = new double [int(E_MAX)];
  probComp  = new double [int(E_MAX)];
  probPpro  = new double [int(E_MAX)];

  checkSuccess = readCrossSection(107, sigmaTotal, photoElectric, compton, pairProduction, int(E_MAX));
  if (checkSuccess!=0)
  {
    printf("There was an error in reading the cross sections! Exiting...\n");
    exit(1);
  }

  for (int i=0; i<int(E_MAX); i++)
  {
    //Assign
    probPelec[i] = photoElectric[i]/sigmaTotal[i];
    probComp[i] = compton[i]/sigmaTotal[i];
    probPpro[i] = pairProduction[i]/sigmaTotal[i];

    //Normalize
    double tempD = probPelec[i] + probComp[i] + probPpro[i];
    probPelec[i] = probPelec[i]/tempD;
    probComp[i] = probComp[i]/tempD;
    probPpro[i] = probPpro[i]/tempD;
  }

  for (int i=0; i<int(E_MAX); i++)
  {
    sigmaTotal[i] = sigmaTotal[i]*IN_LOADED_SCINTILLATOR_DENSITY;
  }

  delete [] photoElectric;  
  delete [] compton;
  delete [] pairProduction;
  
  if (PRINT_INTERACTION_PROBABILITIES)
  {
    for (int i=0; i<int(E_MAX); i++)
    {
      printf("%d %lf %lf %lf\n", i+1, probPelec[i], probComp[i], probPpro[i]);
    }
  }
  
  //Make the output directories
  initString(tagDirName, stringLen_g);
  sprintf(tagDirName, "outputs/indium_tag_events");
  
  checkSuccess=0;
  dir = opendir(tagDirName);
  if (dir) { closedir(dir); } //Directory exists -> close and move on
  else if (errno==ENOENT)     //Directory does not exist -> create
  {
    checkSuccess = mkdir(tagDirName, 0777);
  }
  else                        //opendir() failed for some other reason
  {
    checkSuccess = 1;
  }

  if (checkSuccess)
  {
    printf("Error in opening indium tag event directory! Exiting...\n");
    exit(1);
  }
  strcat(tagDirName, "/");
  
  initString(a2DirName, stringLen_g);
  sprintf(a2DirName, "outputs/a2_background_events");

  checkSuccess=0;
  dir = opendir(a2DirName);    
  if (dir) { closedir(dir); } //Directory exists -> close and move on
  else if (errno==ENOENT)     //Directory does not exist -> create
  {
    checkSuccess = mkdir(a2DirName, 0777);
  }
  else                        //opendir() failed for some other reason
  {
    checkSuccess = 1;
  }

  if (checkSuccess)
  {
    printf("Error in opening A2 background event directory! Exiting...\n");
    exit(1);
  }
  strcat(a2DirName, "/");  
  
  //Initialize other variables
  conversionElectronFraction = ALPHA_613_KEV_STATE/(ALPHA_613_KEV_STATE + 1.0);
    
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Generate the Tag Events &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  nFiles = numTagEvents/1000;
  for (int i=0; i<nFiles; i++)
  {
    initString(outputFileName, stringLen_g);
    sprintf(outputFileName, "%sindium_tag_events_%09d.dat", tagDirName, i);
    
    out = fopen(outputFileName, "w");
    if (out==NULL)
    {
      printf("Could not open the file %s! Skipping to the next file...\n", outputFileName);
      continue;
    }
    
    for (int j=0; j<1000; j++)
    {
      //Print the required headers
      fprintf(out, "#Event\n2\n  %d\n  %d\n", i*1000+j, seed);
      fprintf(out, "#Geant4_properties\n5\n  12\n  612.9\n  0.0 0.0 0.0\n  0.0 0.0\n  0.0\n");
      fprintf(out, "#Energy_deposits\n");

      //Track the 115.6 keV gamma/conversion electron
      if (rng(genFlag, rGen)<=conversionElectronFraction) //It is a conversion electron
      {
        r.setxyz(0.0, 0.0, 0.0);
        Deposit dep;
        dep.setDeposit(r, 115.6, 0.0, 0);
        depList.push_back(dep);
      }
      else                                                //It is a gamma
      {
        //Initialize the gamma
        r.setxyz(0.0, 0.0, 0.0);
        v.setxyz(0.0, 0.0, C_LIGHT);
        
        gammaRay * g0 = new gammaRay(r, v, 115.6, 0.0, 0);
        
        //Track the gamma
        checkSuccess = getGammaDeposits(g0, sigmaTotal, probPelec, probComp, probPpro, genFlag, rGen, &depList);
        
        if (checkSuccess) {}
        delete g0;
      }
      
      //Track the 497.3 gamma
      r.setxyz(0.0, 0.0, 0.0);
      phi = 2.0*PI*rng(genFlag, rGen);
      cPhi = cos(phi);
      sPhi = sin(phi);        
      cTheta = 2.0*rng(genFlag, rGen) - 1.0;
      sTheta = sqrt(1.0 - cTheta*cTheta);
      v.setxyz(cPhi*sTheta, sPhi*sTheta, cTheta);
      v = v*C_LIGHT;
      
      gammaRay * g0 = new gammaRay(r, v, 497.3, 0.0, 0);
      
      //Track the gamma
      checkSuccess = getGammaDeposits(g0, sigmaTotal, probPelec, probComp, probPpro, genFlag, rGen, &depList);
      
      if (checkSuccess) { printf("There was a pair production event! Fire-up the debugger!\n"); }
      delete g0;
      
      //Apply quenching and output
      fprintf(out, "%d\n", int(depList.size()));
      for (unsigned int k=0; k<depList.size(); k++)
      {
        int qIndex = int(depList[k].E) - 1;
        depList[k].E = depList[k].E*quenchData[qIndex];
        fprintf(out, "  2 %lf %lf %lf %lf 0.0 %lf\n", depList[k].r.x*10.0, depList[k].r.y*10.0,
          depList[k].r.z*10.0, depList[k].t, depList[k].E/1000.0);
      }    
      
      //Clean-up for the next event
      depList.clear();
    } //End loop over events
    
    //Clean-up for the next file
    fclose(out);
  } //End loop over files  
    
    
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Generate the A2 Background Events &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  nFiles = numA2Events/1000;
  for (int i=0; i<nFiles; i++)
  {
    initString(outputFileName, stringLen_g);
    sprintf(outputFileName, "%sa2_background_events_%09d.dat", a2DirName, i);
    
    out = fopen(outputFileName, "w");
    if (out==NULL)
    {
      printf("Could not open the file %s! Skipping to the next file...\n", outputFileName);
      continue;
    }
  
    for (int j=0; j<1000; j++)
    {
      //Print the required headers
      fprintf(out, "#Event\n2\n  %d\n  %d\n", i*1000+j, seed);
      fprintf(out, "#Geant4_properties\n5\n  12\n  612.9\n  0.0 0.0 0.0\n  0.0 0.0\n  0.0\n");
      fprintf(out, "#Energy_deposits\n");
      
      //Initialize the 497.3 gamma
      r.setxyz(0.0, 0.0, 0.0);
      v.setxyz(0.0, 0.0, C_LIGHT);
      
      gammaRay * g0 = new gammaRay(r, v, 497.3, 0.0, 0);
      
      //Track the gamma
      checkSuccess = getGammaDeposits(g0, sigmaTotal, probPelec, probComp, probPpro, genFlag, rGen, &depList);
      
      if (checkSuccess) { printf("There was a pair production event! Fire-up the debugger!\n"); }
      delete g0;
      
      //Apply quenching and output
      fprintf(out, "%d\n", int(depList.size()));
      for (unsigned int k=0; k<depList.size(); k++)
      {
        int qIndex = int(depList[k].E) - 1;
        depList[k].E = depList[k].E*quenchData[qIndex];
        fprintf(out, "  2 %lf %lf %lf %lf 0.0 %lf\n", depList[k].r.x*10.0, depList[k].r.y*10.0,
          depList[k].r.z*10.0, depList[k].t, depList[k].E/1000.0);
      }    
          
      //Clean-up for the next event
      depList.clear();
          
    } //End loop over events
    
    //Clean-up for the next file
    fclose(out);    
    
  } //End loop over files
  
  
  //&&&&&&&&&&&&&&&
  //& End Program &
  //&&&&&&&&&&&&&&&
  
  //Free memory
  if (genFlag!=0) { gsl_rng_free(rGen); }
  delete [] quenchData;
  delete [] sigmaTotal;
  delete [] probPelec;
  delete [] probComp;
  delete [] probPpro;
  
  return 0;
}

//*************************************************************************************************

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&


double pathLength(int genFlag, gsl_rng * rGen, double sigmaTotal[], gammaRay * g0)
{
  /*This function samples the path-length of a particle travleing in a medium *
  *  with a given total macroscopic cross-section.                              */

  return -log(1.0 - rng(genFlag, rGen))/sigmaTotal[int(rint(g0->E))-1];
}

//*************************************************************************************************

int getInteraction(int genFlag, gsl_rng * rGen, double * probPelec, double * probComp,
  double * probPpro, gammaRay * g0)
{
  /*This function determines which interaction a particle undergoes */

  //&&&&&&&&&&&&&&&&&&&&&
  //& Declare Variables &
  //&&&&&&&&&&&&&&&&&&&&&

  double cdf1, cdf2;            //CDF values for the first and second interaction
  double rn=rng(genFlag, rGen);  //Random number uniform on [0,1)

  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&

  cdf1 = probPelec[int(rint(g0->E))];
  cdf2 = cdf1 + probComp[int(rint(g0->E))];

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Determine the Interaction &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  if (rn<=cdf1)        //The particle had interaction 1 (photoelectric)
  {
    return 1;
  }
  else if (rn<=cdf2)  //The particle had interaction 2 (Compton scattered)
  {
    return 2;
  }
  else                //The particle had interaction 3 (Pair production)
  {
    return 3;
  }
}

//*************************************************************************************************

void setDirection(int genFlag, gsl_rng * rGen, gammaRay * g0, double cTheta0)
{
  /*This function sets the new propagation direction */

  //&&&&&&&&&&&&&&&&&&&&&
  //& Declare Variables &
  //&&&&&&&&&&&&&&&&&&&&&

  double mag, u, v, w, s;      //Magnitude of the propagation vector, the direction cosines for the x, y, and z directions
  double phi0, cPhi0, sPhi0;  //The azmuthial scattering angle, its cosine and sine
  double sTheta0;              //Sine of the scattering angle

  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&

  u = g0->v.x;
  v = g0->v.y;
  w = g0->v.z;

  mag = u*u + v*v + w*w;
  mag = sqrt(mag);
  
  u = u/mag;
  v = v/mag;
  w = w/mag;

  s = sqrt(1.0 - w*w);

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Determine the Scattering Angles &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  //Isotropic scattering
  phi0 = 2.0*PI*rng(genFlag, rGen);
  cPhi0 = cos(phi0);
  sPhi0 = sin(phi0);

  sTheta0 = sqrt(1.0 - cTheta0*cTheta0);

  //&&&&&&&&&&&&&&&&&&
  //& Set the Values &
  //&&&&&&&&&&&&&&&&&&

  if (s!=0.0 && s!=-0.0)  //Scattering frame not aligned along the global z-axis
  {
    g0->v.setxyz(  ( u*cTheta0 - ((u*w*cPhi0 - v*sPhi0)/s)*sTheta0 )*mag,
                  ( v*cTheta0 - ((v*w*cPhi0 + u*sPhi0)/s)*sTheta0 )*mag,
                  ( s*sTheta0*cPhi0 + w*cTheta0 )*mag
                                                                        );
  }
  else                    //Scattering frame aligned along the global z-axis->phi=phi0 theta=theta0
  {
    g0->v.setxyz(cPhi0*sTheta0*mag, sPhi0*sTheta0*mag, cTheta0*mag);                                    
  }

  return;
}

//*************************************************************************************************

int getGammaDeposits(gammaRay * g0, double * sigmaTotal, double * probPelec, double * probComp,
  double * probPpro, int genFlag, gsl_rng * rGen, vector<Deposit> * deposits)
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Declare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  Deposit dep;          //
  double dt;            //
  double pathLen;       //
  int interactionType;  //
  double depE;          //
  double scattCosine;   //The cosine of the scattering angle
  int errorCnt=0;       //


  //&&&&&&&&&&&&&&&&&
  //& Tracking Loop &
  //&&&&&&&&&&&&&&&&&
  
  while (g0->E>TRACKING_ENERGY_CUTOFF)
  {
    //Determine the interaction point
    pathLen = pathLength(genFlag, rGen, sigmaTotal, g0);
    dt = pathLen/C_LIGHT;
    
    g0->r.x = g0->r.x + g0->v.x*dt;
    g0->r.y = g0->r.y + g0->v.y*dt;
    g0->r.z = g0->r.z + g0->v.z*dt;
    g0->tof += dt;
    
    //Determine the interaction
    interactionType = getInteraction(genFlag, rGen, probPelec, probComp, probPpro, g0);
    
    //Set the deposit and update the state of the gamma
    if (interactionType==1)       //Photoelectric effect
    {
      dep.setDeposit(g0->r, g0->E, g0->tof, g0->region);
      g0->E = 0.0;    
    }
    else if (interactionType==2)  //Compton scattered
    {
      depE = g0->E - comptonScatter(g0, genFlag, rGen, &scattCosine);
      dep.setDeposit(g0->r, depE, g0->tof, g0->region);    
      
      g0->E = g0->E - depE;
      setDirection(genFlag, rGen, g0, scattCosine);
    }
    else                          //Pair produced! -> return error
    {
      dep.setDeposit(g0->r, -1.0, g0->tof, g0->region);    
      g0->E = 0.0;
      errorCnt++;
    }
    
    deposits->push_back(dep);
  } //End tracking loop
  
  if (g0->E>0.0)  //If there is still energy, deposit it at the location of the gamma
  {
    dep.setDeposit(g0->r, g0->E, g0->tof, 0);
    deposits->push_back(dep);
  }


  //&&&&&&&
  //& End &
  //&&&&&&&

  return errorCnt;
}

//*************************************************************************************************
