#ifndef RECON_LIB_H
#define RECON_LIB_H

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

//&&&&&&&&&&&&&
//& Constants &
//&&&&&&&&&&&&&
#define N_SPLIT 4                         //Must be a positive integer and should not be larger than
                                          //the number of cores that are available on the machine.
                                          //Setting to >1 will split both the light transport
                                          //scripts and the reconstruction scripts
#define RECON_PROGRAM 1                   //0 means use Bruce's charge reconstruction (for backwards
                                          //compatibility) 1 means use the new charge/time
                                          //reconstruction
                                          
//TO DO: Implement using Bruce's code and doing a time reconsutrction
                          
#define TEST_FLAG_RECON 0                 //0: no test statements, 1: general test statements, 2:
                                          //timing minimization and waveform checking
#define TIME_INIT 1                       //Cell hit time initialization methods: 0 -> use first
                                          //non-zero bin and average, 1 -> use 50% leading edge and
                                          //average, 2 -> use maximum bin and use the minimum time
#define TIME_RECON 1                      //Flag used for backward compatibility in the time
                                          //reconstruction: 0->only changes of a bin size or greater
                                          //have an effect on the waveform, 1->sub bin sizes change
                                          //the waveform via rebinning.
#define PE_TO_ENERGY_CONVERSION_METHOD 1  //0 -> use a simple factor to convert from photoelectrons
                                          //to MeV's; 1 -> use a linear fit to convert from
                                          //photoelectrons to MeV's 
#define CELL_THRESHOLD 2                  //Flag used to control the threshold with which to output
                                          //cells 0->no threshold, 1->use the correction factor,
                                          //2->use the light transport threshold. 
#define E_MIN 0.01                        //The lowest value for a hit in the light transport. The
                                          //number is in MeV
#define STEP 0.1                          //The initial step size in the charge gradient search. 

//&&&&&&&&&&&&&&&&&&&
//& Data Structures &
//&&&&&&&&&&&&&&&&&&&
typedef struct component_t
{
  int nx, ny, nz;
  double e, t;
} component_t;

typedef struct deposit_t
{
  int id;
  double x, y, z, t, energy;
} deposit_t;

typedef struct detector_t
{
  double component[NCELL+2][NCELL+2][NCELL+2];
} detector_t;

typedef struct pmtTimeBin_t
{
  int nu, nv, nw, bin;
  double count;

} pmtTimeBin_t;

typedef struct vector_t
{
  double x, y, z;
} vector_t;

//&&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Declarations &
//&&&&&&&&&&&&&&&&&&&&&&&&&

void getEnergyScaleFactors(double * slope, double * intercept);
double getEnergy(double pe, double slope, double intercept);

void sortDeposits(deposit_t * list, int len, int sortFlag);
void sortComponents(component_t * c, int len, int sortFlag);
void clearDetector(detector_t * det);
void clearComponents(component_t * c, int len);
void copyDetector(detector_t * in, detector_t * out);

double getChiSqrCharge(double * primary, double * side, component_t * pmtList,
                        component_t * cellList, int pmtCnt, int cellCnt);
                        
double getChiSqrWaveform(double ** primary, double ** side, component_t * pmtList, double ** hitPMT,
                          component_t * cellList, int pmtCnt, int cellCnt);

double getChiSqrTeflon(double * primary, double ** side, component_t * pmtList,
                        component_t * cellList, int pmtCnt, int cellCnt);
  
void gradientSearchCharge(double * primary, double * side, component_t * pmtList,
                           component_t * cellList, int pmtCnt, int cellCnt);
                           
void gradientSearchWaveform(double ** primary, double ** side, component_t * pmtList, double ** hitPMT,
                             component_t * cellList, int pmtCnt, int cellCnt);

void gradientSearchTeflon(double * primary, double ** side, component_t * pmtList,
                           component_t * cellList, int pmtCnt, int cellCnt);

void absComponentEnergy(component_t * c, int len);
void resetTimeZero(component_t * c, int len);
void getWaveform(double ** templateWaveform, component_t cell, int sideIndex, component_t * pmtList,
                  double ** waveform);
                  
int getComponentIndex(component_t * c, int len, int nx, int ny, int nz);
void copyComponents(component_t * in, component_t * out, int len);
int getHighEnergyCellsPmtIndex(component_t * pmtList, component_t * cellList, int pmtCnt,
                                int cellCnt);
double getTotalCharge(component_t * c, int len);
void printDetectorPmts(detector_t * det, int chargeOrTimeFlag);

//**************************************************************************************************

#endif
