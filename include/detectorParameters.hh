#ifndef DETECTOR_PARAMETERS_H
#define DETECTOR_PARAMETERS_H

//&&&&&&&&&&&&&
//& C Headers &
//&&&&&&&&&&&&&
#include <stdio.h>
#include <stdlib.h>

//&&&&&&&&&&&&&&&&&&
//& Custom Headers &
//&&&&&&&&&&&&&&&&&&
#include "parsing.hh"

//&&&&&&&&&&&&&
//& Constants &
//&&&&&&&&&&&&&

//Fundamental constants
#define PI 3.14159265


//Constants for calling the light transport
#define NCELL 15            //Number of cell along a side of NuLat, assumes a cubical detector
#define LATTICE_FLAG 1      //The type of lattice to use: 1->air, 2->perfluooctane, 3->Teflon FEP
#define MIRROR_FLAG 1       //Whether to mirror three sides (1) or not (0). If MIRROR_FLAG==0 and 
                            //GUIDE_FLAG==0 then the region around the PMT is mirrored
#define GUIDE_FLAG 1        //Determines the type of light gudie used for the lattice: 0->no guide, 
                            //1->one lattice channel into a guide, 2->two lattice channels into a 
                            //guide 3->three lattice channels into a guide. There are not 
#define PMT_COUPLING_FLAG 1 //Optically couples (1) or decouples (0) the PMT from the light guide. 
                            //The legacy behavior is to have the PMT optically coupled. 

//Constants for running the light transport
#define PMT_REFLECT 0.0           //The probability of reflecting off the PMT face, default is 0.0.
                                  //If this is set to 0.0 and the PMT is opticaly decoupled from the
                                  //light guide then the code assumes Fresnel reflections between
                                  //the guide index and an air gap.
#define TIME_BIN_SIZE 370         //The number of ps to bin together, default is 370
#define N_TIME_BIN_PRINT 110      //The number of time bins to print out, the maximum value is 200
#define CELL_DIM 6.35             //The length of a cell side in cm, default is 6.35 cm = 2.5 inch
#define PMT_RADIUS 2.3            //The radius of the PMT's photocathode in cm, default is 2.3 cm
                                  //from the Hamamatsu data sheet for the R10533 PMT
                                  
#define GAP_FILM_THICKNESS 1E-5   //The thickness of the film that makes up the gap. This is
                                  //included for legacy purposes and the default value of 1E-5
                                  //corresponds to essentially no film. NOT NOT USE 0, JUST USE THE
                                  //DEFAULT UNLESS SETTING LATTICE_FLAG TO 3.
                                  
#define GUIDE_LENGTH 3.5          //The total light guide length in cm, default is 3.5
#define SQUARE_LENGTH 0.5         //The length of the light guide that has a square cross section,
                                  //default is 0.5
                                  
#define REFLECT_COEFF 0.985       //The reflection coefficent of the mirrors the default is 0.985
#define GUIDE_FILM_THICKNESS 1E-5 //The thickness of the light guide's film, included from a legacy
                                  //application, default is 1E-5 which corresponds to basically no 
                                  //film. NOT NOT USE 0, JUST USE THE DEFAULT
                                  
#define PC_LIGHT_YIELD 12000.0    //Number of blue photons produced by 1 MeV energy deposited in PC,
                                  //default is 12000
                                  
#define QUANTUM_EFF 0.25          //The quantum efficenty of the photocathode, default is 0.25 from
                                  //Hamamatsu data sheet for the R10533 PMT
                                  
#define RELATIVE_YIELD 0.77      //The relative yield for EJ-254-1%, default is 0.77 from the
                                 //EJ-254 data sheet

#define INDEX_SCINT 1.58          //Refractive index of the scintillator, default is 1.58 from the
                                  //EJ-254 data sheet

#define INDEX_ACRYLIC 1.49        //Refractive index of acrylic, the guide material, default is 1.49
#define SCINT_ATTEN_LEN 400.0     //The attenuation length of the scintillator in cm, default is
                                  //400.0 cm. The best fit from the inital tests with a line of 15
                                  //cubes is 93 cm.

#define GUIDE_FILM_ATTEN 500.0    //The attenuation length in cm of the guide film, default is 500.0
                                  //see the note for GUIDE_FILM_THICKNESS

#define ACRYLIC_ATTEN_LEN 700.0   //The attenuation length in cm of the acrylic used for the guide,
                                  //default is 700.0 cm. 

#define GAP_1_ATTEN_LEN 500.0     //The attenuation length in cm of the LATTICE_FLAG==1 gap, default
                                  //is 500.0
                                  
#define N_FILMS_1 1               //The number of films that makeup the gap in LATTICE_FLAG==1,
                                  //included from legacy application, default is 1
                                  
#define INDEX_GAP_1 1.0           //The refractive index of the gap in LATTICE_FLAG==1, default is
                                  //1.0
                                  
#define GAP_WIDTH_1 0.0254        //The width of the gap in cm in LATTICE_FLAG==1, default is 0.0254

#define GAP_2_ATTEN_LEN 500.0     //The attenuation length in cm of the LATTICE_FLAG==2 gap, default
                                  //is 500.0
                                  
#define N_FILMS_2 1               //The number of films that makeup the gap in LATTICE_FLAG==2,
                                  //included from legacy application, default is 1
                                  
#define INDEX_GAP_2 1.225         //The refractive index of the gap in LATTICE_FLAG==2, default is
                                  //1.225
                                  
#define GAP_WIDTH_2 0.0254        //The width of the gap in cm in LATTICE_FLAG==2, default is 0.0254

#define GAP_3_ATTEN_LEN 0.5       //The attenuation length in cm of the LATTICE_FLAG==3 gap, default
                                  //is 0.5
                                  
#define N_FILMS_3 1               //The number of films that makeup the gap in LATTICE_FLAG==3,
                                  //included from legacy application, default is 1
                                  
#define INDEX_GAP_3 1.343         //The refractive index of the gap in LATTICE_FLAG==3, default is
                                  //1.343
                                  
#define GAP_WIDTH_3 0.0127        //The width of the gap in cm in LATTICE_FLAG==1, default is 0.0127


//Other flags

#define MAP_FLAG 1                //If using not standard instrumentation define these here, default
                                  //is 1. NOTE: This is currently not used because all
                                  //instrumentations are either half or full. 

//&&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Declarations &
//&&&&&&&&&&&&&&&&&&&&&&&&&

void printDetectorParameters(FILE * f);
int checkReconData(FILE * f);

//*************************************************************************************************

#endif
