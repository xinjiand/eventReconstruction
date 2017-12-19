#include "detectorParameters.hh"

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&

//*************************************************************************************************

void printDetectorParameters(FILE * f)
{
  fprintf(f, "%lf     !pmtreflect    reflection from pmt faces\n", PMT_REFLECT);
  fprintf(f, "%d          !ntbin         how many ps to bin together (integer)\n", TIME_BIN_SIZE); 
  fprintf(f, "%d          !ntprint       how many time bins to print out (integer)\n", N_TIME_BIN_PRINT);
  fprintf(f, "%lf     !celldim       corresponds to 2.5 inches\n", CELL_DIM);
  fprintf(f, "%lf     !pmtr          estimated photo-cathode of 2 inch tube\n", PMT_RADIUS);                                
  fprintf(f, "%lf     !tfilm         corresponds to basically no film\n", GAP_FILM_THICKNESS);                                     
  fprintf(f, "%lf     !glength       overall length of guide (cm)\n", GUIDE_LENGTH);
  fprintf(f,
    "%lf     !gztrans       detector side to where one transitions from rectangular to condensing guide\n",
      SQUARE_LENGTH);
  fprintf(f, "%lf     !agrcoef       reflection coefficient VM200 film\n", REFLECT_COEFF);
  fprintf(f, "%lf     !tgfilm        thickness of guide film (cm)\n", GUIDE_FILM_THICKNESS); 
  fprintf(f, "%lf !pcyield       photons per MeV\n", PC_LIGHT_YIELD);                                 
  fprintf(f, "%lf     !qe            quantum efficiency\n", QUANTUM_EFF);
  fprintf(f, "%lf     !relativeyield photon yield relative to PC\n", RELATIVE_YIELD);
  fprintf(f, "%lf     !rscint        index of plastic scintillator\n", INDEX_SCINT);
  fprintf(f, "%lf     !rshield       index of refraction for guide material (acrylic)\n", 
    INDEX_ACRYLIC);
  fprintf(f, "%lf   !ascint        1/e attenuation length for scintillator\n",
    SCINT_ATTEN_LEN);
  fprintf(f, "%lf   !agfilm        1/e attenuation length (cm) for guide film\n",
    GUIDE_FILM_ATTEN);
  fprintf(f, "%lf   !agbuffer      1/e attenuation length (cm) for guide fluid\n",
    ACRYLIC_ATTEN_LEN);
  fprintf(f, "%lf   !afilm1        1: 1/e attenuation for acrylic plate\n",
    GAP_1_ATTEN_LEN);
  fprintf(f, "%d            !nfilms1       1: nfilms\n", N_FILMS_1);   
  fprintf(f, "%lf     !rgap1         1: gap index air\n", INDEX_GAP_1);
  fprintf(f, "%lf     !gapwidth      1: region without reflection on 6 sides of cell\n", GAP_WIDTH_1);
  fprintf(f, "%lf   !afilm2        2: 1/e attenuation for film\n", GAP_2_ATTEN_LEN);        
  fprintf(f, "%d            !nfilms2       2: nfilms\n", N_FILMS_2);           
  fprintf(f, "%lf     !rgap2         2: gap index perfluorooctane\n", INDEX_GAP_2);
  fprintf(f, "%lf     !gapwidth      2: region without reflection on 6 sides of cell\n", GAP_WIDTH_2);
  fprintf(f, "%lf     !afilm3        3: 1/e attenuation for film\n", GAP_3_ATTEN_LEN);       
  fprintf(f, "%d            !nfilms3       3: nfilms\n", N_FILMS_3);           
  fprintf(f, "%lf     !rgap3         3: gap index teflon\n", INDEX_GAP_3);                    
  fprintf(f, "%lf     !gapwidth      3: region without reflection on 6 sides of cell\n", GAP_WIDTH_3);
}

//*************************************************************************************************

int checkReconData(FILE * f)
{
  //Declear needed variables
  char line[stringLen_g];
  char description[stringLen_g];
  double parameter;
  int successFlag=0;

  //Read and check the preamble
  for (int i=0; i<30; i++)
  {
    initString(line, stringLen_g);
    fgets(line, stringLen_g, f);
    sscanf(line, "%lf %s", &parameter, description);
    
    if      (i==0)  { if (parameter==PMT_REFLECT)           successFlag++; }
    else if (i==1)  { if (int(parameter)==TIME_BIN_SIZE)    successFlag++; }
    else if (i==2)  { if (int(parameter)==N_TIME_BIN_PRINT) successFlag++; }
    else if (i==3)  { if (parameter==CELL_DIM)              successFlag++; }
    else if (i==4)  { if (parameter==PMT_RADIUS)            successFlag++; }
    else if (i==5)  { if (parameter==GAP_FILM_THICKNESS)    successFlag++; }
    else if (i==6)  { if (parameter==GUIDE_LENGTH)          successFlag++; }
    else if (i==7)  { if (parameter==SQUARE_LENGTH)         successFlag++; }
    else if (i==8)  { if (parameter==REFLECT_COEFF)         successFlag++; }
    else if (i==9)  { if (parameter==GUIDE_FILM_THICKNESS)  successFlag++; }
    else if (i==10) { if (parameter==PC_LIGHT_YIELD)        successFlag++; }
    else if (i==11) { if (parameter==QUANTUM_EFF)           successFlag++; }
    else if (i==12) { if (parameter==RELATIVE_YIELD)        successFlag++; }
    else if (i==13) { if (parameter==INDEX_SCINT)           successFlag++; }
    else if (i==14) { if (parameter==INDEX_ACRYLIC)         successFlag++; }
    else if (i==15) { if (parameter==SCINT_ATTEN_LEN)       successFlag++; }
    else if (i==16) { if (parameter==GUIDE_FILM_ATTEN)      successFlag++; }
    else if (i==17) { if (parameter==ACRYLIC_ATTEN_LEN)     successFlag++; }
    else if (i==18) { if (parameter==GAP_1_ATTEN_LEN)       successFlag++; }
    else if (i==19) { if (int(parameter)==N_FILMS_1)        successFlag++; }
    else if (i==20) { if (parameter==INDEX_GAP_1)           successFlag++; }
    else if (i==21) { if (parameter==GAP_WIDTH_1)           successFlag++; }
    else if (i==22) { if (parameter==GAP_2_ATTEN_LEN)       successFlag++; }
    else if (i==23) { if (int(parameter)==N_FILMS_2)        successFlag++; }
    else if (i==24) { if (parameter==INDEX_GAP_2)           successFlag++; }
    else if (i==25) { if (parameter==GAP_WIDTH_2)           successFlag++; }
    else if (i==26) { if (parameter==GAP_3_ATTEN_LEN)       successFlag++; }
    else if (i==27) { if (int(parameter)==N_FILMS_3)        successFlag++; }
    else if (i==28) { if (parameter==INDEX_GAP_3)           successFlag++; }
    else            { if (parameter==GAP_WIDTH_3)           successFlag++; }
  } 
  
  //Return success or failure
  if (successFlag==30) { return 1; }
  else { return 0; }
}

//*************************************************************************************************
