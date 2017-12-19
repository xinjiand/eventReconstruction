#include "reconLib.hh"

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&

void getEnergyScaleFactors(double * slope, double * intercept)
{
  //This function gives the slope and the intercept from a linear fit to convert photoelectrons to 
  //energy. If the macro PE_TO_ENERGY_CONVERSION_METHOD is set to 0 then only the slope is used.
  //
  //Note that the sum of the slope and/or the intercept basically sets the lowest threshold for the
  //detector since the numbers returned is the energy needed on average to produce one photoelectron
  //
  //These factors were determined by running single cell events with energies of 0.4, 0.5, 1.0, 2.0
  //3.0, and 4.0 MeV and histograming the total reconstructed charge. These centroid for the
  //histograms were determined and a linear fit done to the ordered pairs (PE, Energy). After this
  //positrons with energies of 1, 2, and 3 MeV were run to check that code is correctly
  //reconstructing the physics

  //Declare needed variables
  double slopeOfFit, interceptOfFit;  //Units of MeV/PE and MeV respectivily
  char detectorIdString[100];         //
  int detectorID;                     //
  int printWarningFlag=0;             //
  
  //Get the detector ID
  for (int i=0; i<100; i++) { detectorIdString[i]='\0'; }
  sprintf(detectorIdString, "%d%d%d%d", NCELL, LATTICE_FLAG, MIRROR_FLAG, GUIDE_FLAG);
  detectorID = atoi(detectorIdString);
  
  //Get the correction factor
  slopeOfFit = 1.0;
  interceptOfFit = 0.0;
  
  //TO DO: Switch these to switch statements
  
  if (detectorID==15111)
  {
    slopeOfFit = 0.00311;
    interceptOfFit = 0.0;
    if (SCINT_ATTEN_LEN!=400.0) printWarningFlag=1;
  }
  else if (detectorID==15101)
  {
    slopeOfFit = 0.00311;
    interceptOfFit = 0.0;
    if (SCINT_ATTEN_LEN!=93.0) printWarningFlag=1;  
  }
  else if (detectorID==5111)
  {
    slopeOfFit = 0.00359;
    interceptOfFit = 0.00377;
    if (SCINT_ATTEN_LEN!=93.0) printWarningFlag=1;
  }
  else if (detectorID==5110)
  {
    slopeOfFit = 0.00525;
    interceptOfFit = -0.01595;
    if (SCINT_ATTEN_LEN!=93.0) printWarningFlag=1;
  }
  else if (detectorID==5101)
  {
    slopeOfFit = 0.00455;
    interceptOfFit = 0.00548;
    if (SCINT_ATTEN_LEN!=93.0) printWarningFlag=1;
  }
  else if (detectorID==5100)
  {
    slopeOfFit = 0.00845;
    interceptOfFit = 0.01305;
    if (SCINT_ATTEN_LEN!=93.0) printWarningFlag=1;
  }
  else if (detectorID==9301)
  {
    slopeOfFit = 1.0; //0.0064410551
    interceptOfFit = 0.0;  
  }
  else
  {
    printf("WARNING! The photoelectron to MeV conversion factors have not been determined!\n");
  }

  if (PE_TO_ENERGY_CONVERSION_METHOD==0) interceptOfFit = 0.0;

  //Print needed warnings
  if (printWarningFlag==1)
  {
    printf("WARNING! The energy-PE conversion fit was determined with a different attenuation length.\n");
  }

  //Set values  
  *slope = slopeOfFit;
  *intercept = interceptOfFit;
  
}

//**************************************************************************************************

double getEnergy(double pe, double slope, double intercept)
{
  double energy = slope*pe + intercept;
  return energy;
}

//**************************************************************************************************

void sortDeposits(deposit_t * list, int len, int sortFlag)
{
  //Declare needed variables
  deposit_t temp;
  int swapCnt=1;
  
  //Bubble sort the array
  while (swapCnt)
  {
    swapCnt=0;
    for (int i=0; i<len-1; i++)
    {
      if (sortFlag==0)  //sort by energy
      {
        if (list[i].energy<list[i+1].energy)
        {
          temp.x = list[i].x; temp.y = list[i].y; temp.z = list[i].z;
          temp.t = list[i].t; temp.energy = list[i].energy; temp.id = list[i].id;
          
          list[i].x = list[i+1].x; list[i].y = list[i+1].y; list[i].z = list[i+1].z;
          list[i].t = list[i+1].t; list[i].energy = list[i+1].energy; list[i].id = list[i+1].id;
          
          list[i+1].x = temp.x; list[i+1].y = temp.y; list[i+1].z = temp.z;
          list[i+1].t = temp.t; list[i+1].energy = temp.energy; list[i+1].id = temp.id;
          
          swapCnt++;
        } //End ordering check
      }
      else if (sortFlag==1) //sort by time
      {
        if (list[i].t>list[i+1].t)
        {
          temp.x = list[i].x; temp.y = list[i].y; temp.z = list[i].z;
          temp.t = list[i].t; temp.energy = list[i].energy; temp.id = list[i].id;
          
          list[i].x = list[i+1].x; list[i].y = list[i+1].y; list[i].z = list[i+1].z;
          list[i].t = list[i+1].t; list[i].energy = list[i+1].energy; list[i].id = list[i+1].id;
          
          list[i+1].x = temp.x; list[i+1].y = temp.y; list[i+1].z = temp.z;
          list[i+1].t = temp.t; list[i+1].energy = temp.energy; list[i+1].id = temp.id;
          
          swapCnt++;
        } //End ordering check      
      } //End determine sort time
    } //End loop over elments in the array
  } //End while loop
} //End function

//**************************************************************************************************

void sortComponents(component_t * c, int len, int sortFlag)
{
  //Declare needed variables
  component_t temp;
  int swapCnt=1;
  
  //Bubble sort the array
  while (swapCnt)
  {
    swapCnt=0;
    for (int i=0; i<len-1; i++)
    {
      if (sortFlag==0)  //sort by energy
      {
        if (c[i].e<c[i+1].e)
        {
          temp.nx = c[i].nx; temp.ny = c[i].ny; temp.nz = c[i].nz;
          temp.t = c[i].t; temp.e = c[i].e;
          
          c[i].nx = c[i+1].nx; c[i].ny = c[i+1].ny; c[i].nz = c[i+1].nz;
          c[i].t = c[i+1].t; c[i].e = c[i+1].e;
          
          c[i+1].nx = temp.nx; c[i+1].ny = temp.ny; c[i+1].nz = temp.nz;
          c[i+1].t = temp.t; c[i+1].e = temp.e;
          
          swapCnt++;
        } //End ordering check
      }
      else if (sortFlag==1) //sort by time
      {
        if (c[i].t>c[i+1].t)
        {
          temp.nx = c[i].nx; temp.ny = c[i].ny; temp.nz = c[i].nz;
          temp.t = c[i].t; temp.e = c[i].e;
          
          c[i].nx = c[i+1].nx; c[i].ny = c[i+1].ny; c[i].nz = c[i+1].nz;
          c[i].t = c[i+1].t; c[i].e = c[i+1].e;
          
          c[i+1].nx = temp.nx; c[i+1].ny = temp.ny; c[i+1].nz = temp.nz;
          c[i+1].t = temp.t; c[i+1].e = temp.e;
          
          swapCnt++;
        } //End ordering check      
      } //End determine sort time
    } //End loop over elments in the array
  } //End while loop

} //End function

//**************************************************************************************************

void clearDetector(detector_t * det)
{
  for (int i=0; i<NCELL+2; i++) {
    for (int j=0; j<NCELL+2; j++) {
      for (int k=0; k<NCELL+2; k++) {
        det->component[i][j][k] = 0.0;
  } } }
}

//**************************************************************************************************

void clearComponents(component_t * c, int len)
{
  for (int i=0; i<len; i++)
  { c[i].nx = 0;  c[i].ny = 0;  c[i].nz = 0;  c[i].e = 0.0; c[i].t = 0.0; }
}

//**************************************************************************************************

void copyDetector(detector_t * in, detector_t * out)
{
  clearDetector(out);
  for (int i=0; i<NCELL+2; i++) {
    for (int j=0; j<NCELL+2; j++) {
      for (int k=0; k<NCELL+2; k++) {
        out->component[i][j][k] = in->component[i][j][k];
  } } }  
}

//**************************************************************************************************

double getChiSqrCharge(double * primary, double * side, component_t * pmtList, component_t * cellList,
  int pmtCnt, int cellCnt)
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Delcare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  detector_t * real, * recon;     //
  double chiSqr=0.0, sideCharge;  //
  double realQ, reconQ;           //
  int stop;                       //

  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&

  //Ensure there are only positve energy deposits
  absComponentEnergy(cellList, cellCnt);
  
  //Set the observed PMTs
  real = new detector_t;
  clearDetector(real); 
  for (int i=0; i<pmtCnt; i++)
  {
    real->component[pmtList[i].nx][pmtList[i].ny][pmtList[i].nz] = pmtList[i].e;
  }
  
  //Set the reconstructed hit cells
  recon = new detector_t;
  clearDetector(recon);
  for (int i=0; i<cellCnt; i++)
  {
    recon->component[cellList[i].nx][cellList[i].ny][cellList[i].nz] = cellList[i].e;
  }  

  //Determine in the charges in the PMT's from the hit cells
  for (int i=1; i<=NCELL; i++) {
    for (int j=1; j<=NCELL; j++) {
      for (int k=1; k<=NCELL; k++) {
        if (recon->component[i][j][k]>0.0)
        {
          //The primary PMTs on the x,y,z=0 sides
          recon->component[0][j][k] += recon->component[i][j][k]/primary[i-1];          
          recon->component[i][0][k] += recon->component[i][j][k]/primary[j-1];          
          recon->component[i][j][0] += recon->component[i][j][k]/primary[k-1];
          
          //The primary PMTs on the x,y,z=NCELL+1 sides
          if (MIRROR_FLAG==0)
          {
            recon->component[NCELL+1][j][k] += recon->component[i][j][k]/primary[NCELL-i];
            recon->component[i][NCELL+1][k] += recon->component[i][j][k]/primary[NCELL-j];
            recon->component[i][j][NCELL+1] += recon->component[i][j][k]/primary[NCELL-k];
          }
          
          //The x=0 side PMTs
          sideCharge = (recon->component[i][j][k]/primary[i-1])*side[i-1];
          if (j>1)     recon->component[0][j-1][k] += sideCharge;
          if (j<NCELL) recon->component[0][j+1][k] += sideCharge;
          
          if (k>1)     recon->component[0][j][k-1] += sideCharge;
          if (k<NCELL) recon->component[0][j][k+1] += sideCharge;
          
          //The y=0 side PMTs
          sideCharge = (recon->component[i][j][k]/primary[j-1])*side[j-1];
          if (i>1)     recon->component[i-1][0][k] += sideCharge;
          if (i<NCELL) recon->component[i+1][0][k] += sideCharge;
          
          if (k>1)     recon->component[i][0][k-1] += sideCharge;
          if (k<NCELL) recon->component[i][0][k+1] += sideCharge;          

          //The z=0 side PMTs
          sideCharge = (recon->component[i][j][k]/primary[k-1])*side[k-1];
          if (i>1)     recon->component[i-1][j][0] += sideCharge;
          if (i<NCELL) recon->component[i+1][j][0] += sideCharge;
          
          if (j>1)     recon->component[i][j-1][0] += sideCharge;
          if (j<NCELL) recon->component[i][j+1][0] += sideCharge;
          
          if (MIRROR_FLAG==0)
          {
            //The x=0 side PMTs
            sideCharge = (recon->component[i][j][k]/primary[NCELL-i])*side[NCELL-i];
            if (j>1)     recon->component[NCELL+1][j-1][k] += sideCharge;
            if (j<NCELL) recon->component[NCELL+1][j+1][k] += sideCharge;
            
            if (k>1)     recon->component[NCELL+1][j][k-1] += sideCharge;
            if (k<NCELL) recon->component[NCELL+1][j][k+1] += sideCharge;
            
            //The y=NCELL+1 side PMTs
            sideCharge = (recon->component[i][j][k]/primary[NCELL-j])*side[NCELL-j];
            if (i>1)     recon->component[i-1][NCELL+1][k] += sideCharge;
            if (i<NCELL) recon->component[i+1][NCELL+1][k] += sideCharge;
            
            if (k>1)     recon->component[i][NCELL+1][k-1] += sideCharge;
            if (k<NCELL) recon->component[i][NCELL+1][k+1] += sideCharge;          

            //The z=NCELL+1 side PMTs
            sideCharge = (recon->component[i][j][k]/primary[NCELL-k])*side[NCELL-k];
            if (i>1)     recon->component[i-1][j][NCELL+1] += sideCharge;
            if (i<NCELL) recon->component[i+1][j][NCELL+1] += sideCharge;
            
            if (j>1)     recon->component[i][j-1][NCELL+1] += sideCharge;
            if (j<NCELL) recon->component[i][j+1][NCELL+1] += sideCharge;           
          }          
        } //End check if there is energy in the cell
  } } } //End loop over cells
  
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Determine the Chi-Square &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  stop = (1-MIRROR_FLAG)*(NCELL+1);
  for (int i=0; i<=stop; i+=NCELL+1) {
    for (int j=1; j<=NCELL; j++) {
      for (int k=1; k<=NCELL; k++) {
        realQ = real->component[i][j][k];  //x sides
        reconQ = recon->component[i][j][k];
        
        if (realQ==0.0)
        {
          chiSqr += (reconQ-realQ)*(reconQ-realQ);
        }
        else if (realQ>0.0)
        {
          chiSqr += (((reconQ-realQ)*(reconQ-realQ))/realQ);        
        } 
        
        realQ = real->component[j][i][k];  //y sides
        reconQ = recon->component[j][i][k];

        if (realQ==0.0)
        {
          chiSqr += (reconQ-realQ)*(reconQ-realQ);
        }
        else if (realQ>0.0)
        {
          chiSqr += (((reconQ-realQ)*(reconQ-realQ))/realQ);        
        }
      
        realQ = real->component[j][k][i];  //z sides
        reconQ = recon->component[j][k][i];

        if (realQ==0.0)
        {
          chiSqr += (reconQ-realQ)*(reconQ-realQ);
        }
        else if (realQ>0.0)
        {
          chiSqr += (((reconQ-realQ)*(reconQ-realQ))/realQ);        
        }
  } } }
  
  //Clean up memory and end
  delete real;
  delete recon;
  
  return chiSqr;
}

//**************************************************************************************************

double getChiSqrWaveform(double ** primary, double ** side, component_t * pmtList, double ** hitPMT,
                          component_t * cellList, int pmtCnt, int cellCnt)
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Delcare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  component_t * fullPMT_List;
  component_t tempCell;
  double ** realPMT;  
  double ** reconPMT;
  double chiSqr=0.0, realQ, reconQ;
  int nPMT=3*NCELL*NCELL*(2-MIRROR_FLAG);
  int pmtIndex, stop;
  
  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&

  //Ensure there are only positve energy deposits
  absComponentEnergy(cellList, cellCnt);
  
  //Ensure all times are positive
  if (TIME_RECON) resetTimeZero(cellList, cellCnt);  
  
  //Allocate the PMTs
  realPMT = new double* [nPMT];
  reconPMT = new double* [nPMT];
  for (int i=0; i<nPMT; i++)
  {
    realPMT[i] = new double [N_TIME_BIN_PRINT];
    reconPMT[i] = new double [N_TIME_BIN_PRINT];
  }

  //Set the full PMT list
  fullPMT_List = new component_t [nPMT];
  stop = (1-MIRROR_FLAG)*(NCELL+1);
  pmtIndex=0;
  for (int i=0; i<=stop; i+=NCELL+1) {
    for (int j=1; j<=NCELL; j++) {
      for (int k=1; k<=NCELL; k++) {
        fullPMT_List[pmtIndex].nx=i; fullPMT_List[pmtIndex].ny=j; fullPMT_List[pmtIndex].nz=k;
        fullPMT_List[pmtIndex].e=0.0; fullPMT_List[pmtIndex].t=0.0;  
        
        fullPMT_List[pmtIndex+1].nx=j; fullPMT_List[pmtIndex+1].ny=i; fullPMT_List[pmtIndex+1].nz=k;        
        fullPMT_List[pmtIndex+1].e=0.0; fullPMT_List[pmtIndex+1].t=0.0;  
        
        fullPMT_List[pmtIndex+2].nx=j; fullPMT_List[pmtIndex+2].ny=k; fullPMT_List[pmtIndex+2].nz=i;
        fullPMT_List[pmtIndex+2].e=0.0; fullPMT_List[pmtIndex+2].t=0.0;  
                        
        pmtIndex+=3;
  } } }
  
  //Set the real PMTs
  for (int i=0; i<nPMT; i++)       //First zero out the local waveforms
  {
    for (int j=0; j<N_TIME_BIN_PRINT; j++)
    {
      realPMT[i][j]=0.0;
    }
  }

  for (int i=0; i<pmtCnt; i++)     //Next loop over the hit pmts and assign
  {
    for (int j=0; j<nPMT; j++)
    {
      if (pmtList[i].nx==fullPMT_List[j].nx && pmtList[i].ny==fullPMT_List[j].ny &&
        pmtList[i].nz==fullPMT_List[j].nz)
      {
        //Assign the hit pmt to the full pmt list
        fullPMT_List[j].e = pmtList[i].e;
        fullPMT_List[j].t = pmtList[i].t;
        
        //Copy the waveform
        for (int k=0; k<N_TIME_BIN_PRINT; k++)
        {
          realPMT[j][k] = hitPMT[i][k];        
        }
      }
    } //End loop over the full pmt list
  } //End loop over hit pmts
 
  //Set the reconstructed PMTs 
  for (int i=0; i<nPMT; i++)       //First zero out the local waveforms
  {
    for (int j=0; j<N_TIME_BIN_PRINT; j++)
    {
      reconPMT[i][j]=0.0;
    }
  }

  for (int i=0; i<cellCnt; i++)
  {  
    //Get the primary pmt waveforms
    getWaveform(primary, cellList[i], 0, fullPMT_List, reconPMT); //x=0
    getWaveform(primary, cellList[i], 2, fullPMT_List, reconPMT); //y=0
    getWaveform(primary, cellList[i], 4, fullPMT_List, reconPMT); //z=0
    
    if (MIRROR_FLAG==0)
    {
      getWaveform(primary, cellList[i], 1, fullPMT_List, reconPMT); //x=NCELL+1
      getWaveform(primary, cellList[i], 3, fullPMT_List, reconPMT); //y=NCELL+1
      getWaveform(primary, cellList[i], 5, fullPMT_List, reconPMT); //z=NCELL+1
    } //End add pmts for full instrumentation
    
    //Set the side PMTs
    
    //x=0 and x=NCELL+1 side
    tempCell.nx = cellList[i].nx;
    tempCell.ny = cellList[i].ny;
    tempCell.nz = cellList[i].nz;
    tempCell.e = cellList[i].e;
    tempCell.t = cellList[i].t;

    if (cellList[i].ny>1)
    {
      tempCell.ny--;
      getWaveform(side, tempCell, 0, fullPMT_List, reconPMT);
      if (MIRROR_FLAG==0) getWaveform(side, tempCell, 1, fullPMT_List, reconPMT);
      tempCell.ny++;      
    }
    
    if (cellList[i].ny<NCELL)
    {
      tempCell.ny++;
      getWaveform(side, tempCell, 0, fullPMT_List, reconPMT);          
      if (MIRROR_FLAG==0) getWaveform(side, tempCell, 1, fullPMT_List, reconPMT);
      tempCell.ny--;    
    }
    
    if (cellList[i].nz>1)
    {
      tempCell.nz--;
      getWaveform(side, tempCell, 0, fullPMT_List, reconPMT);          
      if (MIRROR_FLAG==0) getWaveform(side, tempCell, 1, fullPMT_List, reconPMT);
      tempCell.nz++;   
    }
    
    if (cellList[i].nz<NCELL)
    {
      tempCell.nz++;
      getWaveform(side, tempCell, 0, fullPMT_List, reconPMT);          
      if (MIRROR_FLAG==0) getWaveform(side, tempCell, 1, fullPMT_List, reconPMT);
      tempCell.nz--;     
    }
    
    //y=0 and y=NCELL+1 side
    if (cellList[i].nx>1)
    {
      tempCell.nx--;
      getWaveform(side, tempCell, 2, fullPMT_List, reconPMT);
      if (MIRROR_FLAG==0) getWaveform(side, tempCell, 3, fullPMT_List, reconPMT);
      tempCell.nx++;      
    }
    
    if (cellList[i].nx<NCELL)
    {
      tempCell.nx++;
      getWaveform(side, tempCell, 2, fullPMT_List, reconPMT);          
      if (MIRROR_FLAG==0) getWaveform(side, tempCell, 3, fullPMT_List, reconPMT);
      tempCell.nx--;    
    }
    
    if (cellList[i].nz>1)
    {
      tempCell.nz--;
      getWaveform(side, tempCell, 2, fullPMT_List, reconPMT);          
      if (MIRROR_FLAG==0) getWaveform(side, tempCell, 3, fullPMT_List, reconPMT);
      tempCell.nz++;   
    }
    
    if (cellList[i].nz<NCELL)
    {
      tempCell.nz++;
      getWaveform(side, tempCell, 2, fullPMT_List, reconPMT);          
      if (MIRROR_FLAG==0) getWaveform(side, tempCell, 3, fullPMT_List, reconPMT);
      tempCell.nz--;     
    }
    
    //z=0 and z=NCELL+1 side
    if (cellList[i].nx>1)
    {
      tempCell.nx--;
      getWaveform(side, tempCell, 4, fullPMT_List, reconPMT);
      if (MIRROR_FLAG==0) getWaveform(side, tempCell, 5, fullPMT_List, reconPMT);
      tempCell.nx++;      
    }
    
    if (cellList[i].nx<NCELL)
    {
      tempCell.nx++;
      getWaveform(side, tempCell, 4, fullPMT_List, reconPMT);          
      if (MIRROR_FLAG==0) getWaveform(side, tempCell, 5, fullPMT_List, reconPMT);
      tempCell.nx--;    
    }
    
    if (cellList[i].ny>1)
    {
      tempCell.ny--;
      getWaveform(side, tempCell, 4, fullPMT_List, reconPMT);          
      if (MIRROR_FLAG==0) getWaveform(side, tempCell, 5, fullPMT_List, reconPMT);
      tempCell.ny++;   
    }
    
    if (cellList[i].ny<NCELL)
    {
      tempCell.ny++;
      getWaveform(side, tempCell, 4, fullPMT_List, reconPMT);          
      if (MIRROR_FLAG==0) getWaveform(side, tempCell, 5, fullPMT_List, reconPMT);
      tempCell.ny--;     
    }     
  } //End loop over the hit cells

  
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Determine the Chi-Square &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  pmtIndex=0;
  for (int i=0; i<=stop; i+=NCELL+1) {
    for (int j=1; j<=NCELL; j++) {
      for (int k=1; k<=NCELL; k++) {
        //The x-sides
        for (int ii=0; ii<N_TIME_BIN_PRINT; ii++)
        {
          realQ = realPMT[pmtIndex][ii];
          reconQ = reconPMT[pmtIndex][ii];

          if (realQ==0.0)
          {
            chiSqr += (reconQ-realQ)*(reconQ-realQ);
          }
          else
          {
            chiSqr += (((reconQ-realQ)*(reconQ-realQ))/realQ);        
          }        
        }

        //The y-sides
        for (int ii=0; ii<N_TIME_BIN_PRINT; ii++)
        {
          realQ = realPMT[pmtIndex+1][ii];
          reconQ = reconPMT[pmtIndex+1][ii];

          if (realQ==0.0)
          {
            chiSqr += (reconQ-realQ)*(reconQ-realQ);
          }
          else
          {
            chiSqr += (((reconQ-realQ)*(reconQ-realQ))/realQ);        
          }            
        }

        //The z-sides
        for (int ii=0; ii<N_TIME_BIN_PRINT; ii++)
        {
          realQ = realPMT[pmtIndex+2][ii];
          reconQ = reconPMT[pmtIndex+2][ii];

          if (realQ==0.0)
          {
            chiSqr += (reconQ-realQ)*(reconQ-realQ);
          }
          else
          {
            chiSqr += (((reconQ-realQ)*(reconQ-realQ))/realQ);        
          }        
        }

        pmtIndex+=3;
  } } }
 
  
  //Testing statements
  if (TEST_FLAG_RECON==2)
  {
    pmtIndex = getHighEnergyCellsPmtIndex(fullPMT_List, cellList, nPMT, cellCnt);
    
    printf("Reconstructed PMT Waveform:\n");
    for (int bin=0; bin<N_TIME_BIN_PRINT; bin++)
    {
      printf("%d %lf\n", bin, reconPMT[pmtIndex][bin]);
    }
    printf("\n");
  }


  //Clean up memory
  for (int i=0; i<nPMT; i++)
  {
    delete [] realPMT[i];
    delete [] reconPMT[i];
  }
  delete [] realPMT;
  delete [] reconPMT;
  delete [] fullPMT_List;
   
  return chiSqr;
}

//**************************************************************************************************

double getChiSqrTeflon(double * primary, double ** side, component_t * pmtList, component_t * cellList,
  int pmtCnt, int cellCnt)
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Delcare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  detector_t * real, * recon;     //
  double chiSqr=0.0, sideCharge;  //
  double realQ, reconQ;           //
  int nSidePMTs, stop;            //


  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&
  
  //Get the number of side PMTs
  nSidePMTs = NCELL - 1;

  //Ensure there are only positve energy deposits
  absComponentEnergy(cellList, cellCnt);
  
  //Set the observed PMTs
  real = new detector_t;
  clearDetector(real); 
  for (int i=0; i<pmtCnt; i++)
  {
    real->component[pmtList[i].nx][pmtList[i].ny][pmtList[i].nz] = pmtList[i].e;
  }
  
  //Set the reconstructed hit cells
  recon = new detector_t;
  clearDetector(recon);
  for (int i=0; i<cellCnt; i++)
  {
    recon->component[cellList[i].nx][cellList[i].ny][cellList[i].nz] = cellList[i].e;
  }  

  //Determine in the charges in the PMT's from the hit cells
  for (int i=1; i<=NCELL; i++) {
    for (int j=1; j<=NCELL; j++) {
      for (int k=1; k<=NCELL; k++) {
        if (recon->component[i][j][k]>0.0)
        {
          //The primary PMTs on the x,y,z=0 sides
          recon->component[0][j][k] += recon->component[i][j][k]/primary[i-1];          
          recon->component[i][0][k] += recon->component[i][j][k]/primary[j-1];          
          recon->component[i][j][0] += recon->component[i][j][k]/primary[k-1];
          
          //The primary PMTs on the x,y,z=NCELL+1 sides
          if (MIRROR_FLAG==0)
          {
            recon->component[NCELL+1][j][k] += recon->component[i][j][k]/primary[NCELL-i];
            recon->component[i][NCELL+1][k] += recon->component[i][j][k]/primary[NCELL-j];
            recon->component[i][j][NCELL+1] += recon->component[i][j][k]/primary[NCELL-k];
          }

          for (int ii=0; ii<nSidePMTs; ii++)
          {          
            //The x=0 side PMTs
            sideCharge = (recon->component[i][j][k]/primary[i-1])*side[i-1][ii];
            if (j-ii>1)     recon->component[0][j-ii-1][k] += sideCharge;
            if (j+ii<NCELL) recon->component[0][j+ii+1][k] += sideCharge;
            
            if (k-ii>1)     recon->component[0][j][k-ii-1] += sideCharge;
            if (k+ii<NCELL) recon->component[0][j][k+ii+1] += sideCharge;
                        
            //The y=0 side PMTs
            sideCharge = (recon->component[i][j][k]/primary[j-1])*side[j-1][ii];
            if (i-ii>1)     recon->component[i-ii-1][0][k] += sideCharge;
            if (i+ii<NCELL) recon->component[i+ii+1][0][k] += sideCharge;
            
            if (k-ii>1)     recon->component[i][0][k-ii-1] += sideCharge;
            if (k+ii<NCELL) recon->component[i][0][k+ii+1] += sideCharge;          

            //The z=0 side PMTs
            sideCharge = (recon->component[i][j][k]/primary[k-1])*side[k-1][ii];
            if (i-ii>1)     recon->component[i-ii-1][j][0] += sideCharge;
            if (i+ii<NCELL) recon->component[i+ii+1][j][0] += sideCharge;
            
            if (j-ii>1)     recon->component[i][j-ii-1][0] += sideCharge;
            if (j+ii<NCELL) recon->component[i][j+ii+1][0] += sideCharge;
            
            if (MIRROR_FLAG==0)
            {
              //The x=NCELL+1 side PMTs
              sideCharge = (recon->component[i][j][k]/primary[NCELL-i])*side[NCELL-i][ii];
              if (j-ii>1)     recon->component[NCELL+1][j-ii-1][k] += sideCharge;
              if (j+ii<NCELL) recon->component[NCELL+1][j+ii+1][k] += sideCharge;
              
              if (k-ii>1)     recon->component[NCELL+1][j][k-ii-1] += sideCharge;
              if (k+ii<NCELL) recon->component[NCELL+1][j][k+ii+1] += sideCharge;
              
              //The y=NCELL+1 side PMTs
              sideCharge = (recon->component[i][j][k]/primary[NCELL-j])*side[NCELL-j][ii];
              if (i-ii>1)     recon->component[i-ii-1][NCELL+1][k] += sideCharge;
              if (i+ii<NCELL) recon->component[i+ii+1][NCELL+1][k] += sideCharge;
              
              if (k-ii>1)     recon->component[i][NCELL+1][k-ii-1] += sideCharge;
              if (k+ii<NCELL) recon->component[i][NCELL+1][k+ii+1] += sideCharge;          

              //The z=NCELL+1 side PMTs
              sideCharge = (recon->component[i][j][k]/primary[NCELL-k])*side[NCELL-k][ii];
              if (i-ii>1)     recon->component[i-ii-1][j][NCELL+1] += sideCharge;
              if (i+ii<NCELL) recon->component[i+ii+1][j][NCELL+1] += sideCharge;
              
              if (j-ii>1)     recon->component[i][j-ii-1][NCELL+1] += sideCharge;
              if (j+ii<NCELL) recon->component[i][j+ii+1][NCELL+1] += sideCharge;           
            } //End check if using full instrumentation
          } //End loop over side PMTs       
        } //End check if there is energy in the cell
  } } } //End loop over cells

/*
  printf("Real\n");
  printDetectorPmts(real, 0);
  printf("Recon\n");
  printDetectorPmts(recon, 0);
*/  
  
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Determine the Chi-Square &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  stop = (1-MIRROR_FLAG)*(NCELL+1);
  for (int i=0; i<=stop; i+=NCELL+1) {
    for (int j=1; j<=NCELL; j++) {
      for (int k=1; k<=NCELL; k++) {
        realQ = real->component[i][j][k];  //x sides
        reconQ = recon->component[i][j][k];
        
        if (realQ==0.0)
        {
          chiSqr += (reconQ-realQ)*(reconQ-realQ);
        }
        else if (realQ>0.0)
        {
          chiSqr += (((reconQ-realQ)*(reconQ-realQ))/realQ);        
        } 
        
        realQ = real->component[j][i][k];  //y sides
        reconQ = recon->component[j][i][k];

        if (realQ==0.0)
        {
          chiSqr += (reconQ-realQ)*(reconQ-realQ);
        }
        else if (realQ>0.0)
        {
          chiSqr += (((reconQ-realQ)*(reconQ-realQ))/realQ);        
        }
      
        realQ = real->component[j][k][i];  //z sides
        reconQ = recon->component[j][k][i];

        if (realQ==0.0)
        {
          chiSqr += (reconQ-realQ)*(reconQ-realQ);
        }
        else if (realQ>0.0)
        {
          chiSqr += (((reconQ-realQ)*(reconQ-realQ))/realQ);        
        }
  } } }
  
  //Clean up memory and end
  delete real;
  delete recon;
  
  return chiSqr;
}

//**************************************************************************************************

void gradientSearchCharge(double * primary, double * side, component_t * pmtList,
  component_t * cellList, int pmtCnt, int cellCnt)
{
//**************************************************************************************************
// This function is adapted from the one in Bevington. It first determines the gradient and 
// normalizes it. To do this it takes an initial step along the gradient, if the chi square does not
// change then the step length is increased. Once the chi square has changed then the gradient is 
// normalized and multipled by a small length. See: http://en.wikipedia.org/wiki/Gradient_descent.
// Next an initial step is taken and the chi-squared checked to see that if it decreases. Once it
// does the chi-squared is reduced as much as possible by stepping along the negative gradient. 
//**************************************************************************************************

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Declare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  component_t * saveCells;            //
  double step, chiSqr1, norm, delta;  //
  double chiSqr2, chiSqr3;            //
  double * gradient;                  //
  int loopCnt, nNoChange;             //
  
  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&
  
  saveCells = new component_t [cellCnt];
  copyComponents(cellList, saveCells, cellCnt);

  step = STEP;

  gradient = new double [cellCnt];
  for (int i=0; i<cellCnt; i++) { gradient[i] = 0.0; }

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Find Gradient and Normalize &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  norm=0.0;
  loopCnt=0;
  chiSqr1 = getChiSqrCharge(primary, side, pmtList, cellList, pmtCnt, cellCnt); //Loop indep.
  while (norm==0.0 && loopCnt<=10) 
  {
    norm=0.0;

    for (int i=0; i<cellCnt; i++)
    {
      delta = step*sqrt(cellList[i].e/1000.0);
      cellList[i].e += delta;

      //Note the tricky minus sign here (ie this is the negative of how a derivative is defined)
      gradient[i] = chiSqr1 - getChiSqrCharge(primary, side, pmtList, cellList, pmtCnt, cellCnt);
      norm += gradient[i]*gradient[i];
      
      cellList[i].e -= delta;
    }
    norm = sqrt(norm);
    step = 2.0*step;    //If the chi-squared did not change then take a bigger step away from the
                        //current point.
    loopCnt++;
  }  
  if (loopCnt>=10)
  {
    if (TEST_FLAG_RECON==1) printf("Could not define the gradient.\n");
    delete [] saveCells;
    delete [] gradient;
    return;
  }        
  
  for (int i=0; i<cellCnt; i++)
  {
    gradient[i] = sqrt(cellList[i].e/1000.0)*gradient[i]/norm;
  }
  
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Minimize the Chi-Squared &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
  
  //Take the initial step
  chiSqr2 = INFINITY;
  nNoChange=0;
  while (chiSqr2>=chiSqr1 && nNoChange<=10)
  {   
    for (int i=0; i<cellCnt; i++)
    {
      cellList[i].e += gradient[i];
    } //Because the tricky minus sign above                         
                                                 
    chiSqr2 = getChiSqrCharge(primary, side, pmtList, cellList, pmtCnt, cellCnt);
    if (chiSqr2>=chiSqr1)
    {
      nNoChange++;
      
      for (int i=0; i<cellCnt; i++)
      {
        cellList[i].e -= gradient[i];
        gradient[i] = gradient[i]/2.0;  //Decrease the step
      }
    }
  }
  if (nNoChange>=10)
  {
    if (TEST_FLAG_RECON==1) printf("There was no decrease, reverting to passed reconstruction\n");
    copyComponents(saveCells, cellList, cellCnt);
    delete [] saveCells;
    delete [] gradient;
    return;
  }      

  //Step along the gradient until the chi-squared stays the same or goes up  
  while (true)
  {
    for (int i=0; i<cellCnt; i++)
    {
      cellList[i].e += gradient[i];
    }

    chiSqr3 = getChiSqrCharge(primary, side, pmtList, cellList, pmtCnt, cellCnt);
    if (chiSqr3<chiSqr2)
    {
      chiSqr1 = chiSqr2;
      chiSqr2 = chiSqr3;
    }
    else { break; }
  }    
  if (chiSqr2==chiSqr3) //If the chi-squareds are the same then return
  {
    delete [] saveCells;
    delete [] gradient;
    return;
  } 

  //If it went up, then interpolate to find minimum
  copyComponents(cellList, saveCells, cellCnt);
  delta = 1.0/(1.0 + (chiSqr1 - chiSqr2)/(chiSqr3 - chiSqr2)) + 0.5;
  for (int i=0; i<cellCnt; i++)
  {
    cellList[i].e -= delta*gradient[i];
  } 

  //Check that the inerpolation is actually a minimum
  chiSqr3 = getChiSqrCharge(primary, side, pmtList, cellList, pmtCnt, cellCnt);
  if (chiSqr2<chiSqr3)
  {
    copyComponents(saveCells, cellList, cellCnt);
    chiSqr2 = getChiSqrCharge(primary, side, pmtList, cellList, pmtCnt, cellCnt);
  }
  
  //&&&&&&&&&&&&&&&&
  //& End Function &
  //&&&&&&&&&&&&&&&&
  
  delete [] gradient;
  delete [] saveCells;

  return;
}

//**************************************************************************************************

void gradientSearchWaveform(double ** primary, double ** side, component_t * pmtList, double ** hitPMT,
  component_t * cellList, int pmtCnt, int cellCnt)
{
//**************************************************************************************************
// This function is adapted from the one in Bevington. It first determines the gradient and 
// normalizes it. To do this it takes an initial step along the gradient, if the chi square does not
// change then the step length is increased. Once the chi square has changed then the gradient is 
// normalized. See: http://en.wikipedia.org/wiki/Gradient_descent. Next an initial step is taken and
// the chi-squared checked to see that if it decreases. Once it does the chi-squared is reduced as
// much as possible by stepping along the negative gradient.
//**************************************************************************************************

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Declare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  component_t * saveCells;      //
  double step;                  //
  double chiSqr1, norm, delta;  //
  double chiSqr2, chiSqr3;      //
  double * gradient;            //
  int loopCnt, nNoChange;       //
  
  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&
  
  saveCells = new component_t [cellCnt];
  copyComponents(cellList, saveCells, cellCnt);

  if (TIME_RECON==0) step = 1.0;
  else step = 0.1;

  gradient = new double [cellCnt];
  for (int i=0; i<cellCnt; i++) { gradient[i] = 0.0; }

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Find Gradient and Normalize &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  norm=0.0;
  loopCnt=0;
  chiSqr1 = getChiSqrWaveform(primary, side, pmtList, hitPMT, cellList, pmtCnt, cellCnt); //Loop indep.
  while (norm==0.0 && loopCnt<=10) 
  {
    norm=0.0;
              
    for (int i=0; i<cellCnt; i++)
    {
      delta = step;
      cellList[i].t += delta;

      //Note the tricky minus sign here (ie this is the negative of how a derivative is defined)
      gradient[i] = chiSqr1 - getChiSqrWaveform(primary, side, pmtList, hitPMT, cellList, pmtCnt, 
        cellCnt);
      norm += gradient[i]*gradient[i];
 
      cellList[i].t -= delta;      

    }
    norm = sqrt(norm);
    step = 2.0*step;  //If the chi-squared did not change then take a bigger step away 
    
    loopCnt++;
  }  
  if (loopCnt>=10)
  {
    if (TEST_FLAG_RECON==1) printf("Could not define the gradient.\n");
    delete [] saveCells;
    delete [] gradient;
    return;
  }        
  
  if (TIME_RECON==0) step = 1.0;
  else step = 0.1;  
  
  for (int i=0; i<cellCnt; i++)
  {
    gradient[i] = step*(gradient[i]/norm);
  }

  
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Minimize the Chi-Squared &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
  
  //Take the initial step
  chiSqr2 = INFINITY;
  nNoChange=0;
  while (chiSqr2>=chiSqr1 && nNoChange<=10)
  {   
    for (int i=0; i<cellCnt; i++)
    {
      cellList[i].t += gradient[i];
    } //Because the tricky minus sign above                         
                                                 
    chiSqr2 = getChiSqrWaveform(primary, side, pmtList, hitPMT, cellList, pmtCnt, cellCnt);
    if (chiSqr2>=chiSqr1)
    {
      nNoChange++;
      
      for (int i=0; i<cellCnt; i++)
      {
        gradient[i] = gradient[i]/2.0;  //Decrease the step
      }      
    }
  }
  if (nNoChange>=10)
  {
    if (TEST_FLAG_RECON==1) printf("There was no decrease, reverting to passed reconstruction\n");
    copyComponents(saveCells, cellList, cellCnt);
    delete [] saveCells;
    delete [] gradient;
    return;
  }      

  //Step along the gradient until the chi-squared stays the same or goes up  
  while (true)
  {
    for (int i=0; i<cellCnt; i++)
    {
      cellList[i].t += gradient[i];
    }

    chiSqr3 = getChiSqrWaveform(primary, side, pmtList, hitPMT, cellList, pmtCnt, cellCnt);
    if (chiSqr3<chiSqr2)
    {
      chiSqr1 = chiSqr2;
      chiSqr2 = chiSqr3;
    }
    else { break; }
  }    
  if (chiSqr2==chiSqr3) //If the chi-squareds are the same then return
  {
    delete [] saveCells;
    delete [] gradient;
    return;
  } 

  //If it went up, then interpolate to find minimum
  copyComponents(cellList, saveCells, cellCnt);
  delta = 1.0/(1.0 + (chiSqr1 - chiSqr2)/(chiSqr3 - chiSqr2)) + 0.5;
  for (int i=0; i<cellCnt; i++)
  {
    cellList[i].t -= delta*gradient[i];
  } 

  //Check that the inerpolation is actually a minimum
  chiSqr3 = getChiSqrWaveform(primary, side, pmtList, hitPMT, cellList, pmtCnt, cellCnt);
  if (chiSqr2<chiSqr3)
  {
    copyComponents(saveCells, cellList, cellCnt);
    chiSqr2 = getChiSqrWaveform(primary, side, pmtList, hitPMT, cellList, pmtCnt, cellCnt);
  }
  
  //&&&&&&&&&&&&&&&&
  //& End Function &
  //&&&&&&&&&&&&&&&&
  
  delete [] gradient;
  delete [] saveCells;

  return;
}

//**************************************************************************************************

void gradientSearchTeflon(double * primary, double ** side, component_t * pmtList,
  component_t * cellList, int pmtCnt, int cellCnt)
{
//**************************************************************************************************
// This function is adapted from the one in Bevington. It first determines the gradient and 
// normalizes it. To do this it takes an initial step along the gradient, if the chi square does not
// change then the step length is increased. Once the chi square has changed then the gradient is 
// normalized and multipled by a small length. See: http://en.wikipedia.org/wiki/Gradient_descent.
// Next an initial step is taken and the chi-squared checked to see that if it decreases. Once it
// does the chi-squared is reduced as much as possible by stepping along the negative gradient. 
//**************************************************************************************************

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Declare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  component_t * saveCells;            //
  double step, chiSqr1, norm, delta;  //
  double chiSqr2, chiSqr3;            //
  double * gradient;                  //
  int loopCnt, nNoChange;             //
  
  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&
  
  saveCells = new component_t [cellCnt];
  copyComponents(cellList, saveCells, cellCnt);

  step = STEP;

  gradient = new double [cellCnt];
  for (int i=0; i<cellCnt; i++) { gradient[i] = 0.0; }

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Find Gradient and Normalize &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  norm=0.0;
  loopCnt=0;
  chiSqr1 = getChiSqrTeflon(primary, side, pmtList, cellList, pmtCnt, cellCnt); //Loop indep.
  while (norm==0.0 && loopCnt<=10) 
  {
    norm=0.0;

    for (int i=0; i<cellCnt; i++)
    {
      delta = step*sqrt(cellList[i].e/1000.0);
      cellList[i].e += delta;

      //Note the tricky minus sign here (ie this is the negative of how a derivative is defined)
      gradient[i] = chiSqr1 - getChiSqrTeflon(primary, side, pmtList, cellList, pmtCnt, cellCnt);
      norm += gradient[i]*gradient[i];
      
      cellList[i].e -= delta;
    }
    norm = sqrt(norm);
    step = 2.0*step;    //If the chi-squared did not change then take a bigger step away from the
                        //current point.
    loopCnt++;
  }  
  if (loopCnt>=10)
  {
    if (TEST_FLAG_RECON==1) printf("Could not define the gradient.\n");
    delete [] saveCells;
    delete [] gradient;
    return;
  }        
  
  for (int i=0; i<cellCnt; i++)
  {
    gradient[i] = sqrt(cellList[i].e/1000.0)*gradient[i]/norm;
  }
  
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Minimize the Chi-Squared &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
  
  //Take the initial step
  chiSqr2 = INFINITY;
  nNoChange=0;
  while (chiSqr2>=chiSqr1 && nNoChange<=10)
  {   
    for (int i=0; i<cellCnt; i++)
    {
      cellList[i].e += gradient[i];
    } //Because the tricky minus sign above                         
                                                 
    chiSqr2 = getChiSqrTeflon(primary, side, pmtList, cellList, pmtCnt, cellCnt);
    if (chiSqr2>=chiSqr1)
    {
      nNoChange++;
      
      for (int i=0; i<cellCnt; i++)
      {
        cellList[i].e -= gradient[i];
        gradient[i] = gradient[i]/2.0;  //Decrease the step
      }
    }
  }
  if (nNoChange>=10)
  {
    if (TEST_FLAG_RECON==1) printf("There was no decrease, reverting to passed reconstruction\n");
    copyComponents(saveCells, cellList, cellCnt);
    delete [] saveCells;
    delete [] gradient;
    return;
  }      

  //Step along the gradient until the chi-squared stays the same or goes up  
  while (true)
  {
    for (int i=0; i<cellCnt; i++)
    {
      cellList[i].e += gradient[i];
    }

    chiSqr3 = getChiSqrTeflon(primary, side, pmtList, cellList, pmtCnt, cellCnt);
    if (chiSqr3<chiSqr2)
    {
      chiSqr1 = chiSqr2;
      chiSqr2 = chiSqr3;
    }
    else { break; }
  }    
  if (chiSqr2==chiSqr3) //If the chi-squareds are the same then return
  {
    delete [] saveCells;
    delete [] gradient;
    return;
  } 

  //If it went up, then interpolate to find minimum
  copyComponents(cellList, saveCells, cellCnt);
  delta = 1.0/(1.0 + (chiSqr1 - chiSqr2)/(chiSqr3 - chiSqr2)) + 0.5;
  for (int i=0; i<cellCnt; i++)
  {
    cellList[i].e -= delta*gradient[i];
  } 

  //Check that the inerpolation is actually a minimum
  chiSqr3 = getChiSqrTeflon(primary, side, pmtList, cellList, pmtCnt, cellCnt);
  if (chiSqr2<chiSqr3)
  {
    copyComponents(saveCells, cellList, cellCnt);
    chiSqr2 = getChiSqrTeflon(primary, side, pmtList, cellList, pmtCnt, cellCnt);
  }
  
  //&&&&&&&&&&&&&&&&
  //& End Function &
  //&&&&&&&&&&&&&&&&
  
  delete [] gradient;
  delete [] saveCells;

  return;
}

//**************************************************************************************************

void absComponentEnergy(component_t * c, int len)
{
  for (int i=0; i<len; i++) { c[i].e = fabs(c[i].e); }
}

//**************************************************************************************************

void resetTimeZero(component_t * c, int len)
{
  double setTime0=INFINITY;
  for (int i=0; i<len; i++)
  {
    if (c[i].t<setTime0) setTime0 = c[i].t;
  }
  
  if (setTime0<0.0)
  {
    for (int i=0; i<len; i++) { c[i].t = c[i].t - setTime0; }
  }
}
 
//**************************************************************************************************

int getComponentIndex(component_t * c, int len, int nx, int ny, int nz)
{
  int index=-1;
  for (int i=0; i<len; i++)
  {
    if (c[i].nx==nx && c[i].ny==ny && c[i].nz==nz) { index = i; break; }
  }
  return index;
}

//************************************************************************************************** 

void getWaveform(double ** templateWaveform, component_t cell, int sideIndex, component_t * pmtList,
                  double ** waveform)
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Declare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  double slope, intercept;
  double local[N_TIME_BIN_PRINT];
  double sum=0.0;
  int nAway=0, pmtIndex=-1;  
  int nu=0, nv=0, nw=0;
  
  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&
  
  //Get the cell energy
  getEnergyScaleFactors(&slope, &intercept);
  cell.e = getEnergy(cell.e, slope, intercept);
  
  //Set the position related variables 
  if (sideIndex==0)
  {
    nAway = cell.nx-1;
    nu=0; nv=cell.ny; nw=cell.nz;
  }
  else if (sideIndex==1)
  {
    nAway = NCELL-cell.nx;
    nu=NCELL+1; nv=cell.ny; nw=cell.nz;
  }
  else if (sideIndex==2)
  {
    nAway = cell.ny-1;
    nu=cell.nx; nv=0; nw=cell.nz;
  }
  else if (sideIndex==3)
  {
    nAway = NCELL-cell.ny;
    nu=cell.nx; nv=NCELL+1; nw=cell.nz;
  }
  else if (sideIndex==4)
  {
    nAway = cell.nz-1;
    nu=cell.nx; nv=cell.ny; nw=0;
  }
  else if (sideIndex==5)
  {
    nAway = NCELL-cell.nz;
    nu=cell.nx; nv=cell.ny; nw=NCELL+1;
  }
  else
  {
    printf("Incorrect call to getWaveform!\n");
    nAway=0;
    nu=0; nv=cell.ny; nw=cell.nz;
  }
  
  //&&&&&&&&&&&&&&&&&&&&
  //& Get the Waveform &
  //&&&&&&&&&&&&&&&&&&&&
  
  //Zero the waveform
  for (int j=0; j<N_TIME_BIN_PRINT; j++) { local[j]=0.0; }

  //Set the temp waveform
  if (int(cell.t)>0)        //Start time is greater than 0 in integer steps of the bin size
  {
    for (int j=0; j<N_TIME_BIN_PRINT-1; j++)
    {
      if (j>=int(cell.t)) local[j] = templateWaveform[nAway][j-int(cell.t)]*cell.e;
    }
    
    sum=0.0;
    for (int j=N_TIME_BIN_PRINT-1-int(cell.t); j<N_TIME_BIN_PRINT; j++)
    {
      sum += templateWaveform[nAway][j];
    }
    local[N_TIME_BIN_PRINT-1] = sum*cell.e;
  }
  else if (int(cell.t)==0)  //Start time is 0 in integer steps of the bin size
  {
    for (int j=0; j<N_TIME_BIN_PRINT; j++)
    {
      local[j] = templateWaveform[nAway][j]*cell.e;
    }
  }  
  else if (int(cell.t)<0)   //Start time is less than 0 in integer steps of the bin size
  {
    for (int j=-int(cell.t); j<N_TIME_BIN_PRINT; j++)
    {
      local[j+int(cell.t)] = templateWaveform[nAway][j]*cell.e; //Just assume that the tail integral
    }                                                           //bin is in the same relative 
  }                                                             //position and zero counts after
  else
  {
    printf("There is an problem with the cell start time!\n");
    return;
  }

  //Non-Integral start time adjustment
  //  For a better time reconstruction consider that all counts in a template waveform bin are
  //  uniformly distributed in that bin. Therefore, if you have a non-integer start time then you
  //  can adjust the binning the local waveform

  if (TIME_RECON)
  {
    double saveLocal[N_TIME_BIN_PRINT];
    for (int j=0; j<N_TIME_BIN_PRINT; j++) { saveLocal[j] = local[j]; }
            
    double p = cell.t - int(cell.t);
    double q = 1.0 - p;
    
    local[0] = q*saveLocal[0];
    for (int j=1; j<N_TIME_BIN_PRINT-1; j++)
    {
      local[j] = p*saveLocal[j-1] + q*saveLocal[j];
    }
    local[N_TIME_BIN_PRINT-1] = p*saveLocal[N_TIME_BIN_PRINT-2] + saveLocal[N_TIME_BIN_PRINT-1];
                  
  }
  
  //Get the pmt index and add the waveform
  pmtIndex = getComponentIndex(pmtList, NCELL*NCELL*6, nu, nv, nw);
  if (pmtIndex==-1) printf("Error in getWaveform!\n");    
    
  for (int j=0; j<N_TIME_BIN_PRINT; j++) { waveform[pmtIndex][j] += local[j]; }   
}

//**************************************************************************************************

void copyComponents(component_t * in, component_t * out, int len)
{
  for (int i=0; i<len; i++) 
  {
    out[i].nx = in[i].nx;  out[i].ny = in[i].ny;  out[i].nz = in[i].nz;
    out[i].e = in[i].e; out[i].t = in[i].t;
  }
}

//**************************************************************************************************

int getHighEnergyCellsPmtIndex(component_t * pmtList, component_t * cellList, int pmtCnt, 
  int cellCnt)
{
  //Declare needed variables
  double highE=0.0;
  int nv, nw;
  int pmtIndex;
   
  //Get the position of the highest energy cell
  for (int i=0; i<cellCnt; i++)
  {
    if (cellList[i].e>highE) highE = cellList[i].e;
  }
  
  for (int i=0; i<cellCnt; i++)
  {
    if (cellList[i].e==highE)
    {
      nv = cellList[i].ny;
      nw = cellList[i].nz;
      break;
    }
  }
        
  //Get the index for one of the PMTs that views the highest cell
  for (int i=0; i<pmtCnt; i++)
  {
    if (pmtList[i].nx==0 && pmtList[i].ny==nv && pmtList[i].nz==nw) //Just choose the x=0 face
    {                                                               //PMT
      pmtIndex=i;
      break;
    }
  }
  
  return pmtIndex;
}

//**************************************************************************************************

double getTotalCharge(component_t * c, int len)
{
  double sum=0.0;
  for (int i=0; i<len; i++) { sum += c[i].e; }
  return sum;
}

//**************************************************************************************************

void printDetectorPmts(detector_t * det, int chargeOrTimeFlag)
{
  chargeOrTimeFlag = 0;
  
  if (chargeOrTimeFlag==0) { printf("Printing the charge\n"); }
  else { printf("Printing the time\n"); }
  
  int stop = (1-MIRROR_FLAG)*(NCELL+1);
  for (int i=0; i<=stop; i+=NCELL+1) {
    for (int j=1; j<=NCELL; j++) {
      for (int k=1; k<=NCELL; k++) {
        printf("%d %d %d %lf\n", i, j, k, det->component[i][j][k] );
        printf("%d %d %d %lf\n", j, i, k, det->component[j][i][k] );
        printf("%d %d %d %lf\n", j, k, i, det->component[j][k][i] );
  } } }
  printf("\n");
  
}

//**************************************************************************************************
