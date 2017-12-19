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

//&&&&&&&&&&&&&&&&&&
//& Custom Headers &
//&&&&&&&&&&&&&&&&&&
#include "reconLib.hh"
#include "parsing.hh"

//&&&&&&&&&&&&&&&&&&&
//& Data Structures &
//&&&&&&&&&&&&&&&&&&&
typedef struct hexColor_t
{
  int redValue;
  int greenValue;
  int blueValue;
}hexColor_t;

typedef struct colorMap_t
{
  hexColor_t startColor;
  hexColor_t stopColor;
}colorMap_t;

//&&&&&&&&&&&&&
//& Constants &
//&&&&&&&&&&&&&
#define MIN_PMT_CHARGE 0.000  //The minimum normalized charge
#define MIN_CELL_ENERGY 0.000 //The minimum normalized energy

//&&&&&&&&&&&&&&&&&&&&
//& Global Variables &
//&&&&&&&&&&&&&&&&&&&&

//Some useful colors to pre-define
hexColor_t RED = {255, 0, 0};                   //        
hexColor_t ORANGE = {255, 127, 0};              //
hexColor_t YELLOW = {255, 255, 0};              //
hexColor_t CHARTREUSE = {127, 255, 0};          //
hexColor_t GREEN = {0, 255, 0};                 //
hexColor_t SPRING_GREEN = {0, 255, 127};        //
hexColor_t CYAN = {0, 255, 255};                //
hexColor_t BLUE = {0, 0, 255};                  //
hexColor_t VIOLET = {143, 0, 255};              //
hexColor_t WHITE = {255, 255, 255};             //
hexColor_t BLACK = {0, 0, 0};                   //
hexColor_t NO_HIT_GRAY_PMT = {211, 211, 211};   //The original gray used in Bruce's code
hexColor_t NO_HIT_GRAY_CELL = {147, 147, 147};  //

//Some useful color maps to pre-define
colorMap_t ORIGINAL_PMT = {NO_HIT_GRAY_PMT, RED};           //Bruce's original PMT color map
colorMap_t ORIGINAL_CELL = {NO_HIT_GRAY_PMT, BLUE};         //Bruce's original PMT color map
colorMap_t ORIGINAL_PRIME_CELL = {NO_HIT_GRAY_CELL, BLUE};  //The low energy cells go to a darker gray for visual clarity
colorMap_t AUTUMN = {YELLOW, RED};                          //Based off of the MATLAB autumn color map
colorMap_t WINTER = {SPRING_GREEN, BLUE};                   //Based off of the MATLAB winter color map

//TO DO: Add other MATLAB and gnuplot color maps

//&&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Declarations &
//&&&&&&&&&&&&&&&&&&&&&&&&&

int moveToEvent(FILE * event, int eventNumber);
void readEventFileBlock(FILE * event);
void plotTheEvent(char * plotScriptName, char * terminalName, char * scaleType, char * colorMapName,
  double pmtHalfSpacing, double cellHalfSpacing, component_t * pmtList, component_t * cellList,
    int nCells, int nPmts, int plotType, int whatToPlotFlag);
hexColor_t getColor(int pmtOrCell, int scaleTypeFlag, colorMap_t map, double q, int plotType);
void getGnuplotTerminalFileExtension(char * terminalName, char * fileExtension);

//*************************************************************************************************

//&&&&&&&&
//& Main &
//&&&&&&&&
int main ()
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Declare Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  //Configuration variables
  char terminalName[stringLen_g];     //
  char scaleType[stringLen_g];        //
  char colorMapName[stringLen_g];     //
  char plotType[stringLen_g];         //
  double pmtSize, cellSize;           //
  int plotTypeFlag;                   //
  
  //Variables for running the program
  char string0[stringLen_g];          //
  char string1[stringLen_g];          //
  char eventFileName[stringLen_g];    //
  char plotScriptName[stringLen_g];   //
  component_t * pmtList, * cellList;  //
  pmtTimeBin_t * hitBinList;          //
  double pmtHalfSpacing;              //
  double cellHalfSpacing;             //
  double td1;                         //Temp variables
  int whatToPlotFlag;                 //
  int checkNum, numberOfEvents;       //
  int eventNumber, reconOrTrueFlag;   //
  int checkSuccess;                   //
  int numberOfBlocks;                 //
  int nLines, nPmts, nCells, nBins;   //
  
  //Files
  FILE * in, * config, * script, * event; //


  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&
  
  //Open and read the configuration file
  config = fopen("inputs/eventDisplayConfig.txt", "r");
  if (config==NULL)
  {
    printf("Could not open the file inputs/eventDisplayConfig.txt! Exiting...\n");
    exit(1);
  }
  
  initString(terminalName, stringLen_g);
  initString(plotType, stringLen_g);
  initString(scaleType, stringLen_g);
  initString(colorMapName, stringLen_g);
  
  checkNum = getNumEntries(config);
  if (checkNum!=7)
  {
    printf("The configuration file is not properly formatted. Using the defaults...\n");
    
    sprintf(terminalName, "default");
    sprintf(plotType, "energy");
    sprintf(scaleType, "sqrt");
    sprintf(colorMapName, "winter");
    pmtSize = 0.9;
    cellSize = 0.4;
    whatToPlotFlag = 0;
  }
  
  getEntry(config, string0, stringLen_g);
  sscanf(string0, "%s", terminalName);
  
  getEntry(config, string0, stringLen_g);
  sscanf(string0, "%s", plotType);
  
  getEntry(config, string0, stringLen_g);
  sscanf(string0, "%s", scaleType);
  
  getEntry(config, string0, stringLen_g);
  sscanf(string0, "%s", colorMapName);
  
  getEntry(config, string0, stringLen_g);
  sscanf(string0, "%lf", &pmtSize);
  
  getEntry(config, string0, stringLen_g);
  sscanf(string0, "%lf", &cellSize);

  getEntry(config, string0, stringLen_g);
  sscanf(string0, "%d", &whatToPlotFlag);
  
  fclose(config);
  
  pmtHalfSpacing = (1.0 - pmtSize)/2.0;
  cellHalfSpacing = (1.0 - cellSize)/2.0;
  
  if (strstr(plotType, "energy")!=NULL) { plotTypeFlag = 0; }
  else { plotTypeFlag = 1; }
  
  if (whatToPlotFlag<0 || whatToPlotFlag>2) { whatToPlotFlag = 0; }
    
  //Open the input file 
  in = fopen("inputs/eventDisplay.in", "r");
  if (in==NULL)
  {
    printf("Could not open the file inputs/eventDisplay.in! Exiting...\n");
    exit(1);
  } 
  numberOfEvents = getNumEntries(in);
  
  //Open and set the script file
  script = fopen("generateEventDisplay.sh", "w");
  if (script==NULL)
  {
    printf("Could not open the file generateEventDisplay.sh! Exiting...\n");
    fclose(in);
    exit(1);
  }
  
  fprintf(script, "cd eventDisplayFiles\n");


  //&&&&&&&&&&&&&&&&&&&
  //& Read the Events &
  //&&&&&&&&&&&&&&&&&&&
  
  for (int i=0; i<numberOfEvents; i++)
  {
    //Get the needed info
    getEntry(in, string0, stringLen_g);
    
    initString(eventFileName, stringLen_g);
    sscanf(string0, "%s %d %d", eventFileName, &eventNumber, &reconOrTrueFlag);
    
    //Open the event file
    event = fopen(eventFileName, "r");
    if (event==NULL)
    {
      printf("Could not open the file %s! Skipping to the next entry...\n", eventFileName);
      continue;
    }
    
    //Go to the desired event
    checkSuccess = moveToEvent(event, eventNumber);
    
    if (checkSuccess==0)
    {
      printf("There was an error in reading the event file! Skipping to the next entry...\n");
      fclose(event);
      continue;
    }
    
    //Get the desired pmt and cells
    if (reconOrTrueFlag==0) { numberOfBlocks = 4; }
    else { numberOfBlocks = 3; }
    
    for (int j=0; j<numberOfBlocks; j++) { readEventFileBlock(event); }
    
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, event);
    
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, event); 
    sscanf(string0, "%d", &nLines);      
    
    if (reconOrTrueFlag==0) 
    {
      nPmts = nLines;
      pmtList = new component_t [nPmts];
      clearComponents(pmtList, nPmts);
    }
    else
    {
      nCells = nLines;
      cellList = new component_t [nCells];
      clearComponents(cellList, nCells);
    }
    
    for (int j=0; j<nLines; j++)
    {
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event);  
      
      if (reconOrTrueFlag==0)
      {
        sscanf(string0, "%d %d %d %lf", &pmtList[j].nx, &pmtList[j].ny, &pmtList[j].nz, &pmtList[j].e);
      }
      else
      {
        sscanf(string0, "%d %d %d %lf %lf %lf", &cellList[j].nx, &cellList[j].ny, &cellList[j].nz,
          &td1, &cellList[j].t, &cellList[j].e);
      }   
    }
    
    if (reconOrTrueFlag==0 && plotTypeFlag==0) { numberOfBlocks = 1; }
    else { numberOfBlocks = 0; }
    
    for (int j=0; j<numberOfBlocks; j++) { readEventFileBlock(event); }
    
    if (plotTypeFlag==1 && reconOrTrueFlag==0)  //Save the PMT waveforms if needed
    {
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event);
      
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event); 
      sscanf(string0, "%d", &nLines); 
      
      nBins = nLines;
      hitBinList = new pmtTimeBin_t [nBins];
      
      for (int entry=0; entry<nLines; entry++)
      {
        initString(string0, stringLen_g);
        fgets(string0, stringLen_g, event); 
        sscanf(string0, "%d %d %d %d %lf", &hitBinList[entry].nu, &hitBinList[entry].nv,
          &hitBinList[entry].nw, &hitBinList[entry].bin, &hitBinList[entry].count );       
      }
      
    } //End PMT hit times
    
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, event);
    
    initString(string0, stringLen_g);
    fgets(string0, stringLen_g, event); 
    sscanf(string0, "%d", &nLines);      
    
    if (reconOrTrueFlag==0) 
    {
      nCells = nLines;
      cellList = new component_t [nCells];
      clearComponents(cellList, nCells);
    }
    else
    {
      nPmts = nLines;
      pmtList = new component_t [nPmts];
      clearComponents(pmtList, nPmts);
    }
    
    for (int j=0; j<nLines; j++)
    {
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event);  
      
      if (reconOrTrueFlag==0)
      {
        sscanf(string0, "%d %d %d %lf %lf", &cellList[j].nx, &cellList[j].ny, &cellList[j].nz,
          &cellList[j].t, &cellList[j].e);
      }
      else
      {
        sscanf(string0, "%d %d %d %lf", &pmtList[j].nx, &pmtList[j].ny, &pmtList[j].nz, &pmtList[j].e);
      }   
    }
    
    if (plotTypeFlag==1 && reconOrTrueFlag==1)  //Save the PMT waveforms if needed
    {
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event);
      
      initString(string0, stringLen_g);
      fgets(string0, stringLen_g, event); 
      sscanf(string0, "%d", &nLines); 
      
      nBins = nLines;
      hitBinList = new pmtTimeBin_t [nBins];
      
      for (int entry=0; entry<nLines; entry++)
      {
        initString(string0, stringLen_g);
        fgets(string0, stringLen_g, event); 
        sscanf(string0, "%d %d %d %d %lf", &hitBinList[entry].nu, &hitBinList[entry].nv,
          &hitBinList[entry].nw, &hitBinList[entry].bin, &hitBinList[entry].count );       
      }
      
    } //End PMT hit times

    fclose(event);

    //Get the PMT hit times if needed
    if (plotTypeFlag==1)
    {
      for (int j=0; j<nPmts; j++)
      {
        int firstHitBin = N_TIME_BIN_PRINT + 1;
        for (int k=0; k<nBins; k++)
        {
          if (pmtList[j].nx==hitBinList[k].nu && pmtList[j].ny==hitBinList[k].nv &&
            pmtList[j].nz==hitBinList[k].nw )
          {
            if (hitBinList[k].bin<firstHitBin) { firstHitBin = hitBinList[k].bin; }
          }
        }
        
        pmtList[j].t = (firstHitBin - 1)*370.0/1000.0;
      }
    }

    //Get the plot script name
    initString(plotScriptName, stringLen_g);
    sprintf(plotScriptName, "eventDisplayFiles/");
    
    initString(string0, stringLen_g);
    copyString(string0, eventFileName, stringLen_g);
    removePath(string0, 0);
    removeFileExtension(string0, 4);
    
    initString(string1, stringLen_g);
    if (reconOrTrueFlag==0) { sprintf(string1, "_%dPlot", eventNumber); }
    else { sprintf(string1, "_%dTruePlot", eventNumber); }
    
    strcat(string0, string1);
    strcat(plotScriptName, string0);
    
    initString(string1, stringLen_g);
    copyString(string1, plotScriptName, stringLen_g);
    removePath(string1, 0);
    
    //Add line to the script
    if (strstr(terminalName, "default")==NULL)
    {
      fprintf(script, "gnuplot %s\n", string1);
    }
    else
    {
      fprintf(script, "gnuplot -persist %s\n", string1);
    }
    
    //Plot the event
    plotTheEvent(plotScriptName, terminalName, scaleType, colorMapName, pmtHalfSpacing,
      cellHalfSpacing, pmtList, cellList, nCells, nPmts, plotTypeFlag, whatToPlotFlag);

    //Clean-up
    delete [] pmtList;
    delete [] cellList;
    if (plotTypeFlag==1) delete [] hitBinList;
    
  } //End loop over the events to plot


  //&&&&&&&
  //& End &
  //&&&&&&&
  
  fprintf(script, "cd ../\n");
  fprintf(script, "exit\n");
  
  //Close files
  fclose(in);
  fclose(script);
  
  printf("Run the script generateEventDisplay.sh to make the plots.\n");

  return 0;
}

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&

int moveToEvent(FILE * event, int eventNumber)
{
  //Declare needed variables
  char string[stringLen_g];
  int successFlag=0;
  int nLines;
  int currentEventNumber = eventNumber + 1;
  
  //Start fresh
  rewind(event);
  
  //Loop until the correct event is found
  initString(string, stringLen_g);  
  while (currentEventNumber!=eventNumber && fgets(string, stringLen_g, event)!=NULL)
  {
    //Check the event header
    if (strstr(string, "#Event")==NULL)
    {
      return successFlag;
    }
    
    initString(string, stringLen_g);
    fgets(string, stringLen_g, event); 
    sscanf(string, "%d", &nLines);
    
    for (int i=0; i<nLines; i++)
    {
      initString(string, stringLen_g);
      fgets(string, stringLen_g, event);
      if (i==0) sscanf(string, "%d", &currentEventNumber);
    }  
    
    //Check if at destination
    if (currentEventNumber!=eventNumber)
    {
      for (int i=0; i<7; i++)
      {
        readEventFileBlock(event);
      }
    }
    else { successFlag=1; }
  
    initString(string, stringLen_g);  
  } //End search loop

  //End
  return successFlag;  
}

//*************************************************************************************************

void readEventFileBlock(FILE * event)
{
  char string[stringLen_g];
  int nLines;
  
  initString(string, stringLen_g);
  fgets(string, stringLen_g, event);
  
  initString(string, stringLen_g);    
  fgets(string, stringLen_g, event);
  sscanf(string, "%d", &nLines);

  for (int i=0; i<nLines; i++)
  {
    initString(string, stringLen_g);
    fgets(string, stringLen_g, event);
  }
}

//*************************************************************************************************

void plotTheEvent(char * plotScriptName, char * terminalName, char * scaleType, char * colorMapName,
  double pmtHalfSpacing, double cellHalfSpacing, component_t * pmtList, component_t * cellList,
    int nCells, int nPmts, int plotType, int whatToPlotFlag)
{
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //& Declared Needed Variables &
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  colorMap_t map;                     //
  char plotName[stringLen_g];         //
  char fileExtension[stringLen_g];    //
  double conversionSlope;             //
  double conversionIntercept;         //
  double x, y, z, x2, y2, z2, delta;  //
  double cellThres, cellEnergy;       //
  hexColor_t pmtColor, cellColor;     //
  int scaleTypeFlag, polygonCount=0;  //
  int checkComponent;                 //
  FILE * plot;                        //
  
  
  //&&&&&&&&&&&&&&
  //& Initialize &
  //&&&&&&&&&&&&&&
  
  //Open and begin the plot script
  plot = fopen(plotScriptName, "w");
  if (plot==NULL)
  {
    printf("Could not open the file eventDisplay/eventPlot! Exiting function...\n");
    return;
  }

  if (strstr(terminalName, "default")==NULL)
  {
    fprintf(plot, "set terminal %s\n", terminalName);
    
    initString(plotName, stringLen_g);
    copyString(plotName, plotScriptName, stringLen_g);
    removePath(plotName, 0);
    removeFileExtension(plotName, 4);
    
    getGnuplotTerminalFileExtension(terminalName, fileExtension);
    strcat(plotName, fileExtension);
    
    fprintf(plot, "set output '%s'\n", plotName);
  }
  
  //Get the scale function
  
  //TO DO: Actually add the linear and the log functionality
  
  if (strstr(scaleType, "power10")!=NULL) { scaleTypeFlag = 0; }
  else if (strstr(scaleType, "lin")!=NULL) { scaleTypeFlag = 1; }
  else if (strstr(scaleType, "sqrt")!=NULL) { scaleTypeFlag = 2; }
  else if (strstr(scaleType, "log")!=NULL) { scaleTypeFlag = 3; }
  else { scaleTypeFlag = 2; } //Make sqrt the default
  
  //Get the color map name
  if (strstr(colorMapName, "original")!=NULL) { map = ORIGINAL_PMT; }
  else if (strstr(colorMapName, "originalPrime")!=NULL) { map = ORIGINAL_PMT; }
  else if (strstr(colorMapName, "autumn")!=NULL) { map = AUTUMN; }
  else if (strstr(colorMapName, "winter")!=NULL) { map = WINTER; }
  else { map = WINTER; }

  //Get the maximum charge PMT and normalize to it
  sortComponents(pmtList, nPmts, plotType);
  if (plotType==0)
  {
    for (int i=1; i<nPmts; i++)
    {
      pmtList[i].e = pmtList[i].e/pmtList[0].e;
      if (pmtList[i].e<MIN_PMT_CHARGE) {pmtList[i].e = 0.0; }    
    }
    pmtList[0].e = 1.0;
  }
  else
  {
    for (int i=nPmts-2; i>=0; i--)
    {
      pmtList[i].t = pmtList[i].t/pmtList[nPmts-1].t;
    }
    pmtList[nPmts-1].t = 1.0;  
  }
    
  //Get the vertex cell and normalize to it
  getEnergyScaleFactors(&conversionSlope, &conversionIntercept);
  
  if (CELL_THRESHOLD==0) { cellThres = 0.0; }
  else if (CELL_THRESHOLD==1)
  {
    cellThres = conversionSlope + conversionIntercept;
    if (cellThres<0.0)
    {
      printf("The cell hit threshold was negative! Setting to %lf!\n", E_MIN);
      cellThres = E_MIN;
    }
  }
  else { cellThres = E_MIN; }  
  
  sortComponents(cellList, nCells, plotType);
  if (plotType==0)
  {
    for (int i=1; i<nCells; i++)
    {
      if (plotType==0 && cellList[i].e<cellThres) {cellList[i].e = 0.0; }  //Ensures that all hit cells have the same threshold
      cellList[i].e = cellList[i].e/cellList[0].e;
      if (cellList[i].e<MIN_CELL_ENERGY) {cellList[i].e = 0.0; }
    }
    cellList[0].e = 1.0;
  }
  else
  {
    for (int i=nCells-2; i>=0; i--)
    {
      cellList[i].t = cellList[i].t/cellList[nCells-1].t;
    }
    cellList[nCells-1].t = 1.0;  
  }

  //&&&&&&&&&&&&&&&&&&&
  //& Paint the Sides &
  //&&&&&&&&&&&&&&&&&&&
  
  //TO DO: Figure out how to fix the "box flipping" optical illusion
  if (whatToPlotFlag<=1)
  {
    delta = 1.0 - 2.0*pmtHalfSpacing;
    for (int i=0; i<NCELL; i++)
    {
      for (int j=0; j<NCELL; j++)
      {
        //The x=0 side
        x = -pmtHalfSpacing;
        y = i + pmtHalfSpacing;
        z = j + pmtHalfSpacing;
        
        polygonCount++;
        fprintf(plot, "set object %d \\\n", polygonCount);
        fprintf(plot, "polygon from \\\n");
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y+delta, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y+delta, z+delta);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y, z+delta);
        fprintf(plot, "%lf, %lf, %lf \\\n", x, y, z);
        
        checkComponent = getComponentIndex(pmtList, nPmts, 0, i+1, j+1);
        if (checkComponent==-1)
        {
          pmtColor = NO_HIT_GRAY_PMT;
        }
        else
        {
          if (plotType==0)
          {
            pmtColor = getColor(0, scaleTypeFlag, map, pmtList[checkComponent].e, plotType);
          }
          else
          {
            pmtColor = getColor(0, scaleTypeFlag, map, pmtList[checkComponent].t, plotType);        
          }
        }
        
        fprintf(plot, "fillstyle solid noborder fillcolor rgb '#%02x%02x%02x'\n", pmtColor.redValue, 
          pmtColor.greenValue, pmtColor.blueValue);
          
        //y=0 side
        x = i + pmtHalfSpacing;
        y = -pmtHalfSpacing;
        z = j + pmtHalfSpacing;
              
        polygonCount++;
        fprintf(plot, "set object %d \\\n", polygonCount);
        fprintf(plot, "polygon from \\\n");
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x+delta, y, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x+delta, y, z+delta);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y, z+delta);
        fprintf(plot, "%lf, %lf, %lf \\\n", x, y, z);

        checkComponent = getComponentIndex(pmtList, nPmts, i+1, 0, j+1);
        if (checkComponent==-1)
        {
          pmtColor = NO_HIT_GRAY_PMT;
        }
        else
        {
          if (plotType==0)
          {
            pmtColor = getColor(0, scaleTypeFlag, map, pmtList[checkComponent].e, plotType);
          }
          else
          {
            pmtColor = getColor(0, scaleTypeFlag, map, pmtList[checkComponent].t, plotType);        
          }
        }

        fprintf(plot, "fillstyle solid noborder fillcolor rgb '#%02x%02x%02x'\n", pmtColor.redValue, 
          pmtColor.greenValue, pmtColor.blueValue); 
          
        //z=0 side
        x = i + pmtHalfSpacing;
        y = j + pmtHalfSpacing;
        z = -pmtHalfSpacing;
              
        polygonCount++;
        fprintf(plot, "set object %d \\\n", polygonCount);
        fprintf(plot, "polygon from \\\n");
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x+delta, y, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x+delta, y+delta, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y+delta, z);
        fprintf(plot, "%lf, %lf, %lf \\\n", x, y, z);

        checkComponent = getComponentIndex(pmtList, nPmts, i+1, j+1, 0);
        if (checkComponent==-1)
        {
          pmtColor = NO_HIT_GRAY_PMT;
        }
        else
        {
          if (plotType==0)
          {
            pmtColor = getColor(0, scaleTypeFlag, map, pmtList[checkComponent].e, plotType);
          }
          else
          {
            pmtColor = getColor(0, scaleTypeFlag, map, pmtList[checkComponent].t, plotType);        
          }
        }
        
        fprintf(plot, "fillstyle solid noborder fillcolor rgb '#%02x%02x%02x'\n", pmtColor.redValue, 
          pmtColor.greenValue, pmtColor.blueValue); 
      } //End second loop
    } //End first loop
  }

  //&&&&&&&&&&&&&&&&&&&&&&&
  //& Paint the Hit Cells &
  //&&&&&&&&&&&&&&&&&&&&&&&

  if (whatToPlotFlag%2==0)
  {
    if (strstr(colorMapName, "original")!=NULL) { map = ORIGINAL_CELL; }
    else if (strstr(colorMapName, "originalPrime")!=NULL) { map = ORIGINAL_PRIME_CELL; }
    
    
    for (int i=nCells-1; i>=0; i--)
    {
      int paintTheCell=0;
      cellEnergy = cellList[i].e;
      if ( (plotType==0 && cellEnergy>=MIN_CELL_ENERGY) || plotType==1) { paintTheCell = 1; }
     
      if (paintTheCell==1)
      {
        if (plotType==0) { cellColor = getColor(1, scaleTypeFlag, map, cellList[i].e, plotType); }
        else  { cellColor = getColor(1, scaleTypeFlag, map, cellList[i].t, plotType); }

        x = cellList[i].nx - 1.0 + cellHalfSpacing;
        y = cellList[i].ny - 1.0 + cellHalfSpacing;
        z = cellList[i].nz - 1.0 + cellHalfSpacing;
        
        x2 = cellList[i].nx - cellHalfSpacing;
        y2 = cellList[i].ny - cellHalfSpacing;
        z2 = cellList[i].nz - cellHalfSpacing;
        
        polygonCount++;
        fprintf(plot, "set object %d \\\n", polygonCount);
        fprintf(plot, "polygon from \\\n");
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y2, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y2, z2);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y, z2);
        fprintf(plot, "%lf, %lf, %lf \\\n", x, y, z);
        
        fprintf(plot, "fillstyle solid noborder fillcolor rgb '#%02x%02x%02x'\n", cellColor.redValue, 
          cellColor.greenValue, cellColor.blueValue);
        
        polygonCount++;
        fprintf(plot, "set object %d \\\n", polygonCount);
        fprintf(plot, "polygon from \\\n");
        fprintf(plot, "%lf, %lf, %lf to \\\n", x2, y, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x2, y2, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x2, y2, z2);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x2, y, z2);
        fprintf(plot, "%lf, %lf, %lf \\\n", x2, y, z);
        
        fprintf(plot, "fillstyle solid noborder fillcolor rgb '#%02x%02x%02x'\n", cellColor.redValue, 
          cellColor.greenValue, cellColor.blueValue); 
          
        polygonCount++;
        fprintf(plot, "set object %d \\\n", polygonCount);
        fprintf(plot, "polygon from \\\n");
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x2, y, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x2, y, z2);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y, z2);
        fprintf(plot, "%lf, %lf, %lf \\\n", x, y, z);
        
        fprintf(plot, "fillstyle solid noborder fillcolor rgb '#%02x%02x%02x'\n", cellColor.redValue, 
          cellColor.greenValue, cellColor.blueValue); 
          
        polygonCount++;
        fprintf(plot, "set object %d \\\n", polygonCount);
        fprintf(plot, "polygon from \\\n");
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y2, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x2, y2, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x2, y2, z2);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y2, z2);
        fprintf(plot, "%lf, %lf, %lf \\\n", x, y2, z);
        
        fprintf(plot, "fillstyle solid noborder fillcolor rgb '#%02x%02x%02x'\n", cellColor.redValue, 
          cellColor.greenValue, cellColor.blueValue); 

        polygonCount++;
        fprintf(plot, "set object %d \\\n", polygonCount);
        fprintf(plot, "polygon from \\\n");
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x2, y, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x2, y2, z);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y2, z);
        fprintf(plot, "%lf, %lf, %lf \\\n", x, y, z);
        
        fprintf(plot, "fillstyle solid noborder fillcolor rgb '#%02x%02x%02x'\n", cellColor.redValue, 
          cellColor.greenValue, cellColor.blueValue); 
          
        polygonCount++;
        fprintf(plot, "set object %d \\\n", polygonCount);
        fprintf(plot, "polygon from \\\n");
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y, z2);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x2, y, z2);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x2, y2, z2);
        fprintf(plot, "%lf, %lf, %lf to \\\n", x, y2, z2);
        fprintf(plot, "%lf, %lf, %lf \\\n", x, y, z2);
        
        fprintf(plot, "fillstyle solid noborder fillcolor rgb '#%02x%02x%02x'\n", cellColor.redValue, 
          cellColor.greenValue, cellColor.blueValue); 
      } //End check if the cell was hit
    } //End loop over the hit cells
  }

  //&&&&&&&&&&&&&&&&&&&
  //& Finish the Plot &
  //&&&&&&&&&&&&&&&&&&&
  
  fprintf(plot, "unset key\n");
  fprintf(plot, "set ticslevel 0\n");
  fprintf(plot, "unset tics\n");
  fprintf(plot, "set border linecolor rgb '#ffffff'\n");
  fprintf(plot, "set view 70,120\n");
  fprintf(plot, "set xrange [0:15]\n");
  fprintf(plot, "set yrange [0:15]\n");
  fprintf(plot, "set zrange [0:15]\n");
  fprintf(plot, "set view equal xyz\n");
  fprintf(plot, "splot -1\n");
  

  //&&&&&&&
  //& End &
  //&&&&&&&

  fclose(plot);
    
  return;
}

//*************************************************************************************************

hexColor_t getColor(int pmtOrCell, int scaleTypeFlag, colorMap_t map, double q, int plotType)
{
  //Declared needed variables
  hexColor_t componentColor;
  colorMap_t localMap;
  double qMin, qMax;
  double mRed, mGreen, mBlue;
  double bRed, bGreen, bBlue;
  
  //Initialize
  if (scaleTypeFlag==0)       //power 10
  {
    q = pow(10.0, q);
    if (pmtOrCell==0) { qMin = pow(10.0, MIN_PMT_CHARGE); }
    else { qMin = pow(10.0, MIN_CELL_ENERGY); }
    qMax = 10.0;  
  }
  else if (scaleTypeFlag==1)  //lin
  {
    if (pmtOrCell==0) { qMin = MIN_PMT_CHARGE; }
    else { qMin = MIN_CELL_ENERGY; }
    qMax = 1.0;  
  }
  else if (scaleTypeFlag==2)  //Sqrt
  {
    q = sqrt(q);
    if (pmtOrCell==0) { qMin = sqrt(MIN_PMT_CHARGE); }
    else { qMin = sqrt(MIN_CELL_ENERGY); }
    qMax = 1.0; 
  }
  else if (scaleTypeFlag==3)  //log
  {
    q = log(q)/log(10.0);
    if (pmtOrCell==0) { qMin = log(MIN_PMT_CHARGE)/log(10.0); }
    else { qMin = log(MIN_CELL_ENERGY)/log(10.0); }
    qMax = log(1.0);  
  }
  
  if (plotType==0)
  {
    localMap = map;
  }
  else
  {
    localMap.startColor = map.stopColor;
    localMap.stopColor = map.startColor;
  }
      
  mRed = (localMap.startColor.redValue - localMap.stopColor.redValue)/(qMin - qMax);
  mGreen = (localMap.startColor.greenValue - localMap.stopColor.greenValue)/(qMin - qMax);
  mBlue = (localMap.startColor.blueValue - localMap.stopColor.blueValue)/(qMin - qMax);
  
  bRed = localMap.startColor.redValue - mRed*qMin;
  bGreen = localMap.startColor.greenValue - mGreen*qMin;
  bBlue = localMap.startColor.blueValue - mBlue*qMin;
    
  //Set the color
  if (q>=qMin)
  {
    componentColor.redValue = mRed*q + bRed;
    componentColor.greenValue = mGreen*q + bGreen;
    componentColor.blueValue = mBlue*q + bBlue;
  }
  else
  {
    componentColor = NO_HIT_GRAY_PMT;
  }

  return componentColor;
}

//*************************************************************************************************

void getGnuplotTerminalFileExtension(char * terminalName, char * fileExtension)
{
  initString(fileExtension, stringLen_g);
  if (strstr(terminalName, "epscairo")!=NULL)
  {
    sprintf(fileExtension, ".eps");
  }
  else if (strstr(terminalName, "gif")!=NULL)
  {
    sprintf(fileExtension, ".gif");
  }
  else if (strstr(terminalName, "jpeg")!=NULL)
  {
    sprintf(fileExtension, ".jpeg");
  }
  else if (strstr(terminalName, "pngcairo")!=NULL)
  {
    sprintf(fileExtension, ".png");
  }
  else if (strstr(terminalName, "pdfcairo")!=NULL)
  {
    sprintf(fileExtension, ".pdf");
  }
  else if (strstr(terminalName, "png")!=NULL)
  {
    sprintf(fileExtension, ".png");
  }
  else if (strstr(terminalName, "postscript")!=NULL)
  {
    sprintf(fileExtension, ".eps");
  }
  else if (strstr(terminalName, "svg")!=NULL)
  {
    sprintf(fileExtension, ".svg");
  }
}

//*************************************************************************************************
