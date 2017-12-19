#include "parsing.hh"

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&

//*************************************************************************************************

void initString(char str[], int size)
{
  /*  This function initializes the string passed *
  *  to it to an array of null characters.          */
  for (int i=0; i<size; i++) { str[i] = '\0'; }
  return;
}

//*************************************************************************************************

void  clParser(int argc, char *argv[], int clFlags[])
{
  /*   This function parses the command line arguments passed    *
  * to a main function and assigns values to the array clFlags  */
  for (int i=1; i<argc; i++) { printf("Nothing here yet!\n"); }
  return;
}

//*************************************************************************************************

void printfHeader(char str0[], char str1[], int size)
{
  /*  This function generates a header for the     *
  *    log and output files                        */

  //Declare needed variables
  int delimLen=60;
  time_t rawtime;

  //Initialize
  initString(str0, size);
  initString(str1, size);

  //Set run delimiter
  if (size<delimLen) delimLen = size;

  for (int i=0; i<delimLen; i++) { str0[i] = '*'; }

  //Set date
  time(&rawtime);
  sprintf(str1, "%s", ctime(&rawtime));

  return;
}

//*************************************************************************************************

int getNumEntries(FILE * f)
{
  /*  This function returns the number of files    *
  *    listed in the input file. This requires     *
  *    parsing since the input contains comments    */

  char  s[stringLen_g];        //
  int    nEntries=0;          //

  initString(s, stringLen_g);
  while (fgets(s, stringLen_g, f)!=NULL)
  {
    if (s[0]=='#') continue;
    else nEntries++;

    initString(s, stringLen_g);
  }
  rewind(f);

  return nEntries;
}

//*************************************************************************************************

void getEntry(FILE * f, char str[], int size)
{
  /*  This function returns the next file name in  *
  *    the input file.                             */

  char  s0[stringLen_g];      //

  initString(s0, stringLen_g);
  initString(str, size);

  while (fgets(s0, stringLen_g, f)!=NULL)
  {
    if (s0[0]=='#')
    {
      initString(s0, stringLen_g);
      continue;
    }
    else
    {
      memcpy(str, s0, strlen(s0)-1);
      break;
    }
  }

  return;
}

//*************************************************************************************************

void copyString(char out[], char in[], int size)
{
  initString(out, stringLen_g);
  for (int i=0; i<size; i++) { out[i] = in[i]; }
  return;
}

//*************************************************************************************************

void  getDirName(char out[], char in[], int dirNum)
{
  //Get the directory name of the dirNum directory.
  //For example say that you have the string: path/to/file/file.txt.
  //And you want the name of the top level directory, which is path.
  //You would then set dirNum=0.

  initString(out, stringLen_g);
  int dirCnt = getNumDirectories(in);

  if (dirCnt==0) { return; }
  else
  {
    char names[dirCnt][stringLen_g];
    for (int i=0; i<dirCnt; i++) { initString(names[i], stringLen_g); }

    int pos[dirCnt+1]; pos[0]=-1;
    dirCnt=1;
    for (int i=0; i<int(strlen(in)); i++)
    {
      if (in[i]=='/') { pos[dirCnt]=i; dirCnt++; }
    }

    int outPos=0;
    for (int i=pos[dirNum]+1; i<=pos[dirNum+1]-1; i++)
    {
      out[outPos] = in[i];
      outPos++;
    }
  }

  return;
}

//*************************************************************************************************

int getNumDirectories(char s0[])
{
  int dirCnt=0;
  for (int i=0; i<int(strlen(s0)); i++)
  {
    if (s0[i]=='/') dirCnt++;
  }

  return dirCnt;
}

//*************************************************************************************************

void removePath(char s0[], int keepNumPaths)
{
  //This function removes directory paths from the passed file name.
  //keepNumPaths is the number of paths the you want to keep in the file name;
  //therefore:
  // keepNumPaths = 0 -> remove all of the paths, ie the behavior before modifing
  // keepNumPaths = 1 -> remove all but one of the paths
  // keepNumPaths = 2 -> remove all but two of the paths
  // ...
  // keepNumPaths = dirCnt-2 -> remove the top two directories
  // keepNumPaths = dirCnt-1 -> remove the top directory
  // keepNumPaths = dirCnt -> why did you call this function?

  //Declare variables
  char s1[stringLen_g], s2[stringLen_g];
  int dirCnt=0, dirPos;

  //Initialize
  initString(s1, stringLen_g);  initString(s2, stringLen_g);  

  dirCnt = getNumDirectories(s0);
  
  if (keepNumPaths<0) keepNumPaths=0;
  if (keepNumPaths>dirCnt) keepNumPaths=dirCnt;

  dirCnt = dirCnt - keepNumPaths;

  //Cut off the directory names
  copyString(s1, s0, strlen(s0));
  for (int i=0; i<dirCnt; i++)
  {
    //Find the position of the end of the directory
    for (int j=0; j<int(strlen(s0)); j++)
    {
      if (s1[j]=='/') 
      {
        dirPos=j;
        break;
      }
    }

    //Copy just the characters after the position
    for (int j=dirPos+1; j<int(strlen(s1)); j++)
    {
      s2[j-dirPos-1] = s1[j];
    }

    //reinitialize
    initString(s1, stringLen_g);
    copyString(s1, s2, strlen(s2));
    initString(s2, stringLen_g);
  }

  initString(s0, stringLen_g);
  copyString(s0, s1, strlen(s1));

  return;
}

//*************************************************************************************************

void removeFileExtension(char s0[], int extLen)
{
  //Declare variables
  char s1[stringLen_g];
  
  //Remove the file extension
  copyString(s1, s0, strlen(s0)-extLen);
  copyString(s0, s1, strlen(s1));
  
  return;
}

//*************************************************************************************************
