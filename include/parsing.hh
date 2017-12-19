#ifndef PARSING_H
#define PARSING_H

//&&&&&&&&&&&&&
//& C Headers &
//&&&&&&&&&&&&&
#include <time.h>
#include <string.h>

//&&&&&&&&&&&&&&&
//& C++ Headers &
//&&&&&&&&&&&&&&&
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
using namespace std;

//&&&&&&&&&&&&&&&&&&&&
//& Global Variables &
//&&&&&&&&&&&&&&&&&&&&

const int stringLen_g=200;  //Length of strings to use

//&&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Declarations &
//&&&&&&&&&&&&&&&&&&&&&&&&&

void  initString(char str[], int size);
void  clParser(int argc, char *argv[], int clFlags[]);
void  printfHeader(char str0[], char str1[], int size);
int    getNumEntries(FILE * f);
void  getEntry(FILE * f, char str[], int size);
void  copyString(char out[], char in[], int size);
void  getDirName(char out[], char in[], int dirNum);
int   getNumDirectories(char s0[]);
void  removePath(char s0[], int keepNumPaths);
void  removeFileExtension(char s0[], int extLen);

#endif
