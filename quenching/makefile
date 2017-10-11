# Makefile for compling and linking the guide simulation programs
#
# WARNING: This Makefile is cobbled together from my limitted
# knowledge of the make utility. Z. W. Yokley 2013-11-11
#
# Variables to control Makefile operation
#
# Compiler/Linker variables
#
CC = g++
CFLAGS = -Wall -ggdb
LFLAGS = -L/usr/local/lib
ROOTSYS = /home/zwyokley/root
INCPATH = -I/usr/local/include -I$(ROOTSYS)/include/ -I./include
ROOTLIBS = `root-config --cflags --glibs`
GSLFLAGS  = -lgsl -lgslcblas -lm

# Directory variables
#
IN = include/
SRC = src/
LIB = lib/
BACKUP = ./backups/

# Targets
#
all: lib quench

lib: 3dvector.o crossSection.o deposit.o emProcesses.o gammaRay.o parsing.o randomNumberGenerator.o

# Applications
#
quench: quench.o $(LIB)3dvector.o $(LIB)crossSection.o $(LIB)deposit.o \
  $(LIB)emProcesses.o $(LIB)gammaRay.o $(LIB)parsing.o $(LIB)randomNumberGenerator.o 
	$(CC) $(LFLAGS) quench.o $(LIB)3dvector.o $(LIB)crossSection.o $(LIB)deposit.o \
	$(LIB)emProcesses.o $(LIB)gammaRay.o $(LIB)parsing.o $(LIB)randomNumberGenerator.o  \
	$(GSLFLAGS) -o quench
	
quench.o: quench.cc
	$(CC) $(CFLAGS) $(INCPATH) -c quench.cc

# Library object files
#
3dvector.o: $(IN)3dvector.hh $(SRC)3dvector.cc
	$(CC) $(CFLAGS) $(INCPATH) -c $(SRC)3dvector.cc -o $(LIB)3dvector.o

crossSection.o: $(IN)crossSection.hh $(SRC)crossSection.cc
	$(CC) $(CFLAGS) $(INCPATH) -c $(SRC)crossSection.cc -o $(LIB)crossSection.o

deposit.o: $(IN)deposit.hh $(SRC)deposit.cc
	$(CC) $(CFLAGS) $(INCPATH) -c $(SRC)deposit.cc -o $(LIB)deposit.o

emProcesses.o: $(IN)emProcesses.hh $(SRC)emProcesses.cc
	$(CC) $(CFLAGS) $(INCPATH) -c $(SRC)emProcesses.cc -o $(LIB)emProcesses.o

gammaRay.o: $(IN)gammaRay.hh $(SRC)gammaRay.cc
	$(CC) $(CFLAGS) $(INCPATH) -c $(SRC)gammaRay.cc -o $(LIB)gammaRay.o

parsing.o: $(IN)parsing.hh $(SRC)parsing.cc
	$(CC) $(CFLAGS) $(INCPATH) -c $(SRC)parsing.cc -o $(LIB)parsing.o

randomNumberGenerator.o: $(IN)randomNumberGenerator.hh $(SRC)randomNumberGenerator.cc
	$(CC) $(CFLAGS) $(INCPATH) -c $(SRC)randomNumberGenerator.cc -o $(LIB)randomNumberGenerator.o



# Housekeeping
#
cleanLib:
	rm $(LIB)*.o
	
clean:
	rm *.o

cleanAll:
	rm $(LIB)*.o
	rm *.o

backup: 
	tar cfv $(BACKUP)quench_$(shell date +'%Y%m%d').tar crossSections data include inputs old outputs src *.cc makefile README

tar:
	tar cfv quench_$(shell date +'%Y%m%d').tar crossSections data include inputs old outputs src *.cc makefile README









