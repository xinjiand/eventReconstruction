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
FOR = gfortran
FFLAGS = -ggdb -ffixed-line-length-none
CFLAGS = -Wall -ggdb
LFLAGS = -L/usr/local/lib
ROOTSYS = /home/zwyokley/root
INCPATH = -I/usr/local/include -I$(ROOTSYS)/include/ -I./include
ROOTLIBS = `root-config --cflags --glibs`
GSLFLAGS  = -lgsl -lgslcblas -lm
BACKUP = ./backups/
#
# File variables
#

#
# Directory variables
#
IN = include/
SRC = src/
LIB = lib/
INPUT = inputs/

# Targets
#
all: lib readGeantData IBD_Generator trans scriptGenerator addToEventFile removeFromEventFile \
  combineIBD_Events basisGenerator addBasis chargeReconstruction reconstruct reconstructTeflon \
  splitIbdForReconstruction eventListGenerator getTriggerEvents recordTriggerEvents eventDisplay \
  makeIterationsFile cleanAll

lib: analysis.o detectorParameters.o parsing.o randomNumberGenerator.o reconLib.o

# Applications
#
IBD_Generator: IBD_Generator.o $(LIB)detectorParameters.o $(LIB)parsing.o \
  $(LIB)randomNumberGenerator.o $(LIB)reconLib.o
	$(CC) $(LFLAGS) IBD_Generator.o $(LIB)detectorParameters.o $(LIB)parsing.o $(LIB)randomNumberGenerator.o $(LIB)reconLib.o \
	$(GSLFLAGS) -o IBD_Generator

readGeantData: readGeantData.o $(LIB)detectorParameters.o $(LIB)parsing.o \
  $(LIB)randomNumberGenerator.o $(LIB)reconLib.o
	$(CC) $(LFLAGS) readGeantData.o $(LIB)detectorParameters.o $(LIB)parsing.o $(LIB)randomNumberGenerator.o $(LIB)reconLib.o \
	$(GSLFLAGS) -o readGeantData
	  
trans: trans.f
	$(FOR) $(FFLAGS) -o trans trans.f

scriptGenerator: scriptGenerator.o $(LIB)parsing.o $(LIB)reconLib.o
	$(CC) $(LFLAGS) scriptGenerator.o $(LIB)parsing.o $(LIB)reconLib.o $(GSLFLAGS) -o scriptGenerator
	
addToEventFile: addToEventFile.o $(LIB)parsing.o
	$(CC) $(LFLAGS) addToEventFile.o $(LIB)parsing.o $(GSLFLAGS) -o addToEventFile

removeFromEventFile: removeFromEventFile.o $(LIB)parsing.o
	$(CC) $(LFLAGS) removeFromEventFile.o $(LIB)parsing.o $(GSLFLAGS) -o removeFromEventFile

combineIBD_Events: combineIBD_Events.o $(LIB)parsing.o
	$(CC) $(LFLAGS) combineIBD_Events.o $(LIB)parsing.o $(GSLFLAGS) -o combineIBD_Events
	
splitIbdForReconstruction: splitIbdForReconstruction.o $(LIB)parsing.o
	$(CC) $(LFLAGS) splitIbdForReconstruction.o $(LIB)parsing.o $(GSLFLAGS) -o splitIbdForReconstruction
	
basisGenerator: basisGenerator.o $(LIB)detectorParameters.o $(LIB)parsing.o
	$(CC) $(LFLAGS) basisGenerator.o $(LIB)detectorParameters.o $(LIB)parsing.o $(GSLFLAGS) \
	-o basisGenerator
	
addBasis: addBasis.o $(LIB)reconLib.o $(LIB)parsing.o
	$(CC) $(LFLAGS) addBasis.o $(LIB)reconLib.o $(LIB)parsing.o $(GSLFLAGS) \
	-o addBasis

reconstruct: reconstruct.o $(LIB)detectorParameters.o $(LIB)reconLib.o $(LIB)parsing.o $(LIB)randomNumberGenerator.o
	$(CC) $(LFLAGS) reconstruct.o $(LIB)detectorParameters.o $(LIB)reconLib.o \
	$(LIB)parsing.o $(LIB)randomNumberGenerator.o $(GSLFLAGS) \
	-o reconstruct

reconstructTeflon: reconstructTeflon.o $(LIB)detectorParameters.o $(LIB)reconLib.o $(LIB)parsing.o $(LIB)randomNumberGenerator.o
	$(CC) $(LFLAGS) reconstructTeflon.o $(LIB)detectorParameters.o $(LIB)reconLib.o \
	$(LIB)parsing.o $(LIB)randomNumberGenerator.o $(GSLFLAGS) \
	-o reconstructTeflon

chargeReconstruction: chargeReconstruction.f
	$(FOR) $(FFLAGS) -o chargeReconstruction chargeReconstruction.f

eventListGenerator: eventListGenerator.o $(LIB)analysis.o  $(LIB)detectorParameters.o \
  $(LIB)parsing.o $(LIB)randomNumberGenerator.o $(LIB)reconLib.o
	$(CC) $(LFLAGS) eventListGenerator.o $(LIB)analysis.o $(LIB)detectorParameters.o \
	$(LIB)parsing.o $(LIB)randomNumberGenerator.o $(LIB)reconLib.o \
	$(GSLFLAGS) -o eventListGenerator

getTriggerEvents: getTriggerEvents.o $(LIB)analysis.o  $(LIB)detectorParameters.o \
  $(LIB)parsing.o  $(LIB)reconLib.o
	$(CC) $(LFLAGS) getTriggerEvents.o $(LIB)analysis.o $(LIB)detectorParameters.o \
	$(LIB)parsing.o $(LIB)reconLib.o \
	$(GSLFLAGS) -o getTriggerEvents
	
recordTriggerEvents: recordTriggerEvents.o $(LIB)analysis.o  $(LIB)detectorParameters.o \
  $(LIB)parsing.o  $(LIB)reconLib.o
	$(CC) $(LFLAGS) recordTriggerEvents.o $(LIB)analysis.o $(LIB)detectorParameters.o \
	$(LIB)parsing.o $(LIB)reconLib.o \
	$(GSLFLAGS) -o recordTriggerEvents
	
eventDisplay: eventDisplay.o $(LIB)parsing.o $(LIB)reconLib.o
	$(CC) $(LFLAGS) eventDisplay.o $(LIB)parsing.o $(LIB)reconLib.o $(GSLFLAGS) -o eventDisplay

makeIterationsFile: makeIterationsFile.o $(LIB)parsing.o
	$(CC) $(LFLAGS) makeIterationsFile.o $(LIB)parsing.o $(GSLFLAGS) -o makeIterationsFile

readGeantData.o: readGeantData.cc
	$(CC) $(CFLAGS) $(INCPATH) -c readGeantData.cc
	
IBD_Generator.o: IBD_Generator.cc
	$(CC) $(CFLAGS) $(INCPATH) -c IBD_Generator.cc

scriptGenerator.o: scriptGenerator.cc
	$(CC) $(CFLAGS) $(INCPATH) -c scriptGenerator.cc
	
addToEventFile.o: addToEventFile.cc $(IN)detectorParameters.hh
	$(CC) $(CFLAGS) $(INCPATH) -c addToEventFile.cc

removeFromEventFile.o: removeFromEventFile.cc $(IN)detectorParameters.hh
	$(CC) $(CFLAGS) $(INCPATH) -c removeFromEventFile.cc

combineIBD_Events.o: combineIBD_Events.cc
	$(CC) $(CFLAGS) $(INCPATH) -c combineIBD_Events.cc

splitIbdForReconstruction.o: splitIbdForReconstruction.cc
	$(CC) $(CFLAGS) $(INCPATH) -c splitIbdForReconstruction.cc
	
basisGenerator.o: basisGenerator.cc
	$(CC) $(CFLAGS) $(INCPATH) -c basisGenerator.cc
	
addBasis.o: addBasis.cc
	$(CC) $(CFLAGS) $(INCPATH) -c addBasis.cc
	
reconstruct.o: reconstruct.cc
	$(CC) $(CFLAGS) $(INCPATH) -c reconstruct.cc

reconstructTeflon.o: reconstructTeflon.cc
	$(CC) $(CFLAGS) $(INCPATH) -c reconstructTeflon.cc

eventListGenerator.o: eventListGenerator.cc
	$(CC) $(CFLAGS) $(INCPATH) -c eventListGenerator.cc
	
getTriggerEvents.o: getTriggerEvents.cc
	$(CC) $(CFLAGS) $(INCPATH) -c getTriggerEvents.cc
	
recordTriggerEvents.o: recordTriggerEvents.cc
	$(CC) $(CFLAGS) $(INCPATH) -c recordTriggerEvents.cc
	
eventDisplay.o: eventDisplay.cc
	$(CC) $(CFLAGS) $(INCPATH) -c eventDisplay.cc
	
makeIterationsFile.o: makeIterationsFile.cc
	$(CC) $(CFLAGS) $(INCPATH) -c makeIterationsFile.cc	
	
# Library object files
#
analysis.o: $(IN)analysis.hh $(SRC)analysis.cc
	$(CC) $(CFLAGS) $(INCPATH) -c $(SRC)analysis.cc -o $(LIB)analysis.o

detectorParameters.o: $(IN)detectorParameters.hh $(SRC)detectorParameters.cc
	$(CC) $(CFLAGS) $(INCPATH) -c $(SRC)detectorParameters.cc -o $(LIB)detectorParameters.o

parsing.o: $(IN)parsing.hh $(SRC)parsing.cc
	$(CC) $(CFLAGS) $(INCPATH) -c $(SRC)parsing.cc -o $(LIB)parsing.o

randomNumberGenerator.o: $(IN)randomNumberGenerator.hh $(SRC)randomNumberGenerator.cc
	$(CC) $(CFLAGS) $(INCPATH) -c $(SRC)randomNumberGenerator.cc -o $(LIB)randomNumberGenerator.o
	
reconLib.o: $(IN)detectorParameters.hh $(IN)reconLib.hh $(SRC)reconLib.cc
	$(CC) $(CFLAGS) $(INCPATH) -c $(SRC)reconLib.cc -o $(LIB)reconLib.o
#	$(CC) $(CFLAGS) -O3 $(INCPATH) -c $(SRC)reconLib.cc -o $(LIB)reconLib.o
	
# Housekeeping
#
cleanLib:
	rm $(LIB)*.o
	
clean:
	rm *.o

cleanAll:
	rm $(LIB)*.o
	rm *.o

cleanEx:
	rm addBasis addToEventFile basisGenerator IBD_Generator readGeantData scriptGenerator trans
	
backup: 
	tar cfv nuLatFullReconstruction.tar $(IN)*.hh $(SRC)*.cc *.cc *.f README $(INPUT)*.*
	mv nuLatFullReconstruction.tar $(BACKUP)nuLatFullReconstruction_$(shell date +'%Y%m%d').tar

tar:
	tar cfv nuLatFullReconstruction_$(shell date +'%Y%m%d').tar $(IN)*.hh $(SRC)*.cc *.cc *.f README $(INPUT)*.*
