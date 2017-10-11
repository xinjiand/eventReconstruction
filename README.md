# eventReconstruction
========================================================================================================================
This directory contains the following subdirectories:
  lineTests
  nuLatEventReconstruction
  old
  quenching
  
Broad descriptions for these directories are below. More detailed descriptions can be found in the READMEs located in
the directories themselves.
  
  linesTests:               Contains code to model the NuLat line tests with the orthogonal PMT setup. (See Fig. 8.3,
                            for example, in Zach's dissertation.)
              
  nuLatEventReconstruction: Contains Geant4 data, the Fortran-bases light transport code, and primary event
                            reconstruction code for NuLat. This also includes an event reconstruction code for Teflon
                            FEP lattices as well as simulation and analysis results.
                           
  old:                      Contains the directories attenuationLength, shielding, teflonLatticeReconstruction. These 
                            contain unsed (or unfinished) code. There is no README for this directory.
                            
  quenching:                This directory contains the code to simulate the neutrnio tag and A2 background events in 
                            LENS with quenching. The code only gives the true, but quenched, energy deposits. The full
                            quenching analysis includes running the light transport code and doing event reconstruction,
                            which is located in the directory nuLatEventReconstruction.
