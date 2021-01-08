#################################################################################################################################


# Plotting the detector variation xsecs and storing the relevant systematics with respect to the CV sample

root -l Detector_Systematics.cpp

# Plotting the nominal xsecs, adding the 2% POT uncertainty and storing the relevant systematics with respect to the CV sample

root -l POT_Systematics.cpp

# Storing the GEANT4 uncertainties

root -l G4_Systematics.cpp

# Storing the Genie uncertainties

root -l Genie_Systematics.cpp

# Storing the Flux uncertainties

root -l Flux_Systematics.cpp


#################################################################################################################################

# Plotting the final results / total uncertainties / putting everything together

root -l Systematics.cpp

#################################################################################################################################

# Overlay BeamOn / MC results for different runs
# The cross sections should be independent of the runs

root -l OverlayXSec.C

#################################################################################################################################

