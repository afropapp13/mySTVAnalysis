#################################################################################################################################

root -l script_WienerSVD_Systematics.C

root -l WienerSVD_Merge_Covariances.cpp

root -l WienerSVD_XSection_Extraction.cpp

#################################################################################################################################
#################################################################################################################################

# Plotting the detector variation xsecs and storing the relevant systematics with respect to the CV sample

root -l Detector_Systematics_LY.cpp

root -l Detector_Systematics_TPC.cpp

# Plotting the nominal xsecs, adding the 2% POT uncertainty and storing the relevant systematics with respect to the CV sample

root -l POT_Systematics.cpp

# Plotting the nominal xsecs, adding the 1% NTarget uncertainty and storing the relevant systematics with respect to the CV sample

root -l NTarget_Systematics.cpp

# Covariance matrices only with diagonal elements filled with stat uncertainties

root -l Stat_Systematics.cpp

# Storing the GEANT4 uncertainties

root -l G4_Systematics.cpp

# Storing the Genie uncertainties

root -l Genie_Systematics.cpp

# Storing the Flux uncertainties

root -l Flux_Systematics.cpp


#################################################################################################################################

# Plotting the final results / total uncertainties / putting everything together

root -l Systematics.cpp

# Download the relevant xsecs

#################################################################################################################################

# Overlay BeamOn / MC results for different runs
# The cross sections should be independent of the runs

root -l OverlayXSec.C

root -l WienerSVD_OverlayXSec.C

root -l IntegratedXSecs.C

#################################################################################################################################

