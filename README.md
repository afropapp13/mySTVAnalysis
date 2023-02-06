# Don't forget to run the effective efficiencies as well

root -b script_CV.C
root -b script_Detector_Systematics.C
root -b script_G4_Systematics.C
root -b script_Genie_Systematics.C
root -b script_Flux_Systematics.C
root -b script_MC_Stat_Systematics.C
root -b script_NuWro.C

#################################################################################################################################

# Now WienerSVD for GENIE CV / Nominal MC

root -b script_WienerSVD_Systematics.C
root -b WienerSVD_Merge_Covariances.cpp
root -b script_WienerSVD_XSec.C

#################################################################################################################################
#################################################################################################################################

# Plotting the detector variation xsecs and storing the relevant systematics with respect to the CV sample
root -b Detector_Systematics_LY.cpp
root -b Detector_Systematics_TPC.cpp
root -b Detector_Systematics_SCERecomb2.cpp

# Plotting the nominal xsecs, adding the 2% POT uncertainty and storing the relevant systematics with respect to the CV sample
root -b POT_Systematics.cpp

# Plotting the nominal xsecs, adding the 1% NTarget uncertainty and storing the relevant systematics with respect to the CV sample
root -b NTarget_Systematics.cpp

# Covariance matrices only with diagonal elements filled with stat uncertainties
root -b Stat_Systematics.cpp

# Storing the GEANT4 uncertainties
root -b G4_Systematics.cpp

# Storing the Genie uncertainties
root -b Genie_Systematics.cpp

# Storing the Flux uncertainties
root -b Flux_Systematics.cpp

#################################################################################################################################

# Plotting the final results / total uncertainties / putting everything together

root -b Systematics.cpp

#################################################################################################################################

