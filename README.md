# Don't forget to run the effective efficiencies as well

root -l script_CV.C
root -l script_Detector_Systematics.C
root -l script_Flux_Systematics.C
root -l script_G4_Systematics.C
root -l script_Genie_Systematics.C

#################################################################################################################################

# Now WienerSVD

root -l script_WienerSVD_Systematics.C

root -l script_WienerSVD_SmEff_Systematics.C

root -l WienerSVD_Merge_Covariances.cpp

root -l WienerSVD_QuantifyUnc.cpp

root -l script_WienerSVD_XSec.C

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

#################################################################################################################################

# Use the difference between the 2 unfolding techniques as an extra uncertainty

root -l CovarianceMatrices_EEvsSVD.cpp

root -l
.L WienerSVD_Merge_Covariances.cpp
WienerSVD_Merge_Covariances("Overlay9",true)

#################################################################################################################################

# Download the relevant xsecs & plots

# (locally)

./myDownloadScripts/DownloadXSec.sh
./myDownloadScripts/DownloadPlots.sh

cd ../myEvents/
./Syst_DownloadEventRatePlots.sh
cd ../mySTVAnalysis/

#################################################################################################################################

# Overlay BeamOn / MC results for different runs
# The cross sections should be independent of the runs

# (locally)

root -l OverlayGenerators.cpp
root -l OverlayXSecMethods.C

root -l 
.L ../myClasses/Util.C
.x WienerSVD_OverlayGenerators.cpp
.x IntegratedXSecs.cpp
.x WienerSVD_Chi2Covariance.cpp
#.x BinByBinChi2.cpp

#################################################################################################################################

