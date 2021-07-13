# Don't forget to run the effective efficiencies as well

root -l script_CV.C
root -l script_Detector_Systematics.C
root -l script_Flux_Systematics.C
root -l script_G4_Systematics.C
root -l script_Genie_Systematics.C

#################################################################################################################################

# Now WienerSVD for GENIE CV / Nominal MC

root -l script_WienerSVD_Systematics.C
root -l script_WienerSVD_SmEff_Systematics.C
root -l WienerSVD_Merge_Covariances.cpp
#root -l WienerSVD_QuantifyUnc.cpp # most likely to be removed 
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

cd ../myEvents
root -l PeLEE_Create1DPlotsTHStack_SubSpectrum.cpp

#locally
./PeLEE_Syst_DownloadEventRatePlots.sh

cd ../mySTVAnalysis

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

# Quantifying the systematics

root -l
.L WienerSVD_Uncertainties.cpp

WienerSVD_Uncertainties("Stat")
WienerSVD_Uncertainties("POT")
WienerSVD_Uncertainties("NTarget")
WienerSVD_Uncertainties("LY")
WienerSVD_Uncertainties("TPC")
WienerSVD_Uncertainties("Flux")
WienerSVD_Uncertainties("G4")
WienerSVD_Uncertainties("XSec")

# We are not done yet, need to handle the multisims as well and MC / Sm/Eff unc

root -l PrintFinalUnc.cpp

####################################################################

# Fake data studies: we need the MC stat & xsec uncertainties only

# Fake data with CC1p0pi NuWro as data sample # need to be tested 
root -l 
.L WienerSVD_CovarianceMatrices.cpp++
WienerSVD_CovarianceMatrices("Stat","Overlay9","Overlay9NuWro","ExtBNB9","OverlayDirt9")
WienerSVD_CovarianceMatrices("XSec","Overlay9","Overlay9NuWro","ExtBNB9","OverlayDirt9")
# add the rest here after you verify that things are working for NuWro
# do i need to play the same game with MC_Stat?

root -l
.L script_WienerSVD_SmEff_Systematics.cpp++
script_WienerSVD_SmEff_Systematics("SmEff_XSec","Overlay9","Overlay9NuWro","ExtBNB9","OverlayDirt9")
# add the rest here after you verify that things are working for NuWro

# Merge the NuWro covariances # needs testing
root -l
.L WienerSVD_Merge_Covariances.cpp++
WienerSVD_Merge_Covariances("Overlay9","Overlay9NuWro") 

# extract the xsecs # needs testing
root -l 
.L ../../myClasses/Util.C++
.L ../../myClasses/WienerSVD.C++
.L WienerSVD_XSection_Extraction.cpp++
WienerSVD_XSection_Extraction("",false,"Overlay9NuWro")


#################################################################################################################################

