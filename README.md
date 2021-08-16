# Don't forget to run the effective efficiencies as well

root -b script_CV.C
root -b script_Detector_Systematics.C
root -b script_G4_Systematics.C
root -b script_Genie_Systematics.C
root -b script_Flux_Systematics.C

#################################################################################################################################

# Now WienerSVD for GENIE CV / Nominal MC

root -b script_WienerSVD_Systematics.C
root -b WienerSVD_Merge_Covariances.cpp
#root -b WienerSVD_QuantifyUnc.cpp # most likely to be removed 
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

cd ../myEvents
root -b PeLEE_Create1DPlotsTHStack_SubSpectrum.cpp

#locally
./PeLEE_Syst_DownloadEventRatePlots.sh

cd ../mySTVAnalysis

#################################################################################################################################

# Plotting the final results / total uncertainties / putting everything together

root -b Systematics.cpp

root -b
.L ../myClasses/Util.C
.x WienerSVD_DetectorVars.cpp

#################################################################################################################################

# Use the difference between the 2 unfolding techniques as an extra uncertainty

#root -b CovarianceMatrices_EEvsSVD.cpp

#root -b
#.L WienerSVD_Merge_Covariances.cpp
#WienerSVD_Merge_Covariances("Overlay9",true)

#################################################################################################################################

# Download the relevant xsecs & plots

# (locally)

./myDownloadScripts/DownloadXSec.sh
./myDownloadScripts/DownloadPlots.sh

cd ../myEvents/
./PeLEE_Syst_DownloadEventRatePlots.sh
cd ../mySTVAnalysis/

#################################################################################################################################

# Overlay BeamOn / MC results for different runs
# The cross sections should be independent of the runs

# (locally)

root -b OverlayGenerators.cpp
root -b OverlayXSecMethods.C

root -b 
.L ../myClasses/Util.C
.L WienerSVD_OverlayGenerators.cpp
#GENIE versions
WienerSVD_OverlayGenerators()
# AltGen
WienerSVD_OverlayGenerators(false,true)
#GENIE FSI tweaks
WienerSVD_OverlayGenerators(false,false,true)
#GENIE Flag tweaks
WienerSVD_OverlayGenerators(false,false,false,true)
#GENIE Closure Test
WienerSVD_OverlayGenerators(false,false,false,false,true)
#GENIE Nuclear model tweaks
WienerSVD_OverlayGenerators(false,false,false,false,false,true)
.x IntegratedXSecs.cpp
.x WienerSVD_Chi2Covariance.cpp
#.x BinByBinChi2.cpp

#################################################################################################################################

# Fake data studies: we need the MC stat & xsec uncertainties only

# Fake data with CC1p0pi NuWro as data sample # need to be tested 
root -b 
.L WienerSVD_CovarianceMatrices.cpp++
WienerSVD_CovarianceMatrices("Stat","Overlay9","Overlay9NuWro","ExtBNB9","OverlayDirt9")
WienerSVD_CovarianceMatrices("XSec","Overlay9","Overlay9NuWro","ExtBNB9","OverlayDirt9")
# add the rest here after you verify that things are working for NuWro
# do i need to play the same game with MC_Stat?

root -b
.L script_WienerSVD_SmEff_Systematics.cpp++
script_WienerSVD_SmEff_Systematics("SmEff_XSec","Overlay9","Overlay9NuWro","ExtBNB9","OverlayDirt9")
# add the rest here after you verify that things are working for NuWro

# Merge the NuWro covariances # needs testing
root -b
.L WienerSVD_Merge_Covariances.cpp++
WienerSVD_Merge_Covariances("Overlay9","Overlay9NuWro") 

# extract the xsecs # needs testing
root -b 
.L ../../myClasses/Util.C++
.L ../../myClasses/WienerSVD.C++
.L WienerSVD_XSection_Extraction.cpp++
WienerSVD_XSection_Extraction("",false,"Overlay9NuWro")


#################################################################################################################################

