root -b script_CV.C
root -b script_Detector_Systematics.C
root -b script_G4_Systematics.C
root -b script_Genie_Systematics.C
root -b script_Flux_Systematics.C
root -b script_MC_Stat_Systematics.C

#################################################################################################################################

root -b script_WienerSVD_Systematics.C
root -b script_WienerSVD_Systematics_NoTune.C
root -b script_WienerSVD_Systematics_TwiceMEC.C

root -b WienerSVD_Merge_Covariances.cpp

root -b 
.L WienerSVD_Merge_Covariances.cpp
WienerSVD_Merge_Covariances("Overlay9", "","NoTune")

root -b
.L WienerSVD_Merge_Covariances.cpp
WienerSVD_Merge_Covariances("Overlay9", "","TwiceMEC")

root -b script_WienerSVD_XSec.C

#################################################################################################################################

# Plotting the detector variation xsecs and storing the relevant systematics with respect to the corresponding CV sample
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

#cd ../myEvents
#root -b PeLEE_Create1DPlotsTHStack_SubSpectrum.cpp

##locally
#./PeLEE_Syst_DownloadEventRatePlots.sh

#cd ../mySTVAnalysis

#################################################################################################################################

# Plotting the final results / total uncertainties / putting everything together

root -b Systematics.cpp

# Detector variation breakdown
root -b
.L ../../myClasses/Util.C
.x WienerSVD_DetectorVars.cpp

#################################################################################################################################

# Use the difference between the 2 unfolding techniques as an extra uncertainty

#root -b CovarianceMatrices_EEvsSVD.cpp

#root -b
#.L WienerSVD_Merge_Covariances.cpp
#WienerSVD_Merge_Covariances("Overlay9",true)

#################################################################################################################################

# (locally).L ../../myClasses/Util.C

#cd ../myEvents/
#./PeLEE_Syst_DownloadEventRatePlots.sh
#cd ../mySTVAnalysis/

#################################################################################################################################

# Fake data studies

# We need the stat & xsec uncertainties only

root -b
.L FakeData_WienerSVD_StatCovarianceMatrices.cpp
FakeData_WienerSVD_StatCovarianceMatrices("Stat","Overlay9","Overlay9NuWro","ExtBNB9","OverlayDirt9")
FakeData_WienerSVD_StatCovarianceMatrices("Stat","Overlay9","NoTuneOverlay9","ExtBNB9","OverlayDirt9")
FakeData_WienerSVD_StatCovarianceMatrices("Stat","Overlay9","TwiceMECOverlay9","ExtBNB9","OverlayDirt9")

# Merge the alternative MC covariances
root -b
.L WienerSVD_Merge_Covariances.cpp++
WienerSVD_Merge_Covariances("Overlay9","Overlay9NuWro","Overlay9NuWro")
WienerSVD_Merge_Covariances("Overlay9","NoTuneOverlay9","NoTuneOverlay9")
WienerSVD_Merge_Covariances("Overlay9","TwiceMECOverlay9","TwiceMECOverlay9") 

#################################################################################################################################

# (locally)
./myDownloadScripts/DownloadXSec.sh
# Unfolding Uncertainty 
cd Playground
root -b ModelIndepedent_XSecMethod.C
#cd ../myXSec/v08_00_00_52/
#scp WienerSVD_UnfoldingUnc_Combined_v08_00_00_52.root apapadop@uboonegpvm05.fnal.gov:/uboone/data/users/apapadop/mySTVAnalysis/myXSec/v08_00_00_52
#cd ../..
cd ..

#################################################################################################################################

# Fake data studies
root -b script_FakeData.C

# Fake data with effective efficiences
root -b FakeData_XSection_Extraction.cpp

# (locally)
./myDownloadScripts/DownloadXSec.sh
cd Playground
root -b NuWro_ModelIndepedent_XSecMethod.C
root -b CompareUnc.C
cd ..

#locally
# run ALL the event generators !!!!!!!!!!!!!!!

#################################################################################################################################

# (locally)
# Overlay BeamOn / MC results for different runs
# The cross sections should be independent of the runs

root -b script_WienerSVD_OverlayGenerators.C
root -b script_TwoDimWienerSVD_OverlayGenerators.C
root -b script_ThreeDimWienerSVD_OverlayGenerators.C

hadd -f myXSec/v08_00_00_52/GenXSec/All_XSecs_Combined_v08_00_00_52.root myXSec/v08_00_00_52/GenXSec/All_XSecs_*_Combined_v08_00_00_52.root

root -b GeneratorBand.cxx
root -b DataGeneratorBreakdown.cxx

./myDownloadScripts/CopyBreakDown.sh

root -b WienerXSecUncertainty.cxx

#.x IntegratedXSecs.cpp

#root -b
#.L ../myClasses/Util.C
#.L WienerSVD_Chi2Covariance.cpp
#WienerSVD_Chi2Covariance("Kine")
#WienerSVD_Chi2Covariance("STV")
#WienerSVD_Chi2Covariance("Long")

#.x BinByBinChi2.cpp

# (locally) Pure generator comparisons
#root -b
#.L CompareGenerators.cpp
#CompareGenerators()
#CompareGenerators(false,true)

#################################################################################################################################

# (locally) Effective efficiency
#root -b OverlayGenerators.cpp
#root -b OverlayXSecMethods.C

# Comparison to WC
#root -b WC_OverlayXSecMethods.C

# Comparison to MCC8
#root -b MCC8_OverlayXSecMethods.C

# Download the plots
./myDownloadScripts/DownloadPlots.sh

#################################################################################################################################

