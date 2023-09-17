root -b script_CV.C
root -b script_Detector_Systematics.C
root -b script_G4_Systematics.C
root -b script_Genie_Systematics.C
#root -b script_DetailedGenie_Systematics.C
root -b script_Flux_Systematics.C
root -b script_MC_Stat_Systematics.C
root -b script_NuWro.C

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

## Event rate covariances
#root -b script_WienerSVD_Systematics.C

#################################################################################################################################

# Plotting the final results / total uncertainties / putting everything together

root -b Systematics.cpp

# Detector variation breakdown
root -b
.L ../../myClasses/Util.C
.x WienerSVD_DetectorVars.cpp
# XSec variation breakdown
.x WienerSVD_XSecVars.cpp

#################################################################################################################################

# Fake data studies

# We need the stat & xsec uncertainties only

root -b
.L FakeData_WienerSVD_StatCovarianceMatrices.cpp
FakeData_WienerSVD_StatCovarianceMatrices("Stat","Overlay9","Overlay9NuWro","ExtBNB9","OverlayDirt9")
FakeData_WienerSVD_StatCovarianceMatrices("Stat","Overlay9","NoTuneOverlay9","ExtBNB9","OverlayDirt9")
FakeData_WienerSVD_StatCovarianceMatrices("Stat","Overlay9","TwiceMECOverlay9","ExtBNB9","OverlayDirt9")
FakeData_WienerSVD_StatCovarianceMatrices("Stat","Overlay9","GENIEv2Overlay9","ExtBNB9","OverlayDirt9")

# Merge the alternative MC covariances
root -b
.L WienerSVD_Merge_Covariances.cpp++
WienerSVD_Merge_Covariances("Overlay9","Overlay9NuWro","Overlay9NuWro")
WienerSVD_Merge_Covariances("Overlay9","NoTuneOverlay9","NoTuneOverlay9")
WienerSVD_Merge_Covariances("Overlay9","TwiceMECOverlay9","TwiceMECOverlay9") 
WienerSVD_Merge_Covariances("Overlay9","GENIEv2Overlay9","GENIEv2Overlay9") 

#################################################################################################################################

# Fake data studies with Wiener SVD
root -b script_FakeData.C

#################################################################################################################################

# (locally)
# Overlay BeamOn / MC results for different runs
# The cross sections should be independent of the runs

root -b script_ThreeDimWienerSVD_OverlayGenerators.C

hadd -f myXSec/v08_00_00_52/GenXSec/All_XSecs_Combined_v08_00_00_52.root myXSec/v08_00_00_52/GenXSec/All_XSecs_*_Combined_v08_00_00_52.root

root -b DataGeneratorBreakdown.cxx

root -b WienerXSecUncertainty.cxx

#################################################################################################################################

