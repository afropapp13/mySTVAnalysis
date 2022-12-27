export TopoBreakDown=/home/afroditi/Dropbox/PhD/myCode/21th_assignment_CalibratedProducts/CodeRootFiles/uboonecode_v08/myEvents/myPlots/pdf/1D/v08_00_00_52/_NoCuts_PID_NuScore/TopologicalBreakDown/
export InteBreakDown=/home/afroditi/Dropbox/PhD/myCode/21th_assignment_CalibratedProducts/CodeRootFiles/uboonecode_v08/myEvents/myPlots/pdf/1D/v08_00_00_52/_NoCuts_PID_NuScore/InteractionBreakDown/

export BeamOnXSec=/home/afroditi/Dropbox/PhD/myCode/21th_assignment_CalibratedProducts/CodeRootFiles/uboonecode_v08/mySTVAnalysis/myPlots/pdf/v08_00_00_52/BeamOn9/
export OverlayXSec=/home/afroditi/Dropbox/PhD/myCode/21th_assignment_CalibratedProducts/CodeRootFiles/uboonecode_v08/mySTVAnalysis/myPlots/pdf/v08_00_00_52/Overlay9/

export OutPutDir=/home/afroditi/Dropbox/Apps/Overleaf/General\ Imbalance\ Variables\ for\ Measuring\ Nuclear\ Effects\ and\ Demonstration\ with\ MicroBooNE\ Data/figures/gen/Afro/

declare -a arrPlots=(
"DeltaPhi3DPlot"
"DeltaPnPlot"
"DeltaPnPerpPlot"
"DeltaPnParPlot"
"DeltaAlpha3DqPlot"
"DeltaAlpha3DMuPlot"
"SerialDeltaPn_DeltaAlpha3DqPlot"
"SerialDeltaPn_DeltaAlpha3DMuPlot"
"SerialDeltaAlpha3Dq_DeltaPnPlot"
"SerialDeltaAlpha3DMu_DeltaPnPlot"
"DeltaAlpha3Dq_DeltaPn_0_00To0_20Plot"
"DeltaAlpha3Dq_DeltaPn_0_20To0_40Plot"
"DeltaAlpha3Dq_DeltaPn_0_40To1_00Plot"
"DeltaAlpha3DMu_DeltaPn_0_00To0_20Plot"
"DeltaAlpha3DMu_DeltaPn_0_20To0_40Plot"
"DeltaAlpha3DMu_DeltaPn_0_40To1_00Plot"
"DeltaPn_DeltaAlpha3Dq_0_00To45_00Plot"
"DeltaPn_DeltaAlpha3Dq_45_00To90_00Plot"
"DeltaPn_DeltaAlpha3Dq_90_00To135_00Plot"
"DeltaPn_DeltaAlpha3Dq_135_00To180_00Plot"
"DeltaPn_DeltaAlpha3DMu_0_00To45_00Plot"
"DeltaPn_DeltaAlpha3DMu_45_00To90_00Plot"
"DeltaPn_DeltaAlpha3DMu_90_00To135_00Plot"
"DeltaPn_DeltaAlpha3DMu_135_00To180_00Plot"
"SerialDeltaPn_DeltaAlpha3DqPlot"
"SerialDeltaPn_DeltaAlpha3DMuPlot"
"SerialDeltaAlpha3Dq_DeltaPnPlot"
"SerialDeltaAlpha3DMu_DeltaPnPlot"
"SerialDeltaPn_DeltaAlpha3DqPlot_Slice_0"
"SerialDeltaPn_DeltaAlpha3DqPlot_Slice_1"
"SerialDeltaPn_DeltaAlpha3DqPlot_Slice_2"
"SerialDeltaPn_DeltaAlpha3DqPlot_Slice_3"
"SerialDeltaAlpha3Dq_DeltaPnPlot_Slice_0"
"SerialDeltaAlpha3Dq_DeltaPnPlot_Slice_1"
"SerialDeltaAlpha3Dq_DeltaPnPlot_Slice_2"
"SerialDeltaPn_DeltaAlpha3DMuPlot_Slice_0"
"SerialDeltaPn_DeltaAlpha3DMuPlot_Slice_1"
"SerialDeltaPn_DeltaAlpha3DMuPlot_Slice_2"
"SerialDeltaPn_DeltaAlpha3DMuPlot_Slice_3"
"SerialDeltaAlpha3DMu_DeltaPnPlot_Slice_0"
"SerialDeltaAlpha3DMu_DeltaPnPlot_Slice_1"
"SerialDeltaAlpha3DMu_DeltaPnPlot_Slice_2"
)

for plot in "${arrPlots[@]}"
do

	##############################################################################
	
	#if [[ ${plot} != *"Slice"* ]]; then

		cp ${OverlayXSec}/WienerSVD_Total_CovarianceMatrices_${plot}Overlay9_Combined_v08_00_00_52.pdf	"${OutPutDir}"
		cp ${OverlayXSec}/WienerSVD_Total_FracCovarianceMatrices_${plot}Overlay9_Combined_v08_00_00_52.pdf	"${OutPutDir}"
		
		cp ${TopoBreakDown}/THStack_BreakDown_Reco${plot}_Combined_v08_00_00_52_NoCuts_PID_NuScore.pdf "${OutPutDir}"/THStack_BreakDown_Reco${plot}_Combined_v08_00_00_52_NoCuts_PID_NuScore_Topo.pdf	
	
		cp ${InteBreakDown}/THStack_BreakDown_Reco${plot}_Combined_v08_00_00_52_NoCuts_PID_NuScore.pdf "${OutPutDir}"/THStack_BreakDown_Reco${plot}_Combined_v08_00_00_52_NoCuts_PID_NuScore_Inte.pdf	
		
		cp ${OverlayXSec}/StandardEff${plot}_Combined_v08_00_00_52.pdf	"${OutPutDir}"	
		
		cp ${OverlayXSec}/ResponseMatrices_${plot}Overlay9_Combined_v08_00_00_52.pdf	"${OutPutDir}"
		
		cp ${OverlayXSec}/ClosureTest_WienerSVD_XSections_${plot}_Combined_v08_00_00_52.pdf	"${OutPutDir}"	
		
		cp ${OverlayXSec}/MCERSyst_${plot}Overlay9_Combined.pdf	"${OutPutDir}"	
		
		cp ${BeamOnXSec}/WienerSVD_Generator_TotalUnc_Data_XSections_${plot}_Combined_v08_00_00_52.pdf	"${OutPutDir}"
		
		cp ${BeamOnXSec}/OtherGenWienerSVD_Generator_TotalUnc_Data_XSections_${plot}_Combined_v08_00_00_52.pdf	"${OutPutDir}"
		
		cp ${OverlayXSec}/Overlay9NuWroWienerSVD_XSections_${plot}_Combined_v08_00_00_52.pdf	"${OutPutDir}"	
		
		cp ${OverlayXSec}/NoTuneOverlay9WienerSVD_XSections_${plot}_Combined_v08_00_00_52.pdf	"${OutPutDir}"	
		
		cp ${OverlayXSec}/TwiceMECOverlay9WienerSVD_XSections_${plot}_Combined_v08_00_00_52.pdf	"${OutPutDir}"						
				
	#fi
	
	
	#if [[ ${plot} == *"Serial"*] && [${plot} == *"Slice"* ]]; then
		
		cp ${BeamOnXSec}/MultiDimWienerSVD_Generator_TotalUnc_Data_2DXSections_${plot}_Combined_v08_00_00_52.pdf	"${OutPutDir}"
		cp ${BeamOnXSec}/OtherGenMultiDimWienerSVD_Generator_TotalUnc_Data_2DXSections_${plot}_Combined_v08_00_00_52.pdf	"${OutPutDir}"			

	#fi
	
	##############################################################################

done
