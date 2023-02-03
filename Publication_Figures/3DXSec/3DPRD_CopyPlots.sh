export TopoBreakDown=/home/afroditi/Dropbox/PhD/myCode/21th_assignment_CalibratedProducts/CodeRootFiles/uboonecode_v08/myEvents/myPlots/pdf/1D/v08_00_00_52/_NoCuts_PID_NuScore/TopologicalBreakDown/
export InteBreakDown=/home/afroditi/Dropbox/PhD/myCode/21th_assignment_CalibratedProducts/CodeRootFiles/uboonecode_v08/myEvents/myPlots/pdf/1D/v08_00_00_52/_NoCuts_PID_NuScore/InteractionBreakDown/

export BeamOnXSec=/home/afroditi/Dropbox/PhD/myCode/21th_assignment_CalibratedProducts/CodeRootFiles/uboonecode_v08/mySTVAnalysis/myPlots/pdf/v08_00_00_52/BeamOn9/
export OverlayXSec=/home/afroditi/Dropbox/PhD/myCode/21th_assignment_CalibratedProducts/CodeRootFiles/uboonecode_v08/mySTVAnalysis/myPlots/pdf/v08_00_00_52/Overlay9/

export OutPutDir=/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_3DXSec_ECal_PRD/Figures/

declare -a arrPlots=(
"SerialECal_DeltaPTDeltaAlphaTPlot"
"SerialECal_DeltaPtxDeltaPtyPlot"					
"SerialECal_MuonCosThetaMuonMomentumPlot"
"SerialECal_ProtonCosThetaProtonMomentumPlot"
"ECal_DeltaPT_0_00To0_20_DeltaAlphaT_0_00To45_00Plot"
"ECal_DeltaPT_0_00To0_20_DeltaAlphaT_45_00To90_00Plot"
"ECal_DeltaPT_0_00To0_20_DeltaAlphaT_90_00To135_00Plot"	
"ECal_DeltaPT_0_00To0_20_DeltaAlphaT_135_00To180_00Plot"
"ECal_DeltaPT_0_20To0_40_DeltaAlphaT_0_00To45_00Plot"
"ECal_DeltaPT_0_20To0_40_DeltaAlphaT_45_00To90_00Plot"
"ECal_DeltaPT_0_20To0_40_DeltaAlphaT_90_00To135_00Plot"	
"ECal_DeltaPT_0_20To0_40_DeltaAlphaT_135_00To180_00Plot"
"ECal_DeltaPT_0_40To1_00_DeltaAlphaT_0_00To45_00Plot"
"ECal_DeltaPT_0_40To1_00_DeltaAlphaT_45_00To90_00Plot"
"ECal_DeltaPT_0_40To1_00_DeltaAlphaT_90_00To135_00Plot"	
"ECal_DeltaPT_0_40To1_00_DeltaAlphaT_135_00To180_00Plot"
"ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_75ToMinus0_15Plot"		
"ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_15To0_15Plot"
"ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_0_15To0_45Plot"	
"ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_75ToMinus0_15Plot"		
"ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_15To0_15Plot"
"ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_0_15To0_45Plot"
"ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_75ToMinus0_15Plot"		
"ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_15To0_15Plot"
"ECal_DeltaPtx_0_15To0_55_DeltaPty_0_15To0_45Plot"							
"ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_10To0_40Plot"
"ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_10To0_40Plot"
"ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_10To0_40Plot"
"ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_10To0_40Plot"
"ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_40To0_60Plot"
"ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_40To0_60Plot"
"ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_40To0_60Plot"
"ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_40To0_60Plot"
"ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_60To1_20Plot"
"ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_60To1_20Plot"
"ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_60To1_20Plot"
"ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_60To1_20Plot"
"ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_30To0_50Plot"
"ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_30To0_50Plot"
"ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_30To0_50Plot"
"ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_30To0_50Plot"
"ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_50To0_70Plot"
"ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_50To0_70Plot"
"ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_50To0_70Plot"
"ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_50To0_70Plot"
"ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_70To1_00Plot"
"ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_70To1_00Plot"
"ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_70To1_00Plot"
"ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_70To1_00Plot"
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
		
		cp ${OverlayXSec}/Smear_WienerSVD_XSections_${plot}_Combined_v08_00_00_52.pdf	"${OutPutDir}"			
		
		cp ${OverlayXSec}/ClosureTest_WienerSVD_XSections_${plot}_Combined_v08_00_00_52.pdf	"${OutPutDir}"	
		
		cp ${OverlayXSec}/MCERSyst_${plot}Overlay9_Combined.pdf	"${OutPutDir}"	
		
		cp ${OverlayXSec}/Overlay9NuWroWienerSVD_XSections_${plot}_Combined_v08_00_00_52.pdf	"${OutPutDir}"	
		
		cp ${OverlayXSec}/NoTuneOverlay9WienerSVD_XSections_${plot}_Combined_v08_00_00_52.pdf	"${OutPutDir}"	
		
		cp ${OverlayXSec}/TwiceMECOverlay9WienerSVD_XSections_${plot}_Combined_v08_00_00_52.pdf	"${OutPutDir}"						
				
	#fi
	
	##############################################################################

done

declare -a slicesPlots=(
"ECal_DeltaPTDeltaAlphaTPlot_0"
"ECal_DeltaPTDeltaAlphaTPlot_1"
"ECal_DeltaPTDeltaAlphaTPlot_2"
"ECal_DeltaPTDeltaAlphaTPlot_3"
"ECal_DeltaPTDeltaAlphaTPlot_4"
"ECal_DeltaPTDeltaAlphaTPlot_5"
"ECal_DeltaPTDeltaAlphaTPlot_6"
"ECal_DeltaPTDeltaAlphaTPlot_7"
"ECal_DeltaPTDeltaAlphaTPlot_8"
"ECal_DeltaPTDeltaAlphaTPlot_9"
"ECal_DeltaPTDeltaAlphaTPlot_10"
"ECal_DeltaPTDeltaAlphaTPlot_11"
"ECal_DeltaPtxDeltaPtyPlot_0"
"ECal_DeltaPtxDeltaPtyPlot_1"
"ECal_DeltaPtxDeltaPtyPlot_2"
"ECal_DeltaPtxDeltaPtyPlot_3"
"ECal_DeltaPtxDeltaPtyPlot_4"
"ECal_DeltaPtxDeltaPtyPlot_5"
"ECal_DeltaPtxDeltaPtyPlot_6"
"ECal_DeltaPtxDeltaPtyPlot_7"
"ECal_DeltaPtxDeltaPtyPlot_8"
"ECal_MuonCosThetaMuonMomentumPlot_0"
"ECal_MuonCosThetaMuonMomentumPlot_1"
"ECal_MuonCosThetaMuonMomentumPlot_2"
"ECal_MuonCosThetaMuonMomentumPlot_3"
"ECal_MuonCosThetaMuonMomentumPlot_4"
"ECal_MuonCosThetaMuonMomentumPlot_5"
"ECal_MuonCosThetaMuonMomentumPlot_6"
"ECal_MuonCosThetaMuonMomentumPlot_7"
"ECal_MuonCosThetaMuonMomentumPlot_8"
"ECal_MuonCosThetaMuonMomentumPlot_9"
"ECal_MuonCosThetaMuonMomentumPlot_10"
"ECal_MuonCosThetaMuonMomentumPlot_11"
"ECal_ProtonCosThetaProtonMomentumPlot_0"
"ECal_ProtonCosThetaProtonMomentumPlot_1"
"ECal_ProtonCosThetaProtonMomentumPlot_2"
"ECal_ProtonCosThetaProtonMomentumPlot_3"
"ECal_ProtonCosThetaProtonMomentumPlot_4"
"ECal_ProtonCosThetaProtonMomentumPlot_5"
"ECal_ProtonCosThetaProtonMomentumPlot_6"
"ECal_ProtonCosThetaProtonMomentumPlot_7"
"ECal_ProtonCosThetaProtonMomentumPlot_8"
"ECal_ProtonCosThetaProtonMomentumPlot_9"
"ECal_ProtonCosThetaProtonMomentumPlot_10"
"ECal_ProtonCosThetaProtonMomentumPlot_11"
)

for sliceplot in "${slicesPlots[@]}"
do

		cp ${BeamOnXSec}/InteractionBreakdown/InteractionBreakDown_OverlayGENIE_Serial${sliceplot}_Combined_v08_00_00_52.pdf	"${OutPutDir}"

done
