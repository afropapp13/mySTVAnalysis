. ../myClasses/Constants.sh

export InPutDir=/home/afroditi/Dropbox/PhD/myCode/21th_assignment_CalibratedProducts/CodeRootFiles/uboonecode_v08/mySTVAnalysis/myPlots/pdf/${UBCode}/BeamOn9/InteractionBreakdown
export OutPutDir=/home/afroditi/Dropbox/Apps/Overleaf/XSec_Playground/Figures

declare -a arrPlots=(
"InteractionBreakDown_OverlayGENIE_DeltaPTPlot_Combined_"${UBCode}".pdf"
"InteractionBreakDown_Genie_v3_0_6_NoFSI_DeltaPTPlot_Combined_"${UBCode}".pdf"
"InteractionBreakDown_GiBUU_DeltaPTPlot_Combined_"${UBCode}".pdf"
"InteractionBreakDown_GiBUUNoFSI_DeltaPTPlot_Combined_"${UBCode}".pdf"

"InteractionBreakDown_OverlayGENIE_SerialDeltaPT_DeltaAlphaTPlot_0_Combined_"${UBCode}".pdf"
"InteractionBreakDown_Genie_v3_0_6_NoFSI_SerialDeltaPT_DeltaAlphaTPlot_0_Combined_"${UBCode}".pdf"
"InteractionBreakDown_GiBUU_SerialDeltaPT_DeltaAlphaTPlot_0_Combined_"${UBCode}".pdf"
"InteractionBreakDown_GiBUUNoFSI_SerialDeltaPT_DeltaAlphaTPlot_0_Combined_"${UBCode}".pdf"

"InteractionBreakDown_OverlayGENIE_SerialDeltaPT_DeltaAlphaTPlot_3_Combined_"${UBCode}".pdf"
"InteractionBreakDown_Genie_v3_0_6_NoFSI_SerialDeltaPT_DeltaAlphaTPlot_3_Combined_"${UBCode}".pdf"
"InteractionBreakDown_GiBUU_SerialDeltaPT_DeltaAlphaTPlot_3_Combined_"${UBCode}".pdf"
"InteractionBreakDown_GiBUUNoFSI_SerialDeltaPT_DeltaAlphaTPlot_3_Combined_"${UBCode}".pdf"

"InteractionBreakDown_G21NoFSI_DeltaAlphaTPlot_Combined_"${UBCode}".pdf"
"InteractionBreakDown_SuSav2_DeltaAlphaTPlot_Combined_"${UBCode}".pdf"
"InteractionBreakDown_G21hA_DeltaAlphaTPlot_Combined_"${UBCode}".pdf"
"InteractionBreakDown_G21G4_DeltaAlphaTPlot_Combined_"${UBCode}".pdf"

"InteractionBreakDown_G21NoFSI_SerialDeltaAlphaT_DeltaPTPlot_0_Combined_"${UBCode}".pdf"
"InteractionBreakDown_SuSav2_SerialDeltaAlphaT_DeltaPTPlot_0_Combined_"${UBCode}".pdf"
"InteractionBreakDown_G21hA_SerialDeltaAlphaT_DeltaPTPlot_0_Combined_"${UBCode}".pdf"
"InteractionBreakDown_G21G4_SerialDeltaAlphaT_DeltaPTPlot_0_Combined_"${UBCode}".pdf"

"InteractionBreakDown_G21NoFSI_SerialDeltaAlphaT_DeltaPTPlot_2_Combined_"${UBCode}".pdf"
"InteractionBreakDown_SuSav2_SerialDeltaAlphaT_DeltaPTPlot_2_Combined_"${UBCode}".pdf"
"InteractionBreakDown_G21hA_SerialDeltaAlphaT_DeltaPTPlot_2_Combined_"${UBCode}".pdf"
"InteractionBreakDown_G21G4_SerialDeltaAlphaT_DeltaPTPlot_2_Combined_"${UBCode}".pdf"

)

# Loop over the plots

for plot in "${arrPlots[@]}"
do

	##############################################################################

	cp ${InPutDir}/${plot}	${OutPutDir}
	
	##############################################################################

done

# End of the loop over the plots


