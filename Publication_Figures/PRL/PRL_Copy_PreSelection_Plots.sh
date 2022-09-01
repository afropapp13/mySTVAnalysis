export InPutDir=/home/afroditi/Dropbox/PhD/myCode/21th_assignment_CalibratedProducts/CodeRootFiles/uboonecode_v08/myPreSelection/myPlots/v08_00_00_52
export OutPutDir=/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance/Figures

declare -a arrPlots=(
"MuonMomentum2DRangeCanvasNoQC_Combined.pdf"
"MuonMomentum2DRangeCanvas_Combined.pdf"
)

for plot in "${arrPlots[@]}"
do

	##############################################################################

	cp ${InPutDir}/${plot}	${OutPutDir}
	
	##############################################################################

done
