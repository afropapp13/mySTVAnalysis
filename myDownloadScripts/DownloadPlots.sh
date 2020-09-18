export UserID=apapadop

export UBCode=v08_00_00_43

export OutPutDir=/uboone/data/users/$UserID/mySTVAnalysis/myPlots/$UBCode

declare -a arrRun=("Run1")

# Loop over the run numbers

for RunNumber in "${arrRun[@]}"
do

	##############################################################################

	scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/BeamOn9/*.pdf ./myPlots/BeamOn9/$UBCode/
	scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/Overlay9/*.pdf ./myPlots/Overlay9/$UBCode/

	##############################################################################

done

# End of the loop over the run numbers
