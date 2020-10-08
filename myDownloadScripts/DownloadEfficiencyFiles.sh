export UserID=apapadop

export UBCode=v08_00_00_43

export OutPutDir=/uboone/data/users/$UserID/mySTVAnalysis/myEfficiencies/$UBCode

declare -a arrRun=("Run1")

# Loop over the run numbers

for RunNumber in "${arrRun[@]}"
do

	##############################################################################

	scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/FileEfficiences_Overlay9_Run1_${UBCode}.root ./myEfficiencies/$UBCode/
	#scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/FileEfficiences_Overlay9_Run2_All_UBGenie_0_v08_00_00_43.root ./myEfficiencies/$UBCode/
	scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/FileEfficiences_Overlay9_Run3_${UBCode}.root ./myEfficiencies/$UBCode/
	#scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/FileEfficiences_Overlay9_Run4_All_UBGenie_0_v08_00_00_43.root ./myEfficiencies/$UBCode/
	#scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/FileEfficiences_Overlay9_Run5_All_UBGenie_0_v08_00_00_43.root ./myEfficiencies/$UBCode/

	##############################################################################

done

# End of the loop over the run numbers
