export UserID=apapadop

export UBCode=v08_00_00_52

export OutPutDir=/uboone/data/users/$UserID/mySTVAnalysis/myXSec/$UBCode

declare -a arrRun=("Run1")

# Loop over the run numbers

for RunNumber in "${arrRun[@]}"
do

	##############################################################################

	scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/ExtractedXSec_Overlay9_Run1_${UBCode}.root ./myXSec/$UBCode/
	#scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/ExtractedXSec_Overlay9_Run2_${UBCode}.root ./myXSec/$UBCode/
	scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/ExtractedXSec_Overlay9_Run3_${UBCode}.root ./myXSec/$UBCode/
	#scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/ExtractedXSec_Overlay9_Run4_${UBCode}.root ./myXSec/$UBCode/
	#scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/ExtractedXSec_Overlay9_Run5_${UBCode}.root ./myXSec/$UBCode/

	##############################################################################

done

# End of the loop over the run numbers
