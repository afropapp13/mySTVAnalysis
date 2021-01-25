. ../myClasses/Constants.sh

export OutPutDir=/uboone/data/users/$UserID/mySTVAnalysis/myXSec/$UBCode

declare -a arrRun=("Run1" "Run3")

# Loop over the run numbers

for RunNumber in "${arrRun[@]}"
do

	##############################################################################

	scp $UserID@$UBgpvm:$OutPutDir/ExtractedXSec_Overlay9_${RunNumber}_${UBCode}.root ./myXSec/$UBCode/

	##############################################################################

done

# End of the loop over the run numbers
