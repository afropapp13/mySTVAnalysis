. ../myClasses/Constants.sh

export OutPutDir=/uboone/data/users/$UserID/mySTVAnalysis/myEfficiencies/$UBCode

declare -a arrRun=("Run1")

# Loop over the run numbers

for RunNumber in "${arrRun[@]}"
do

	##############################################################################

	scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/FileEfficiences_Overlay9_${RunNumber}_${UBCode}.root ./myEfficiencies/$UBCode/

	##############################################################################

done

# End of the loop over the run numbers
