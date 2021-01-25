. ../myClasses/Constants.sh

export OutPutDir=/uboone/data/users/$UserID/mySTVAnalysis/myMigrationMatrices/$UBCode

declare -a arrRun=("Run1" "Run3")

# Loop over the run numbers

for RunNumber in "${arrRun[@]}"
do

	##############################################################################

	scp $UserID@$UBgpvm:$OutPutDir/FileMigrationMatrices_Overlay9_${RunNumber}_${UBCode}.root ./myMigrationMatrices/$UBCode/

	##############################################################################

done

# End of the loop over the run numbers
