export UserID=apapadop

export UBCode=v08_00_00_43

export OutPutDir=/uboone/data/users/$UserID/mySTVAnalysis/mySystematics/$UBCode

declare -a arrRun=("Run1")

# Loop over the run numbers

for RunNumber in "${arrRun[@]}"
do

	##############################################################################

	scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/*.root ./mySystematics/$UBCode/

	##############################################################################

done

# End of the loop over the run numbers
