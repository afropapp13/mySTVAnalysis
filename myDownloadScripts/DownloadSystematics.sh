. ../myClasses/Constants.sh

export OutPutDir=/uboone/data/users/$UserID/mySTVAnalysis/mySystematics/$UBCode

scp $UserID@$UBgpvm:$OutPutDir/*.root ./mySystematics/$UBCode/
