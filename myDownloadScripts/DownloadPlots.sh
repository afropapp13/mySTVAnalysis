. ../myClasses/Constants.sh

export OutPutDir=/uboone/data/users/$UserID/mySTVAnalysis/myPlots/$UBCode

##############################################################################

scp -r $UserID@${UBgpvm}:$OutPutDir/BeamOn9 ./myPlots/pdf/$UBCode/
scp -r $UserID@${UBgpvm}:$OutPutDir/Overlay9 ./myPlots/pdf/$UBCode/

##############################################################################
