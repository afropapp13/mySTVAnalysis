. ../myClasses/Constants.sh

export OutPutDir=/uboone/data/users/$UserID/mySTVAnalysis/myPlots/$UBCode

##############################################################################

scp $UserID@${UBgpvm}:$OutPutDir/BeamOn9/*.pdf ./myPlots/pdf/$UBCode/BeamOn9
scp $UserID@${UBgpvm}:$OutPutDir/Overlay9/*.pdf ./myPlots/pdf/$UBCode/Overlay9

##############################################################################
