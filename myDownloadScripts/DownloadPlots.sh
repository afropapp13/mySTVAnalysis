. ../myClasses/Constants.sh

export OutPutDir=/uboone/data/users/$UserID/mySTVAnalysis/myPlots/$UBCode

##############################################################################

scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/BeamOn9/*.pdf ./myPlots/pdf/$UBCode/BeamOn9
scp $UserID@uboonegpvm05.fnal.gov:$OutPutDir/Overlay9/*.pdf ./myPlots/pdf/$UBCode/Overlay9

##############################################################################

# End of the loop over the run numbers
