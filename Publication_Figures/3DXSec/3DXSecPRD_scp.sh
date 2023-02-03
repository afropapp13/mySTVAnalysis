export OutPutDir=/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_3DXSec_ECal_PRD/Figures/
export InPutDir=/uboone/data/users/apapadop/FlatTTreePlots/ThreeDXSec/
export ANDir=/home/afroditi/Dropbox/PhD/myCode/21th_assignment_CalibratedProducts/CodeRootFiles/uboonecode_v08/mySTVAnalysis/myPlots/pdf/v08_00_00_52/BeamOn9

scp apapadop@uboonegpvm07.fnal.gov:"${InPutDir}"*.pdf "${OutPutDir}"
#scp apapadop@uboonegpvm07.fnal.gov:"${InPutDir}"*.pdf "${ANDir}"
