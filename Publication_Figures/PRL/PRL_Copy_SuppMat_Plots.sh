export InPutDirNoCuts=/home/afroditi/Dropbox/PhD/myCode/21th_assignment_CalibratedProducts/CodeRootFiles/uboonecode_v08/myEvents/myPlots/pdf/1D/v08_00_00_52/
export InPutDirXSec=/home/afroditi/Dropbox/PhD/myCode/21th_assignment_CalibratedProducts/CodeRootFiles/uboonecode_v08/mySTVAnalysis/myPlots/pdf/v08_00_00_52/BeamOn9/
export OutPutDir=/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance/Figures

cp ${InPutDirNoCuts}/PRL_SuppMat_RecoProtonLLRPIDPlot_Combined_v08_00_00_52_NoCuts.pdf	${OutPutDir}

cp ${InPutDirXSec}/OtherGenWienerSVD_Generator_TotalUnc_Data_XSections_Serial{DeltaPT_DeltaAlphaT,DeltaAlphaT_DeltaPT,DeltaPtx_DeltaPty}Plot_Combined_v08_00_00_52.pdf	${OutPutDir}
