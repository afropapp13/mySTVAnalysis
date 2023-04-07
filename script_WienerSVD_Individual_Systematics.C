#include <iostream>
#include <vector>

void script_WienerSVD_Individual_Systematics() {

	// -----------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../../myClasses/Util.C++");

	gROOT->ProcessLine(".L WienerSVD_IndividualCovarianceMatrices.cpp++");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"All_UBGenie\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"AxFFCCQEshape_UBGenie\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"DecayAngMEC_UBGenie\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NormCCCOH_UBGenie\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NormNCCOH_UBGenie\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"RPA_CCQE_UBGenie\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"ThetaDelta2NRad_UBGenie\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"Theta_Delta2Npi_UBGenie\")");		
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"VecFFCCQEshape_UBGenie\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"XSecShape_CCMEC_UBGenie\")");	

	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"AGKYpT1pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"AGKYxF1pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"AhtBY_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"BhtBY_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"CV1uBY_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"CV2uBY_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"EtaNCEL_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"FrAbs_N_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"FrAbs_pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"FrCEx_N_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"FrCEx_pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"FrInel_N_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"FrInel_pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"FrPiProd_N_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"FrPiProd_pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"FracDelta_CCMEC_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"FracPN_CCMEC_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"MFP_N_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"MFP_pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"MaCCQE_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"MaCCRES_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"MaNCEL_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"MaNCRES_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"MvCCRES_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"MvNCRES_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvbarnCC1pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvbarnCC2pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvbarnNC1pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvbarnNC2pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvbarpCC1pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvbarpCC2pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvbarpNC1pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvbarpNC2pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvnCC1pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvnCC2pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvnNC1pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvnNC2pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvpCC1pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvpCC2pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvpNC1pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NonRESBGvpNC2pi_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NormCCMEC_UBGenie\",\"Run1_DecompXSecUnc\")");		
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"NormNCMEC_UBGenie\",\"Run1_DecompXSecUnc\")");
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"RDecBR1eta_UBGenie\",\"Run1_DecompXSecUnc\")");		
	gROOT->ProcessLine("WienerSVD_IndividualCovarianceMatrices(\"RDecBR1gamma_UBGenie\",\"Run1_DecompXSecUnc\")");			

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
