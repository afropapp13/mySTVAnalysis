#include <iostream>
#include <vector>

void script_WienerSVD_SmEff_Systematics() {

	// -----------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../../myClasses/Util.C++");
	gROOT->ProcessLine(".L ../../myClasses/WienerSVD.C++");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L WienerSVD_SmearEff_CovarianceMatrices.cpp++");

	// -----------------------------------------------------------------------------------------

//	gROOT->ProcessLine("WienerSVD_SmearEff_CovarianceMatrices(\"SmEff_Stat\")");

	gROOT->ProcessLine("WienerSVD_SmearEff_CovarianceMatrices(\"SmEff_POT\")");

	gROOT->ProcessLine("WienerSVD_SmearEff_CovarianceMatrices(\"SmEff_NTarget\")");

	gROOT->ProcessLine("WienerSVD_SmearEff_CovarianceMatrices(\"SmEff_LY\")");

	gROOT->ProcessLine("WienerSVD_SmearEff_CovarianceMatrices(\"SmEff_TPC\")");

	gROOT->ProcessLine("WienerSVD_SmearEff_CovarianceMatrices(\"SmEff_XSec\")");

	gROOT->ProcessLine("WienerSVD_SmearEff_CovarianceMatrices(\"SmEff_G4\")");

	gROOT->ProcessLine("WienerSVD_SmearEff_CovarianceMatrices(\"SmEff_Flux\")");

//	gROOT->ProcessLine("WienerSVD_SmearEff_CovarianceMatrices(\"SmEff_Dirt\")");

	// -----------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
