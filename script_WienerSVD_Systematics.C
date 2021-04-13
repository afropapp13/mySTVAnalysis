#include <iostream>
#include <vector>

void script_WienerSVD_Systematics() {

	// -----------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L WienerSVD_CovarianceMatrices.cpp++");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"NTarget\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"POT\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"Stat\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"LY\",\"Overlay9_CV\")");

	// -----------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
