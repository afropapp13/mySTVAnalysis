#include <iostream>
#include <vector>

void script_WienerSVD_Systematics() {

	// -----------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../../myClasses/Util.C++");

	gROOT->ProcessLine(".L WienerSVD_CovarianceMatrices.cpp++");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"Stat\",\"Overlay9\",\"BeamOn9\",\"ExtBNB9\",\"OverlayDirt9\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"POT\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"NTarget\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"XSec\")");

//	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"DetailedXSec\")");	

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"G4\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"Flux\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"Dirt\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"LY\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"TPC\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"SCERecomb2\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"MC_Stat\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"NuWro\")");		

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
