#include <iostream>
#include <vector>

void script_WienerSVD_EventRateSystematics() {

	// -----------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../../myClasses/Util.C++");

	gROOT->ProcessLine(".L WienerSVD_EventRateCovarianceMatrices.cpp++");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine("WienerSVD_EventRateCovarianceMatrices(\"Stat\",\"Overlay9\",\"BeamOn9\",\"ExtBNB9\",\"OverlayDirt9\")");

	gROOT->ProcessLine("WienerSVD_EventRateCovarianceMatrices(\"POT\")");

	gROOT->ProcessLine("WienerSVD_EventRateCovarianceMatrices(\"NTarget\")");

	gROOT->ProcessLine("WienerSVD_EventRateCovarianceMatrices(\"XSec\")");

//	gROOT->ProcessLine("WienerSVD_EventRateCovarianceMatrices(\"DetailedXSec\")");	

	gROOT->ProcessLine("WienerSVD_EventRateCovarianceMatrices(\"G4\")");

	gROOT->ProcessLine("WienerSVD_EventRateCovarianceMatrices(\"Flux\")");

	gROOT->ProcessLine("WienerSVD_EventRateCovarianceMatrices(\"Dirt\")");

	gROOT->ProcessLine("WienerSVD_EventRateCovarianceMatrices(\"LY\")");

	gROOT->ProcessLine("WienerSVD_EventRateCovarianceMatrices(\"TPC\")");

	gROOT->ProcessLine("WienerSVD_EventRateCovarianceMatrices(\"SCERecomb2\")");

	gROOT->ProcessLine("WienerSVD_EventRateCovarianceMatrices(\"MC_Stat\")");

	gROOT->ProcessLine("WienerSVD_EventRateCovarianceMatrices(\"NuWro\")");		

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
