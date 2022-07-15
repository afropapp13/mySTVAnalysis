#include <iostream>
#include <vector>

void script_WienerSVD_Systematics_NoTune() {

	// -----------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../../myClasses/Util.C++");
//	gROOT->ProcessLine(".L ../../myClasses/WienerSVD.C++");

	gROOT->ProcessLine(".L WienerSVD_CovarianceMatrices.cpp++");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"Stat\",\"Overlay9\",\"BeamOn9\",\"ExtBNB9\",\"OverlayDirt9\",\"NoTune\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"POT\",\"Overlay9\",\"BeamOn9\",\"ExtBNB9\",\"OverlayDirt9\",\"NoTune\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"NTarget\",\"Overlay9\",\"BeamOn9\",\"ExtBNB9\",\"OverlayDirt9\",\"NoTune\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"XSec\",\"Overlay9\",\"BeamOn9\",\"ExtBNB9\",\"OverlayDirt9\",\"NoTune\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"G4\",\"Overlay9\",\"BeamOn9\",\"ExtBNB9\",\"OverlayDirt9\",\"NoTune\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"Flux\",\"Overlay9\",\"BeamOn9\",\"ExtBNB9\",\"OverlayDirt9\",\"NoTune\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"Dirt\",\"Overlay9\",\"BeamOn9\",\"ExtBNB9\",\"OverlayDirt9\",\"NoTune\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"LY\",\"Overlay9\",\"BeamOn9\",\"ExtBNB9\",\"OverlayDirt9\",\"NoTune\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"TPC\",\"Overlay9\",\"BeamOn9\",\"ExtBNB9\",\"OverlayDirt9\",\"NoTune\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"SCERecomb2\",\"Overlay9\",\"BeamOn9\",\"ExtBNB9\",\"OverlayDirt9\",\"NoTune\")");

	gROOT->ProcessLine("WienerSVD_CovarianceMatrices(\"MC_Stat\",\"Overlay9\",\"BeamOn9\",\"ExtBNB9\",\"OverlayDirt9\",\"NoTune\")");

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
