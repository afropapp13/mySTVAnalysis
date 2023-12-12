#include <iostream>
#include <vector>

void systematics() {

	// -----------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../myClasses/Util.C++");

	gROOT->ProcessLine(".L covariances.cxx++");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine("covariances(\"Stat\",\"Overlay9\",\"BeamOn9\",\"ExtBNB9\",\"OverlayDirt9\")");

	gROOT->ProcessLine("covariances(\"POT\")");

	gROOT->ProcessLine("covariances(\"NTarget\")");

	gROOT->ProcessLine("covariances(\"XSec\")");

	gROOT->ProcessLine("covariances(\"G4\")");

	gROOT->ProcessLine("covariances(\"Flux\")");

	gROOT->ProcessLine("covariances(\"Dirt\")");

	gROOT->ProcessLine("covariances(\"LY\")");

	gROOT->ProcessLine("covariances(\"TPC\")");

	gROOT->ProcessLine("covariances(\"SCERecomb2\")");

	gROOT->ProcessLine("covariances(\"MC_Stat\")");

	gROOT->ProcessLine("covariances(\"NuWro\")");		

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
