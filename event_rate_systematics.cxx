#include <iostream>
#include <vector>

void event_rate_systematics() {

	// -----------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../myClasses/Util.C++");

	gROOT->ProcessLine(".L event_rate_covariances.cxx++");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine("event_rate_covariances(\"Stat\",\"Overlay9\",\"BeamOn9\",\"ExtBNB9\",\"OverlayDirt9\")");

	gROOT->ProcessLine("event_rate_covariances(\"POT\")");

	gROOT->ProcessLine("event_rate_covariances(\"NTarget\")");

	gROOT->ProcessLine("event_rate_covariances(\"XSec\")");

	gROOT->ProcessLine("event_rate_covariances(\"G4\")");

	gROOT->ProcessLine("event_rate_covariances(\"Flux\")");

	gROOT->ProcessLine("event_rate_covariances(\"Dirt\")");

	gROOT->ProcessLine("event_rate_covariances(\"LY\")");

	gROOT->ProcessLine("event_rate_covariances(\"TPC\")");

	gROOT->ProcessLine("event_rate_covariances(\"SCERecomb2\")");

	gROOT->ProcessLine("event_rate_covariances(\"MC_Stat\")");

	gROOT->ProcessLine("event_rate_covariances(\"NuWro\")");		

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
