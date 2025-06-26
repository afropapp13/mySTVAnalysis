#include <iostream>
#include <vector>

void mcc9_10_event_rate_systematics() {

	// -----------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../../../generators/Util.C++");

	gROOT->ProcessLine(".L mcc9_10_event_rate_covariances.cxx++");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine("mcc9_10_event_rate_covariances(\"Stat\",\"mcc9_10_Overlay9\",\"mcc9_10_BeamOn9\",\"mcc9_10_ExtBNB9\",\"mcc9_10_OverlayDirt9\")");

	gROOT->ProcessLine("mcc9_10_event_rate_covariances(\"POT\")");

	gROOT->ProcessLine("mcc9_10_event_rate_covariances(\"NTarget\")");

	gROOT->ProcessLine("mcc9_10_event_rate_covariances(\"XSec\")");

	gROOT->ProcessLine("mcc9_10_event_rate_covariances(\"G4\")");

	gROOT->ProcessLine("mcc9_10_event_rate_covariances(\"Flux\")");

	gROOT->ProcessLine("mcc9_10_event_rate_covariances(\"Dirt\")");

	/*gROOT->ProcessLine("mcc9_10_event_rate_covariances(\"LY\")");

	gROOT->ProcessLine("mcc9_10_event_rate_covariances(\"TPC\")");

	gROOT->ProcessLine("mcc9_10_event_rate_covariances(\"SCERecomb2\")");*/

	gROOT->ProcessLine("mcc9_10_event_rate_covariances(\"MC_Stat\")");

	gROOT->ProcessLine("mcc9_10_event_rate_covariances(\"test_det\")");	

	/*gROOT->ProcessLine("mcc9_10_event_rate_covariances(\"NuWro\")");*/

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
