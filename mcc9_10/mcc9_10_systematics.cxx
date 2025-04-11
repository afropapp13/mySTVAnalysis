#include <iostream>
#include <vector>

void mcc9_10_systematics() {

	// -----------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../../myClasses/Util.C++");

	gROOT->ProcessLine(".L mcc9_10_covariances.cxx++");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine("mcc9_10_covariances(\"Stat\",\"mcc9_10_Overlay9\",\"mcc9_10_BeamOn9\",\"mcc9_10_ExtBNB9\",\"mcc9_10_OverlayDirt9\")");

	gROOT->ProcessLine("mcc9_10_covariances(\"POT\")");

	gROOT->ProcessLine("mcc9_10_covariances(\"NTarget\")");

	gROOT->ProcessLine("mcc9_10_covariances(\"XSec\")");

	gROOT->ProcessLine("mcc9_10_covariances(\"G4\")");

	gROOT->ProcessLine("mcc9_10_covariances(\"Flux\")");

	gROOT->ProcessLine("mcc9_10_covariances(\"Dirt\")");

	/*gROOT->ProcessLine("mcc9_10_covariances(\"LY\")");

	gROOT->ProcessLine("mcc9_10_covariances(\"TPC\")");

	gROOT->ProcessLine("mcc9_10_covariances(\"SCERecomb2\")");*/

	gROOT->ProcessLine("mcc9_10_covariances(\"MC_Stat\")");

	/*gROOT->ProcessLine("mcc9_10_covariances(\"NuWro\")");*/	

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
