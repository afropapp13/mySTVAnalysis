#include <iostream>
#include <vector>

void run_cv() {

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Nominal Overlay

	WhichSampleArray.push_back("");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L efficiency.cxx++");

	gROOT->ProcessLine(".L migration_matrices.cxx++");

	gROOT->ProcessLine(".L response_matrices.cxx++");

	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("efficiency(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("migration_matrices(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("response_matrices(\""+WhichSampleArray[i]+"\")");

		//gROOT->ProcessLine("response_matrices(\""+WhichSampleArray[i]+"\",false,\"NoTune\")");

		//gROOT->ProcessLine("response_matrices(\""+WhichSampleArray[i]+"\",false,\"TwiceMEC\")");				

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
