#include <iostream>
#include <vector>

void mcc9_10_run_cv() {

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Nominal Overlay

	WhichSampleArray.push_back("");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L mcc9_10_efficiency.cxx++");

	gROOT->ProcessLine(".L mcc9_10_migration_matrices.cxx++");

	gROOT->ProcessLine(".L mcc9_10_response_matrices.cxx++");

	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("mcc9_10_efficiency(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("mcc9_10_migration_matrices(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("mcc9_10_response_matrices(\""+WhichSampleArray[i]+"\")");

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
