#include <iostream>
#include <vector>

void run_nuwro() {

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Nominal Overlay

	WhichSampleArray.push_back("");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L response_matrices.cxx++");

	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("response_matrices(\""+WhichSampleArray[i]+"\",false,\"\",true)");

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
