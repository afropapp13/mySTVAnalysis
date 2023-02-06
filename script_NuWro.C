#include <iostream>
#include <vector>

void script_NuWro() {

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Nominal Overlay

	WhichSampleArray.push_back("");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ResponseMatrices.cpp++");

	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("ResponseMatrices(\""+WhichSampleArray[i]+"\",false,true)");

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
