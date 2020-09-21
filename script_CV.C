#include <iostream>
#include <vector>

void script_CV(){

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Nominal Overlay

	WhichSampleArray.push_back("");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L Efficiencies.cpp++");

	gROOT->ProcessLine(".L MigrationMatrices.cpp++");

	gROOT->ProcessLine(".L DataDistributions.cpp++");

	gROOT->ProcessLine(".L XSection_Extraction.cpp++");


	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("Efficiencies(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("MigrationMatrices(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("DataDistributions(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("XSection_Extraction(\""+WhichSampleArray[i]+"\")");

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
