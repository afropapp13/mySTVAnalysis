#include <iostream>
#include <vector>

void script_CV() {

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Nominal Overlay

	WhichSampleArray.push_back("");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L StandardEfficiencies.cpp++");

	gROOT->ProcessLine(".L EffectiveEfficiencies.cpp++");

	gROOT->ProcessLine(".L MigrationMatrices.cpp++");

	gROOT->ProcessLine(".L ResponseMatrices.cpp++");

//	gROOT->ProcessLine(".L FFEfficiencies.cpp++");

	gROOT->ProcessLine(".L DataDistributions.cpp++");

	gROOT->ProcessLine(".L XSection_Extraction.cpp++");

//	gROOT->ProcessLine(".L FFXSection_Extraction.cpp++");

	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("StarndardEfficiencies(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("EffectiveEfficiencies(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("MigrationMatrices(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("ResponseMatrices(\""+WhichSampleArray[i]+"\")");

//		gROOT->ProcessLine("FFEfficiencies(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("DataDistributions(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("XSection_Extraction(\""+WhichSampleArray[i]+"\")");

//		gROOT->ProcessLine("FFXSection_Extraction(\""+WhichSampleArray[i]+"\")");

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
