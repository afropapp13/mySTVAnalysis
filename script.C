#include <iostream>
#include <vector>

//#include  "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
//#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
//#include "/home/afroditi/Dropbox/PhD/Secondary_Code/ToString.cpp"
//#include "../myClasses/Constants.h"

void script(){

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Nominal Overlay

	WhichSampleArray.push_back("");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L Create1DPlotsTotal.cpp++");

	gROOT->ProcessLine(".L MigrationMatrices.cpp++");

	gROOT->ProcessLine(".L DataDistributions.cpp++");

	gROOT->ProcessLine(".L XSection_Extraction.cpp++");


	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("Create1DPlotsTotal(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("MigrationMatrices(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("DataDistributions(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("XSection_Extraction(\""+WhichSampleArray[i]+"\")");

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
