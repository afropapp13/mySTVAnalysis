#include <iostream>
#include <vector>

void script_WienerSVD_XSec() {

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Nominal Overlay

	WhichSampleArray.push_back("");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../../myClasses/Util.C++");
	gROOT->ProcessLine(".L ../../myClasses/WienerSVD.C++");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L WienerSVD_XSection_Extraction.cpp++");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("WienerSVD_XSection_Extraction(\""+WhichSampleArray[i]+"\")");

		// Closure test

		gROOT->ProcessLine("WienerSVD_XSection_Extraction(\""+WhichSampleArray[i]+"\",true)");		

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
