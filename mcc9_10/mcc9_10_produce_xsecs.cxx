#include <iostream>
#include <vector>

void mcc9_10_produce_xsecs() {

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Nominal Overlay

	WhichSampleArray.push_back("");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../../myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L ../../myClasses/Util.C++");
	gROOT->ProcessLine(".L ../../myClasses/WienerSVD.C++");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L mcc9_10_extract_xsec.cxx++");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("mcc9_10_extract_xsec(\""+WhichSampleArray[i]+"\")");

		// Closure test

		gROOT->ProcessLine("mcc9_10_extract_xsec(\""+WhichSampleArray[i]+"\",true)");

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
