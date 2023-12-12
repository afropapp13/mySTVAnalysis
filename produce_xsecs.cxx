#include <iostream>
#include <vector>

void produce_xsecs() {

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Nominal Overlay

	WhichSampleArray.push_back("");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L ../myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L ../myClasses/Util.C++");
	gROOT->ProcessLine(".L ../myClasses/WienerSVD.C++");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L extract_xsec.cxx++");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("extract_xsec(\""+WhichSampleArray[i]+"\")");

		// Closure test

		gROOT->ProcessLine("extract_xsec(\""+WhichSampleArray[i]+"\",true)");

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
