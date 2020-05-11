{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------

	// Run 1 Systematics

	WhichSampleArray.push_back("_CV");
	WhichSampleArray.push_back("_X");
	WhichSampleArray.push_back("_YZ");
	WhichSampleArray.push_back("_LY");
	WhichSampleArray.push_back("_LYRayleigh");

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L Create1DPlotsTotal.cpp++");

	gROOT->ProcessLine(".L XSection_Extraction.cpp++");


	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("Create1DPlotsTotal(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("XSection_Extraction(\""+WhichSampleArray[i]+"\")");

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
