{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------
	
	// Run 1 & 3
	
	WhichSampleArray.push_back("_CV");
	WhichSampleArray.push_back("_LYDown");
	WhichSampleArray.push_back("_LYRayleigh");


	// Only Run 3 for now

	WhichSampleArray.push_back("_LYAttenuation");		
	WhichSampleArray.push_back("_WireModX");
	WhichSampleArray.push_back("_WireModYZ");
	WhichSampleArray.push_back("_WireModThetaYZ");
	WhichSampleArray.push_back("_WireModThetaXZ");
//	WhichSampleArray.push_back("_dEdx"); // Buggy, remove it from the analysis for now
	WhichSampleArray.push_back("_Recombination2");
	WhichSampleArray.push_back("_SCE");		

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L Efficiencies.cpp++");

	gROOT->ProcessLine(".L XSection_Extraction.cpp++");


	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("Efficiencies(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("XSection_Extraction(\""+WhichSampleArray[i]+"\")");

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
