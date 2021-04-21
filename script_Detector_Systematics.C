{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------
	
	WhichSampleArray.push_back("_CV");

	WhichSampleArray.push_back("_LYDown");
	WhichSampleArray.push_back("_LYRayleigh");
	WhichSampleArray.push_back("_LYAttenuation");	
	
	WhichSampleArray.push_back("_WireModX");
	WhichSampleArray.push_back("_WireModYZ");
	WhichSampleArray.push_back("_WireModThetaYZ");
	WhichSampleArray.push_back("_WireModThetaXZ");
	WhichSampleArray.push_back("_dEdx");
	WhichSampleArray.push_back("_Recombination2");
	WhichSampleArray.push_back("_SCE");		

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L StandardEfficiencies.cpp++");

	gROOT->ProcessLine(".L EffectiveEfficiencies.cpp++");

	//gROOT->ProcessLine(".L MigrationMatrices.cpp++");

	gROOT->ProcessLine(".L ResponseMatrices.cpp++");

	//gROOT->ProcessLine(".L CovarianceMatrices.cpp++");

//	gROOT->ProcessLine(".L FFEfficiencies.cpp++");

	gROOT->ProcessLine(".L XSection_Extraction.cpp++");

//	gROOT->ProcessLine(".L FFXSection_Extraction.cpp++");

	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("StandardEfficiencies(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("EffectiveEfficiencies(\""+WhichSampleArray[i]+"\")");

		//gROOT->ProcessLine("MigrationMatrices(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("ResponseMatrices(\""+WhichSampleArray[i]+"\")");

		//gROOT->ProcessLine("CovarianceMatrices(\""+WhichSampleArray[i]+"\")");

//		gROOT->ProcessLine("FFEfficiencies(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("XSection_Extraction(\""+WhichSampleArray[i]+"\")");

//		gROOT->ProcessLine("FFXSection_Extraction(\""+WhichSampleArray[i]+"\")");

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
