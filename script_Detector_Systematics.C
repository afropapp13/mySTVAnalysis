{

	vector<TString> WhichSampleArray;

	// -----------------------------------------------------------------------------------------
	
	// 1M events for Run1 & most of Run3 

	WhichSampleArray.push_back("_CV");

	WhichSampleArray.push_back("_LYDown");
	WhichSampleArray.push_back("_LYRayleigh");
	WhichSampleArray.push_back("_LYAttenuation");	

	WhichSampleArray.push_back("_X");
	WhichSampleArray.push_back("_YZ");
	WhichSampleArray.push_back("_ThetaYZ");
	WhichSampleArray.push_back("_ThetaXZ");

	// 500k for SCE / Recombination2 / dEdx 

	WhichSampleArray.push_back("_CV");

	WhichSampleArray.push_back("_Recombination2");
	WhichSampleArray.push_back("_SCE");
//	WhichSampleArray.push_back("_dEdx");		

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

		gROOT->ProcessLine("StandardEfficiencies(\""+WhichSampleArray[i]+"\",true)");

		gROOT->ProcessLine("EffectiveEfficiencies(\""+WhichSampleArray[i]+"\",true)");

		//gROOT->ProcessLine("MigrationMatrices(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("ResponseMatrices(\""+WhichSampleArray[i]+"\",true)");

		//gROOT->ProcessLine("CovarianceMatrices(\""+WhichSampleArray[i]+"\")");

//		gROOT->ProcessLine("FFEfficiencies(\""+WhichSampleArray[i]+"\")");

		gROOT->ProcessLine("XSection_Extraction(\""+WhichSampleArray[i]+"\",-1,true)");

//		gROOT->ProcessLine("FFXSection_Extraction(\""+WhichSampleArray[i]+"\")");

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
