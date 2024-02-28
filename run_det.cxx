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

	WhichSampleArray.push_back("_CVextra");

	WhichSampleArray.push_back("_Recombination2");
	WhichSampleArray.push_back("_SCE");
//	WhichSampleArray.push_back("_dEdx");		

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L efficiency.cxx++");

	gROOT->ProcessLine(".L response_matrices.cxx++");

	for (int i =0;i < (int)(WhichSampleArray.size()); i++) {

		gROOT->ProcessLine("efficiency(\""+WhichSampleArray[i]+"\",true)");

		gROOT->ProcessLine("response_matrices(\""+WhichSampleArray[i]+"\",true)");

		//gROOT->ProcessLine("response_matrices(\""+WhichSampleArray[i]+"\",true,\"NoTune\")");
		//gROOT->ProcessLine("response_matrices(\""+WhichSampleArray[i]+"\",true,\"TwiceMEC\")");				

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
