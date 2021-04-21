{

	vector<TString> WhichSampleArray;
	vector<int> Universes;
	int NGENIEUniverses = 100;
	
	// -----------------------------------------------------------------------------------------
	
	// Run 1 for now
	
	WhichSampleArray.push_back("_AxFFCCQEshape_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_DecayAngMEC_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NormCCCOH_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NormNCCOH_UBGenie"); Universes.push_back(2);
//	WhichSampleArray.push_back("_RPA_CCQE_Reduced_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_RPA_CCQE_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_ThetaDelta2NRad_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_Theta_Delta2Npi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_VecFFCCQEshape_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_XSecShape_CCMEC_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_All_UBGenie"); Universes.push_back(NGENIEUniverses);		

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L StandardEfficiencies.cpp++");

	gROOT->ProcessLine(".L EffectiveEfficiencies.cpp++");

	//gROOT->ProcessLine(".L MigrationMatrices.cpp++");

	gROOT->ProcessLine(".L ResponseMatrices.cpp++");

	//gROOT->ProcessLine(".L CovarianceMatrices.cpp++");

	gROOT->ProcessLine(".L XSection_Extraction.cpp++");

//	gROOT->ProcessLine(".L FFEfficiencies.cpp++");

//	gROOT->ProcessLine(".L FFXSection_Extraction.cpp++");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {
	
		for (int k = 0; k < Universes[i]; k++) {

			gROOT->ProcessLine("StandardEfficiencies(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

			gROOT->ProcessLine("EffectiveEfficiencies(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

//			gROOT->ProcessLine("MigrationMatrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

			gROOT->ProcessLine("ResponseMatrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

//			gROOT->ProcessLine("CovarianceMatrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

//			gROOT->ProcessLine("FFEfficiencies(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

			gROOT->ProcessLine("XSection_Extraction(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

//			gROOT->ProcessLine("FFXSection_Extraction(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

		}

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
