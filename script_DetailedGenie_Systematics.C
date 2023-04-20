{

	vector<TString> WhichSampleArray;
	vector<int> Universes;
	int NGENIEUniverses = 2;
	
	// -----------------------------------------------------------------------------------------
	
	WhichSampleArray.push_back("_AGKYpT1pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_AGKYxF1pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_AhtBY_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_BhtBY_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_CV1uBY_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_CV2uBY_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_EtaNCEL_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_FrAbs_N_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_FrAbs_pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_FrCEx_N_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_FrCEx_pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_FrInel_N_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_FrInel_pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_FrPiProd_N_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_FrPiProd_pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_FracDelta_CCMEC_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_FracPN_CCMEC_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_MFP_N_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_MFP_pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_MaCCQE_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_MaCCRES_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_MaNCEL_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_MaNCRES_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_MvCCRES_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_MvNCRES_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvbarnCC1pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvbarnCC2pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvbarnNC1pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvbarnNC2pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvbarpCC1pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvbarpCC2pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvbarpNC1pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvbarpNC2pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvnCC1pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvnCC2pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvnNC1pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvnNC2pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvpCC1pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvpCC2pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvpNC1pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NonRESBGvpNC2pi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NormCCMEC_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_NormNCMEC_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_RDecBR1eta_UBGenie"); Universes.push_back(2);				
	WhichSampleArray.push_back("_RDecBR1gamma_UBGenie"); Universes.push_back(2);

	// Unisims

	WhichSampleArray.push_back("_UnShortAxFFCCQEshape_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_UnShortDecayAngMEC_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_UnShortRPA_CCQE_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_UnShortTheta_Delta2Npi_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_UnShortVecFFCCQEshape_UBGenie"); Universes.push_back(2);
	WhichSampleArray.push_back("_UnShortXSecShape_CCMEC_UBGenie"); Universes.push_back(2);	

	// -----------------------------------------------------------------------------------------

//	gROOT->ProcessLine(".L StandardEfficiencies.cpp++");

//	gROOT->ProcessLine(".L EffectiveEfficiencies.cpp++");

	//gROOT->ProcessLine(".L MigrationMatrices.cpp++");

	gROOT->ProcessLine(".L ResponseMatrices.cpp++");

	//gROOT->ProcessLine(".L CovarianceMatrices.cpp++");

//	gROOT->ProcessLine(".L XSection_Extraction.cpp++");

//	gROOT->ProcessLine(".L FFEfficiencies.cpp++");

//	gROOT->ProcessLine(".L FFXSection_Extraction.cpp++");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {
	
		for (int k = 0; k < Universes[i]; k++) {

//			gROOT->ProcessLine("StandardEfficiencies(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

//			gROOT->ProcessLine("EffectiveEfficiencies(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

//			gROOT->ProcessLine("MigrationMatrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

			gROOT->ProcessLine("ResponseMatrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

//			gROOT->ProcessLine("ResponseMatrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\",false,\"NoTune\")");
//			gROOT->ProcessLine("ResponseMatrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\",false,\"TwiceMEC\")");			

//			gROOT->ProcessLine("CovarianceMatrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

//			gROOT->ProcessLine("FFEfficiencies(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

//			gROOT->ProcessLine("XSection_Extraction(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

//			gROOT->ProcessLine("FFXSection_Extraction(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

		}

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
