{

	vector<TString> WhichSampleArray;
	vector<int> Universes;
	int NGENIEUniverses = 100;
	
	// -----------------------------------------------------------------------------------------
	
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

	gROOT->ProcessLine(".L efficiency.cxx++");

	gROOT->ProcessLine(".L response_matrices.cxx++");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {
	
		for (int k = 0; k < Universes[i]; k++) {

			gROOT->ProcessLine("efficiency(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

			gROOT->ProcessLine("response_matrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

			//gROOT->ProcessLine("response_matrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\",false,\"NoTune\")");
			//gROOT->ProcessLine("response_matrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\",false,\"TwiceMEC\")");			

		}

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
