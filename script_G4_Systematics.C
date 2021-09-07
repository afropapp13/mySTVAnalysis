{

	vector<TString> WhichSampleArray;
	vector<int> Universes;
	int NG4Universes = 100;
	
	// -----------------------------------------------------------------------------------------

	WhichSampleArray.push_back("_reinteractions"); Universes.push_back(NG4Universes);	
//	WhichSampleArray.push_back("_reinteractions_piminus_Geant4"); Universes.push_back(NG4Universes);
//	WhichSampleArray.push_back("_reinteractions_piplus_Geant4"); Universes.push_back(NG4Universes);
//	WhichSampleArray.push_back("_reinteractions_proton_Geant4"); Universes.push_back(NG4Universes);		

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

			//gROOT->ProcessLine("MigrationMatrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

			gROOT->ProcessLine("ResponseMatrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

			gROOT->ProcessLine("ResponseMatrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\",false,\"NoTune\")");

			//gROOT->ProcessLine("CovarianceMatrices(\""+WhichSampleArray[i]+"\")");

//			gROOT->ProcessLine("FFEfficiencies(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

			gROOT->ProcessLine("XSection_Extraction(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

//			gROOT->ProcessLine("FFXSection_Extraction(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

		}

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
