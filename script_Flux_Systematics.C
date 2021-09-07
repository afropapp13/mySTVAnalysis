{

	vector<TString> WhichSampleArray;
	vector<int> Universes;
	int NFluxUniverses = 100;
	
	// -----------------------------------------------------------------------------------------

	WhichSampleArray.push_back("_fluxes"); Universes.push_back(NFluxUniverses);	
//	WhichSampleArray.push_back("_horncurrent_FluxUnisim"); Universes.push_back(NFluxUniverses);
//	WhichSampleArray.push_back("_kminus_PrimaryHadronNormalization"); Universes.push_back(NFluxUniverses);
//	WhichSampleArray.push_back("_kplus_PrimaryHadronFeynmanScaling"); Universes.push_back(NFluxUniverses);
//	WhichSampleArray.push_back("_kzero_PrimaryHadronSanfordWang"); Universes.push_back(NFluxUniverses);
//	WhichSampleArray.push_back("_nucleoninexsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
//	WhichSampleArray.push_back("_nucleonqexsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
//	WhichSampleArray.push_back("_nucleontotxsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
//	WhichSampleArray.push_back("_piminus_PrimaryHadronSWCentralSplineVariation"); Universes.push_back(NFluxUniverses);
//	WhichSampleArray.push_back("_pioninexsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
//	WhichSampleArray.push_back("_pionqexsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
//	WhichSampleArray.push_back("_piontotxsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
//	WhichSampleArray.push_back("_piplus_PrimaryHadronSWCentralSplineVariation"); Universes.push_back(NFluxUniverses);
//	WhichSampleArray.push_back("_expskin_FluxUnisim"); Universes.push_back(10); // Watch out, different models			

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

			//gROOT->ProcessLine("CovarianceMatrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

//			gROOT->ProcessLine("FFEfficiencies(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

			gROOT->ProcessLine("XSection_Extraction(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\","+TString(std::to_string(k))+")");

//			gROOT->ProcessLine("FFXSection_Extraction(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\","+TString(std::to_string(k))+")");

		}

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
