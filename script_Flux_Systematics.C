{

	vector<TString> WhichSampleArray;
	vector<int> Universes;
	int NFluxUniverses = 100;
	//int NFluxUniverses = 1;	
	
	// -----------------------------------------------------------------------------------------
	
	// Run 1
	
	WhichSampleArray.push_back("_horncurrent_FluxUnisim"); Universes.push_back(NFluxUniverses);
	WhichSampleArray.push_back("_kminus_PrimaryHadronNormalization"); Universes.push_back(NFluxUniverses);
	WhichSampleArray.push_back("_kplus_PrimaryHadronFeynmanScaling"); Universes.push_back(NFluxUniverses);
	WhichSampleArray.push_back("_kzero_PrimaryHadronSanfordWang"); Universes.push_back(NFluxUniverses);
	WhichSampleArray.push_back("_nucleoninexsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
	WhichSampleArray.push_back("_nucleonqexsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
	WhichSampleArray.push_back("_nucleontotxsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
	WhichSampleArray.push_back("_piminus_PrimaryHadronSWCentralSplineVariation"); Universes.push_back(NFluxUniverses);
	WhichSampleArray.push_back("_pioninexsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
	WhichSampleArray.push_back("_pionqexsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
	WhichSampleArray.push_back("_piontotxsec_FluxUnisim"); Universes.push_back(NFluxUniverses);
	WhichSampleArray.push_back("_piplus_PrimaryHadronSWCentralSplineVariation"); Universes.push_back(NFluxUniverses);
	WhichSampleArray.push_back("_expskin_FluxUnisim"); Universes.push_back(10);  Watch out, different models			

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L Create1DPlotsTotal.cpp++");

	gROOT->ProcessLine(".L XSection_Extraction.cpp++");


	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {
	
		//for (int k = 0; k < Universes[i]; k++) {
		for (int k = 0; k < 50; k++) {		
		//for (int k = 50; k < Universes[i]; k++) {		

			gROOT->ProcessLine("Create1DPlotsTotal(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

			gROOT->ProcessLine("XSection_Extraction(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\","+TString(std::to_string(k))+")");

		}

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
