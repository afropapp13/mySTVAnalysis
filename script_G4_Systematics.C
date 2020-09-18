{

	vector<TString> WhichSampleArray;
	vector<int> Universes;
	int NG4Universes = 100;
	
	// -----------------------------------------------------------------------------------------
	
	// Run 1 & 3
	
	WhichSampleArray.push_back("_reinteractions_piminus_Geant4"); Universes.push_back(NG4Universes);
//	WhichSampleArray.push_back("_reinteractions_piplus_Geant4"); Universes.push_back(NG4Universes);
//	WhichSampleArray.push_back("_reinteractions_proton_Geant4"); Universes.push_back(NG4Universes);		

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L Efficiencies.cpp++");

	gROOT->ProcessLine(".L XSection_Extraction.cpp++");


	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {
	
		for (int k = 0; k < Universes[i]; k++) {
		//for (int k = 10; k < 35; k++) {
		//for (int k = 50; k < Universes[i]; k++) {		

			gROOT->ProcessLine("Efficiencies(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

			gROOT->ProcessLine("XSection_Extraction(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

		}

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
