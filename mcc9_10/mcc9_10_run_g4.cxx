{

	vector<TString> WhichSampleArray;
	vector<int> Universes;
	int NG4Universes = 100;
	
	// -----------------------------------------------------------------------------------------

	WhichSampleArray.push_back("_reinteractions"); Universes.push_back(NG4Universes);	

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L mcc9_10_efficiency.cxx++");

	gROOT->ProcessLine(".L mcc9_10_response_matrices.cxx++");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {
	
		for (int k = 0; k < Universes[i]; k++) {

			gROOT->ProcessLine("mcc9_10_efficiency(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");		

			gROOT->ProcessLine("mcc9_10_response_matrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

		}

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
