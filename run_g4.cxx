{

	vector<TString> WhichSampleArray;
	vector<int> Universes;
	int NG4Universes = 100;
	
	// -----------------------------------------------------------------------------------------

	WhichSampleArray.push_back("_reinteractions"); Universes.push_back(NG4Universes);	

	// -----------------------------------------------------------------------------------------

	gROOT->ProcessLine(".L efficiency.cxx++");

	gROOT->ProcessLine(".L response_matrices.cxx++");

	for (int i = 0;i < (int)(WhichSampleArray.size()); i++) {
	
		for (int k = 0; k < Universes[i]; k++) {

			gROOT->ProcessLine("efficiency(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");		

			gROOT->ProcessLine("response_matrices(\""+WhichSampleArray[i]+"_"+TString(std::to_string(k))+"\")");

		}

	}

	// -----------------------------------------------------------------------------------------

	//gROOT->ProcessLine(".q");

}
