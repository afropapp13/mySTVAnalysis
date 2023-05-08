#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

using namespace std;

//----------------------------------------//

void DataReleaseRootFile() {

	//----------------------------------------//

	vector<TString> PlotNames;
	PlotNames.push_back("DeltaPn");	
	PlotNames.push_back("SerialDeltaPn_DeltaAlpha3Dq");		
	PlotNames.push_back("DeltaAlpha3Dq");	
	PlotNames.push_back("SerialDeltaAlpha3Dq_DeltaPn");
	PlotNames.push_back("DeltaPhi3D");	
	PlotNames.push_back("DeltaPnPar");
	PlotNames.push_back("DeltaPnPerp");	
	PlotNames.push_back("DeltaPnPerpx");
	PlotNames.push_back("DeltaPnPerpy");		

	const int NPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << NPlots << endl;

	//----------------------------------------//

		// Open the file that contains all the xsecs

		TString xsecname = "../../myXSec/v08_00_00_52/GenXSec/All_XSecs_Combined_v08_00_00_52.root";
		TFile* fxsec = new TFile(xsecname,"readonly");

		//----------------------------------------//

		// Data release

		TString drname = "/home/afroditi/Dropbox/Apps/Overleaf/General Imbalance Variables for Measuring Nuclear Effects and Demonstration with MicroBooNE Data/figures/gen/Afro/DataRelease.root";
		TFile* fdr = new TFile(drname,"recreate");

		//----------------------------------------//		

		// Loop over the plots

		for (int iplot = 0; iplot < NPlots; iplot ++) {								

			TH1D* FullUnc = (TH1D*)( fxsec->Get("FullUnc_" + PlotNames[iplot]+"Plot") );
			TH2D* Cov = (TH2D*)fxsec->Get("UnfCov_"+PlotNames[iplot]+"Plot");				
			TH2D* Ac = (TH2D*)fxsec->Get("Ac_"+PlotNames[iplot]+"Plot");		

			fdr->cd();
			FullUnc->Write("TotalUnc_" + PlotNames[iplot]);
			Cov->Write("Cov_" + PlotNames[iplot]);
			Ac->Write("Ac_" + PlotNames[iplot]);						

		} // End of the loop over the plots

		fdr->Close();

} // End of the program 