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
	PlotNames.push_back("MuonCosTheta");
	PlotNames.push_back("DeltaPT");	
	PlotNames.push_back("SerialDeltaPT_DeltaAlphaT");
	PlotNames.push_back("SerialDeltaPT_MuonCosTheta");
	PlotNames.push_back("SerialDeltaPT_ProtonCosTheta");		
	PlotNames.push_back("DeltaAlphaT");	
	PlotNames.push_back("SerialDeltaAlphaT_DeltaPT");
	PlotNames.push_back("SerialDeltaAlphaT_MuonCosTheta");
	PlotNames.push_back("SerialDeltaAlphaT_ProtonCosTheta");	
	PlotNames.push_back("DeltaPhiT");	
	PlotNames.push_back("SerialDeltaPhiT_DeltaPT");	
	PlotNames.push_back("DeltaPtx");	
	PlotNames.push_back("SerialDeltaPtx_DeltaPty");		
	PlotNames.push_back("ECal");	
	PlotNames.push_back("SerialECal_DeltaPT");
	PlotNames.push_back("SerialECal_DeltaAlphaT");	
	PlotNames.push_back("SerialECal_DeltaPty");

	const int NPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << NPlots << endl;

	//----------------------------------------//

		// Open the file that contains all the xsecs

		TString xsecname = "../../myXSec/v08_00_00_52/GenXSec/All_XSecs_Combined_v08_00_00_52.root";
		TFile* fxsec = new TFile(xsecname,"readonly");

		//----------------------------------------//

		// Data release

		TString drname = "/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance_PRD_Rename/DataRelease.root";
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