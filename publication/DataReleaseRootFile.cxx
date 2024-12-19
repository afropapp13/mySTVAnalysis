#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

//----------------------------------------//

void DataReleaseRootFile() {

	//----------------------------------------//

	vector<TString> PlotNames;
	PlotNames.push_back("ThetaVis");
	PlotNames.push_back("SerialThetaVis_ECal");
	PlotNames.push_back("SerialThetaVis_DeltaPn");
	PlotNames.push_back("SerialThetaVis_PMiss");		

	const int NPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << NPlots << endl;

	//----------------------------------------//

		// Open the file that contains all the xsecs

		TString xsecname = PathToExtractedXSec+"/WienerSVD_ExtractedXSec_Overlay9_Combined_"+UBCodeVersion+".root";
		TFile* fxsec = new TFile(xsecname,"readonly");

		//----------------------------------------//

		// Data release

		TString drname = "DataRelease.root";
		TFile* fdr = new TFile(drname,"recreate");

		//----------------------------------------//		

		// Loop over the plots

		for (int iplot = 0; iplot < NPlots; iplot ++) {								

			TH1D* FullUnc = (TH1D*)( fxsec->Get("Reco" + PlotNames[iplot]+"Plot") );
			TH2D* Cov = (TH2D*)fxsec->Get("UnfCov"+PlotNames[iplot]+"Plot");				
			TH2D* Ac = (TH2D*)fxsec->Get("Ac"+PlotNames[iplot]+"Plot");		

			fdr->cd();
			FullUnc->Write("TotalUnc_" + PlotNames[iplot]);
			Cov->Write("Cov_" + PlotNames[iplot]);
			Ac->Write("Ac_" + PlotNames[iplot]);						

		} // End of the loop over the plots

		fdr->Close();

} // End of the program 
