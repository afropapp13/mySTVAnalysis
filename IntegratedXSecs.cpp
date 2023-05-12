#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "/home/afroditi/Dropbox/PhD/Secondary_Code/mySimFunctions.cpp"
#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

//----------------------------------------//

void IntegratedXSecs() {

	//----------------------------------------//

	TH1D::SetDefaultSumw2();
	std::cout.precision(4);

	//----------------------------------------//

	vector<TString> PlotNames;

//	PlotNames.push_back("DeltaPtxPlot"); 
	PlotNames.push_back("DeltaPTPlot"); 
	PlotNames.push_back("DeltaAlphaTPlot"); 
	PlotNames.push_back("DeltaPhiTPlot");
//	PlotNames.push_back("MuonMomentumPlot"); 
//	PlotNames.push_back("MuonCosThetaPlot");
//	PlotNames.push_back("MuonCosThetaSingleBinPlot");	 
//	PlotNames.push_back("MuonPhiPlot");
//	PlotNames.push_back("ProtonMomentumPlot"); 
//	PlotNames.push_back("ProtonCosThetaPlot");
//	PlotNames.push_back("ProtonPhiPlot");
//	PlotNames.push_back("VertexXPlot");
//	PlotNames.push_back("VertexYPlot");	
//	PlotNames.push_back("VertexZPlot");

//	PlotNames.push_back("SerialDeltaPT_DeltaAlphaTPlot_0");
//	PlotNames.push_back("SerialDeltaPT_DeltaAlphaTPlot_1");
//	PlotNames.push_back("SerialDeltaPT_DeltaAlphaTPlot_2");
//	PlotNames.push_back("SerialDeltaPT_DeltaAlphaTPlot_3");			

//	PlotNames.push_back("SerialDeltaAlphaT_DeltaPTPlot_0");
//	PlotNames.push_back("SerialDeltaAlphaT_DeltaPTPlot_1");	
//	PlotNames.push_back("SerialDeltaAlphaT_DeltaPTPlot_2");

//	PlotNames.push_back("SerialDeltaPtx_DeltaPtyPlot_0");
//	PlotNames.push_back("SerialDeltaPtx_DeltaPtyPlot_1");
//	PlotNames.push_back("SerialDeltaPtx_DeltaPtyPlot_2");		

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	//----------------------------------------//

	vector<TString> Runs;
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	//----------------------------------------//

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

	//----------------------------------------//

		cout << Runs[WhichRun] << endl << endl;

		vector<TString> NameOfSamples; NameOfSamples.clear();
	
		NameOfSamples.push_back("Data");
		//NameOfSamples.push_back("OverlayGENIE");

		//NameOfSamples.push_back("Genie_v3_0_6_Out_Of_The_Box");
		//NameOfSamples.push_back("Genie_v3_0_6_uB_Tune_1");
		//NameOfSamples.push_back("SuSav2");
		//NameOfSamples.push_back("NuWro");
		//NameOfSamples.push_back("GiBUU");
		//NameOfSamples.push_back("GENIEv2");
		//NameOfSamples.push_back("NEUT");
		//NameOfSamples.push_back("GENIEv3_0_4");

		const int NSamples = NameOfSamples.size();

		vector< vector<TH1D*> > Plots; Plots.clear(); Plots.resize(N1DPlots);

		//----------------------------------------//

		TFile* FileSample = new TFile("myXSec/v08_00_00_52/GenXSec/All_XSecs_Combined_v08_00_00_52.root","readonly");

		for (int iplot = 0; iplot < N1DPlots; iplot++) {

			Plots[iplot].resize(NSamples);
			cout << PlotNames[iplot] << endl << endl;

			for (int isample = 0; isample < NSamples; isample++) {

				if (NameOfSamples[isample] == "Data") { Plots[iplot][isample] = (TH1D*)(FileSample->Get( "FullUnc_" + PlotNames[iplot] )); }
				else { Plots[iplot][isample] = (TH1D*)(FileSample->Get( NameOfSamples[isample] + "_" + PlotNames[iplot] )); }

				cout << NameOfSamples[isample] << " & " << IntegratedXSec(Plots[iplot][isample]) << endl;

			} // end of the loop over the samples

			cout << endl;

		} // End of the loop over the plots

	} // End of the loop over the runs	

} // End of the program 
