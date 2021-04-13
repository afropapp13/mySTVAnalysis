#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TEfficiency.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "ubana/myClasses/Constants.h"
#include "ubana/myClasses/Util.h"
#include "ubana/myClasses/WienerSVD.h"

using namespace std;
using namespace Constants;

// -------------------------------------------------------------------------------------------------------------------------------------

void WienerSVD_Merge_Covariances(TString OverlaySample = "Overlay9") {

	// -------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();

	gStyle->SetOptStat(0);			

	// -------------------------------------------------------------------------------------

	vector<TString> PlotNames;
	PlotNames.push_back("DeltaPTPlot"); 
//	PlotNames.push_back("DeltaAlphaTPlot"); 
//	PlotNames.push_back("DeltaPhiTPlot");
//	PlotNames.push_back("MuonMomentumPlot"); 
//	PlotNames.push_back("MuonCosThetaPlot"); 
//	PlotNames.push_back("MuonPhiPlot");
//	PlotNames.push_back("ProtonMomentumPlot"); 
//	PlotNames.push_back("ProtonCosThetaPlot");
//	PlotNames.push_back("ProtonPhiPlot");
////	PlotNames.push_back("ECalPlot");
////	PlotNames.push_back("EQEPlot"); 
////	PlotNames.push_back("Q2Plot");
////	PlotNames.push_back("kMissPlot");
////	PlotNames.push_back("PMissPlot");
////	PlotNames.push_back("PMissMinusPlot");

	const int NPlots = PlotNames.size();
	//cout << "Number of 1D Plots = " << NPlots << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");
////	Runs.push_back("Run2");
//	Runs.push_back("Run3");
////	Runs.push_back("Run4");
////	Runs.push_back("Run5");				

	int NRuns = (int)(Runs.size());
	//cout << "Number of Runs = " << NRuns << endl;

	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> UncSources;
	UncSources.push_back("Stat");
	UncSources.push_back("POT");
	UncSources.push_back("NTargets");

	int NSamples = UncSources.size();

	vector<TFile*> CovFiles;
	CovFiles.resize(NSamples);

	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<vector<TH2D*> > Covariances;
	Covariances.resize(NPlots,vector<TH2D*>(NSamples));

	// -------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {					

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString TotalFileCovarianceSpecName = "WienerSVD_Total_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TString TotalFileCovarianceName = MigrationMatrixPath + TotalFileCovarianceSpecName;
		TFile* TotalFileCovarianceMatrices = new TFile(TotalFileCovarianceName,"recreate");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

			// Loop over the samples

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				TString FileCovarianceSpecName = "WienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
				TString FileCovarianceName = MigrationMatrixPath + FileCovarianceSpecName;
				CovFiles[WhichSample] = new TFile(FileCovarianceName,"readonly");

				if (WhichSample == 0) { Covariances[WhichPlot][WhichSample] = (TH2D*)(CovFiles[WhichSample]->Get("Covariance_"+PlotNames[WhichPlot]) ); }
				else { Covariances[WhichPlot][0]->Add( (TH2D*)(CovFiles[WhichSample]->Get("Covariance_"+PlotNames[WhichPlot]) ) ); }

			} // End of the loop over the samples

			TotalFileCovarianceMatrices->cd();
			Covariances[WhichPlot][0]->Write("TotalCovariance_"+PlotNames[WhichPlot]);
		
		} // End of the loop over the plots

		TotalFileCovarianceMatrices->Close();

	} // End of the loop over the runs	

} // End of the program 
