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

#include "ubana/AnalysisCode/Secondary_Code/GlobalSettings.cpp"

// -----------------------------------------------------------------------------------------------

double IntegratedXSecError(TH2D* LocalCovMatrix) {

	int n = LocalCovMatrix->GetNbinsX();
	double Nuedges[n+1];
			    
	for (int i = 0; i < n+1; i++) { Nuedges[i] = LocalCovMatrix->GetXaxis()->GetBinLowEdge(i+1); }

	TH1D* unc = new TH1D("unc","",n,Nuedges);

	for (int i = 1; i <= n; i++) { 

		double CovValue = LocalCovMatrix->GetBinContent(i,i);	
		unc->SetBinContent(i,TMath::Sqrt(CovValue)*100.);	
		
	}

	double IntegratedXSecErrorSquared = 0;

	for (int WhichXBin = 0; WhichXBin < n; WhichXBin++) {

		double BinWidth = unc->GetBinWidth(WhichXBin+1);
		double BinError = unc->GetBinContent(WhichXBin+1);

		cout << "Bin " << WhichXBin + 1 << "  BinWidth = " << BinWidth << "  BinError = " << BinError << endl;
		IntegratedXSecErrorSquared += TMath::Power(BinError * BinWidth,2.);
		
	}

	double IntegratedXSecErrorValue = TMath::Sqrt(IntegratedXSecErrorSquared);

	return IntegratedXSecErrorValue;

}

// -----------------------------------------------------------------------------------------------

void WienerSVD_QuantifyUnc(TString OverlaySample = "Overlay9") {

	// -------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();

	gStyle->SetOptStat(0);	

	// -------------------------------------------------------------------------------------

	GlobalSettings();
	TGaxis::SetMaxDigits(3);			

	// -------------------------------------------------------------------------------------

	vector<TString> PlotNames;
//	PlotNames.push_back("DeltaPTPlot"); 
	PlotNames.push_back("DeltaAlphaTPlot"); 
//	PlotNames.push_back("DeltaPhiTPlot");
//	PlotNames.push_back("MuonMomentumPlot"); 
//	PlotNames.push_back("MuonCosThetaPlot"); 
//	PlotNames.push_back("MuonPhiPlot");
//	PlotNames.push_back("ProtonMomentumPlot"); 
//	PlotNames.push_back("ProtonCosThetaPlot");
//	PlotNames.push_back("ProtonPhiPlot");
//	PlotNames.push_back("ECalPlot");
//	PlotNames.push_back("EQEPlot"); 
//	PlotNames.push_back("Q2Plot");

	const int NPlots = PlotNames.size();
	//cout << "Number of 1D Plots = " << NPlots << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");
//	Runs.push_back("Run2");
//	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");				

	int NRuns = (int)(Runs.size());
	//cout << "Number of Runs = " << NRuns << endl;

	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> UncSources;

	UncSources.push_back("Stat");
	UncSources.push_back("POT");
	UncSources.push_back("NTarget");
	UncSources.push_back("LY");
	UncSources.push_back("TPC");
	UncSources.push_back("XSec");
	UncSources.push_back("G4");
	UncSources.push_back("Flux");
	UncSources.push_back("Dirt");

	UncSources.push_back("MC_Stat");
	UncSources.push_back("MC_POT");
	UncSources.push_back("MC_NTarget");
	UncSources.push_back("MC_LY");
	UncSources.push_back("MC_TPC");
	UncSources.push_back("MC_XSec");
	UncSources.push_back("MC_G4");
	UncSources.push_back("MC_Flux");
	UncSources.push_back("MC_Dirt");

	int NSamples = UncSources.size();

	vector<TFile*> CovFiles;
	CovFiles.resize(NSamples);

	vector<TFile*> XSecFiles;
	XSecFiles.resize(NSamples);	

	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TH2D*> Covariances;
	Covariances.resize(NPlots);

	vector<TH2D*> StatCovariances;
	StatCovariances.resize(NPlots);

	vector<TH2D*> SystCovariances;
	SystCovariances.resize(NPlots);		

	// -------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {					

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString TotalFileCovarianceSpecName = "WienerSVD_Total_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TString TotalFileCovarianceName = MigrationMatrixPath + TotalFileCovarianceSpecName;

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

			// Loop over the samples

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				TString FileCovarianceSpecName = "WienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
				TString FileCovarianceName = MigrationMatrixPath + FileCovarianceSpecName;
				CovFiles[WhichSample] = new TFile(FileCovarianceName,"readonly");

				TH2D* LocalCovMatrix = (TH2D*)(CovFiles[WhichSample]->Get(UncSources[WhichSample]+"_Covariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]));

				if (string(UncSources[WhichSample]).find("Stat") != std::string::npos) { 

					if (UncSources[WhichSample] == "Stat") { StatCovariances[WhichPlot] = LocalCovMatrix; }
					else { StatCovariances[WhichPlot]->Add(LocalCovMatrix); }
				
				} else {

					if (WhichSample == 1) { 

						SystCovariances[WhichPlot] = LocalCovMatrix; 

					}
					else if (WhichSample > 1) { 
						
						SystCovariances[WhichPlot]->Add(LocalCovMatrix); 
						
					}
					
				}				

			} // End of the loop over the samples

			// ------------------------------------------------------------------

			Covariances[WhichPlot] = (TH2D*) (StatCovariances[WhichPlot]->Clone());
			Covariances[WhichPlot]->Add(SystCovariances[WhichPlot]);				

			// ------------------------------------------------------------------

			// Print out the contribution for each one of the uncertainties using the DeltaAlphaT plot

			//if (PlotNames[WhichPlot] == "DeltaAlphaTPlot") {

				double StatUnc = IntegratedXSecError(StatCovariances[WhichPlot]);

				cout << Runs[WhichRun] << "  " << PlotNames[WhichPlot] << endl;

				cout << "Stat Unc = " << StatUnc << endl;

			//}			

			cout << endl << endl;

			// ------------------------------------------------------------------

		
		} // End of the loop over the plots

	} // End of the loop over the runs	

} // End of the program 
