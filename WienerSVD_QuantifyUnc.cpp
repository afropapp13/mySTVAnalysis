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

double IntegratedXSecError(TH2D* FracCovMatrix,TH1D* CV) {

	int n = FracCovMatrix->GetNbinsX();

	double IntegratedXSecErrorSquared = 0;
	double IntegratedXSec = 0;

	for (int i = 1; i <= n; i++) { 

		double BinWidth = CV->GetBinWidth(i);
		double BinCV = CV->GetBinContent(i);

		//double TotalBinError = CV->GetBinError(i); // Total uncertainty in that bin	
		double FracError = TMath::Sqrt(FracCovMatrix->GetBinContent(i,i)); // Fractional contribution from a given source of uncertainty
		double BinError = TMath::Abs(BinCV*FracError);

		IntegratedXSecErrorSquared += TMath::Power(BinError * BinWidth,2.);
		IntegratedXSec += BinCV * BinWidth;

		//cout << "Bin " << i << "  BinWidth = " << BinWidth  << " CV = " << BinCV << "  BinError = " << BinError << endl;
		
	}

	double IntegratedXSecErrorValue = TMath::Sqrt(IntegratedXSecErrorSquared) / IntegratedXSec * 100.;

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
	std::cout.precision(3);			

	// -------------------------------------------------------------------------------------

	vector<TString> PlotNames;
//	PlotNames.push_back("DeltaPTPlot"); 
//	PlotNames.push_back("DeltaAlphaTPlot"); 
//	PlotNames.push_back("DeltaPhiTPlot");
//	PlotNames.push_back("MuonMomentumPlot"); 
	PlotNames.push_back("MuonCosThetaPlot"); 
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
	Runs.push_back("Run3");
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
//	UncSources.push_back("MC_POT");
//	UncSources.push_back("MC_NTarget");
//	UncSources.push_back("MC_LY");
//	UncSources.push_back("MC_TPC");
//	UncSources.push_back("MC_XSec");
//	UncSources.push_back("MC_G4");
//	UncSources.push_back("MC_Flux");
//	UncSources.push_back("MC_Dirt");

//	UncSources.push_back("SmEff_Stat");
	UncSources.push_back("SmEff_POT");
	UncSources.push_back("SmEff_NTarget");
	UncSources.push_back("SmEff_LY");
	UncSources.push_back("SmEff_TPC");
	UncSources.push_back("SmEff_XSec");
	UncSources.push_back("SmEff_G4");
	UncSources.push_back("SmEff_Flux");
//	UncSources.push_back("SmEff_Dirt");

	int NSamples = UncSources.size();

	vector<TFile*> CovFiles;
	CovFiles.resize(NSamples);

	vector<TFile*> XSecFiles;
	XSecFiles.resize(NSamples);	

	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TH2D*> FracCovariances;
	FracCovariances.resize(NPlots);

	vector<TH2D*> TotalCovariances;
	TotalCovariances.resize(NPlots);

	vector<TH2D*> StatFracCovariances;
	StatFracCovariances.resize(NPlots);

	vector<TH2D*> POTFracCovariances;
	POTFracCovariances.resize(NPlots);

	vector<TH2D*> NTargetFracCovariances;
	NTargetFracCovariances.resize(NPlots);	

	vector<TH2D*> LYFracCovariances;
	LYFracCovariances.resize(NPlots);		

	vector<TH2D*> TPCFracCovariances;
	TPCFracCovariances.resize(NPlots);	

	vector<TH2D*> XSecFracCovariances;
	XSecFracCovariances.resize(NPlots);

	vector<TH2D*> G4FracCovariances;
	G4FracCovariances.resize(NPlots);	

	vector<TH2D*> FluxFracCovariances;
	FluxFracCovariances.resize(NPlots);

	vector<TH2D*> DirtFracCovariances;
	DirtFracCovariances.resize(NPlots);

	// -------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {					

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString XSecFileName = PathToExtractedXSec+"WienerSVD_ExtractedXSec_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TFile* ExtractedXSec = TFile::Open(XSecFileName,"readonly");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

			// Total Covariance Matrix (not fractional covariance matrix)

			TotalCovariances[WhichPlot] = (TH2D*)(ExtractedXSec->Get("UnfCov"+PlotNames[WhichPlot]));

			// Loop over the samples

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				// Fractional Covariances

				TString FileCovarianceSpecName = "WienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
				TString FileCovarianceName = MigrationMatrixPath + FileCovarianceSpecName;
				CovFiles[WhichSample] = new TFile(FileCovarianceName,"readonly");

				TH2D* LocalFracCovMatrix = (TH2D*)(CovFiles[WhichSample]->Get(UncSources[WhichSample]+"_Covariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]));

				if (string(UncSources[WhichSample]).find("Stat") != std::string::npos) { 

					if (UncSources[WhichSample] == "Stat") { StatFracCovariances[WhichPlot] = LocalFracCovMatrix; }
					else { StatFracCovariances[WhichPlot]->Add(LocalFracCovMatrix); }
				
				} else if (string(UncSources[WhichSample]).find("POT") != std::string::npos) { 

					if (UncSources[WhichSample] == "POT") { POTFracCovariances[WhichPlot] = LocalFracCovMatrix; }
					else { POTFracCovariances[WhichPlot]->Add(LocalFracCovMatrix); }
					
				} else if (string(UncSources[WhichSample]).find("NTarget") != std::string::npos) { 

					if (UncSources[WhichSample] == "NTarget") { NTargetFracCovariances[WhichPlot] = LocalFracCovMatrix; }
					else { NTargetFracCovariances[WhichPlot]->Add(LocalFracCovMatrix); }
					
				} else if (string(UncSources[WhichSample]).find("LY") != std::string::npos) { 

					if (UncSources[WhichSample] == "LY") { LYFracCovariances[WhichPlot] = LocalFracCovMatrix; }
					else { LYFracCovariances[WhichPlot]->Add(LocalFracCovMatrix); }
					
				} else if (string(UncSources[WhichSample]).find("TPC") != std::string::npos) { 

					if (UncSources[WhichSample] == "TPC") { TPCFracCovariances[WhichPlot] = LocalFracCovMatrix; }
					else { TPCFracCovariances[WhichPlot]->Add(LocalFracCovMatrix); }
					
				} else if (string(UncSources[WhichSample]).find("XSec") != std::string::npos) { 

					if (UncSources[WhichSample] == "XSec") { XSecFracCovariances[WhichPlot] = LocalFracCovMatrix; }
					else { XSecFracCovariances[WhichPlot]->Add(LocalFracCovMatrix); }
					
				} else if (string(UncSources[WhichSample]).find("G4") != std::string::npos) { 

					if (UncSources[WhichSample] == "G4") { G4FracCovariances[WhichPlot] = LocalFracCovMatrix; }
					else { G4FracCovariances[WhichPlot]->Add(LocalFracCovMatrix); }
					
				} else if (string(UncSources[WhichSample]).find("Flux") != std::string::npos) { 

					if (UncSources[WhichSample] == "Flux") { FluxFracCovariances[WhichPlot] = LocalFracCovMatrix; }
					else { FluxFracCovariances[WhichPlot]->Add(LocalFracCovMatrix); }
					
				} else if (string(UncSources[WhichSample]).find("Dirt") != std::string::npos) { 

					if (UncSources[WhichSample] == "Dirt") { DirtFracCovariances[WhichPlot] = LocalFracCovMatrix; }
					else { DirtFracCovariances[WhichPlot]->Add(LocalFracCovMatrix); }
					
				}

			} // End of the loop over the samples				

			// ------------------------------------------------------------------

			TH1D* CV = (TH1D*)(ExtractedXSec->Get("Reco"+PlotNames[WhichPlot]));

			// ------------------------------------------------------------------

			// Print out the contribution for each one of the uncertainties using the DeltaAlphaT plot

			//if (PlotNames[WhichPlot] == "DeltaAlphaTPlot") {

				double StatUnc = IntegratedXSecError(StatFracCovariances[WhichPlot],CV);
				double POTUnc = IntegratedXSecError(POTFracCovariances[WhichPlot],CV);
				double NTargetUnc = IntegratedXSecError(NTargetFracCovariances[WhichPlot],CV);
				double LYUnc = IntegratedXSecError(LYFracCovariances[WhichPlot],CV);
				double TPCUnc = IntegratedXSecError(TPCFracCovariances[WhichPlot],CV);
				double XSecUnc = IntegratedXSecError(XSecFracCovariances[WhichPlot],CV);
				double G4Unc = IntegratedXSecError(G4FracCovariances[WhichPlot],CV);
				double FluxUnc = IntegratedXSecError(FluxFracCovariances[WhichPlot],CV);
				double DirtUnc = IntegratedXSecError(DirtFracCovariances[WhichPlot],CV);

				cout << Runs[WhichRun] << "  " << PlotNames[WhichPlot] << endl << endl;

				cout << "Stat & " << StatUnc << "\\\\" << endl;
				cout << "%POT & " << POTUnc << "\\\\" << endl;
				cout << "%NTarget & " << NTargetUnc << "\\\\" << endl;
				cout << "LY & " << LYUnc  << "\\\\" << endl;
				cout << "TPC & " << TPCUnc  << "\\\\" << endl;
				cout << "XSec & " << XSecUnc  << "\\\\" << endl;
				cout << "G4 & " << G4Unc  << "\\\\" << endl;
				cout << "Flux & " << FluxUnc  << "\\\\" << endl;
				cout << "Dirt & " << DirtUnc  << "\\\\" << endl;

				double TotalUnc = TMath::Sqrt(TMath::Power(StatUnc,2.) + TMath::Power(POTUnc,2.) + TMath::Power(NTargetUnc,2.) + TMath::Power(LYUnc,2.) + TMath::Power(TPCUnc,2.) + TMath::Power(XSecUnc,2.) + TMath::Power(G4Unc,2.)  + TMath::Power(FluxUnc,2.) + TMath::Power(DirtUnc,2.));
				cout << "\\hline" << endl << "\\hline" << endl;
				cout << "Total & " << TotalUnc  << "\\\\" << endl << endl;

			//}			

			// ------------------------------------------------------------------

		
		} // End of the loop over the plots

	} // End of the loop over the runs	

} // End of the program 
