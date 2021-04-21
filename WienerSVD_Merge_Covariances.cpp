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

// -------------------------------------------------------------------------------------------------------------------------------------

void WienerSVD_Merge_Covariances(TString OverlaySample = "Overlay9") {

	// -------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();

	gStyle->SetOptStat(0);	

	// -------------------------------------------------------------------------------------

	GlobalSettings();
	TGaxis::SetMaxDigits(3);			

	// -------------------------------------------------------------------------------------

	vector<TString> PlotNames;
	PlotNames.push_back("DeltaPTPlot"); 
	PlotNames.push_back("DeltaAlphaTPlot"); 
	PlotNames.push_back("DeltaPhiTPlot");
	PlotNames.push_back("MuonMomentumPlot"); 
	PlotNames.push_back("MuonCosThetaPlot"); 
	PlotNames.push_back("MuonPhiPlot");
	PlotNames.push_back("ProtonMomentumPlot"); 
	PlotNames.push_back("ProtonCosThetaPlot");
	PlotNames.push_back("ProtonPhiPlot");
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

	int NSamples = UncSources.size();

	vector<TFile*> CovFiles;
	CovFiles.resize(NSamples);

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
		TFile* TotalFileCovarianceMatrices = new TFile(TotalFileCovarianceName,"recreate");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

			// Loop over the samples

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				TString FileCovarianceSpecName = "WienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
				TString FileCovarianceName = MigrationMatrixPath + FileCovarianceSpecName;
				CovFiles[WhichSample] = new TFile(FileCovarianceName,"readonly");

				if (WhichSample == 0) { 

					StatCovariances[WhichPlot] = (TH2D*)(CovFiles[WhichSample]->Get(UncSources[WhichSample]+"_Covariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]) ); 
				
				} else {

					if (WhichSample == 1) { 

						SystCovariances[WhichPlot] = (TH2D*)(CovFiles[WhichSample]->Get(UncSources[WhichSample]+"_Covariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]) ); 
						TotalFileCovarianceMatrices->cd();
						SystCovariances[WhichPlot]->Write(UncSources[WhichSample]+"Covariance_"+PlotNames[WhichPlot]);

					}
					else if (WhichSample > 1) { 
						
						TH2D* LocalCovMatrix = (TH2D*)(CovFiles[WhichSample]->Get(UncSources[WhichSample]+"_Covariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]));
						TotalFileCovarianceMatrices->cd();
						LocalCovMatrix->Write(UncSources[WhichSample]+"Covariance_"+PlotNames[WhichPlot]);
						SystCovariances[WhichPlot]->Add(LocalCovMatrix); 
						
					}
					
				}

			} // End of the loop over the samples

			TotalFileCovarianceMatrices->cd();
			StatCovariances[WhichPlot]->Write("StatCovariance_"+PlotNames[WhichPlot]);
			SystCovariances[WhichPlot]->Write("SystCovariance_"+PlotNames[WhichPlot]);						

			Covariances[WhichPlot] = (TH2D*) (StatCovariances[WhichPlot]);
			Covariances[WhichPlot]->Add(SystCovariances[WhichPlot]);			
			Covariances[WhichPlot]->Write("TotalCovariance_"+PlotNames[WhichPlot]);

			TString CanvasName = "Total_"+PlotNames[WhichPlot]+OverlaySample+"_"+Runs[WhichRun];
			TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.16);
			PlotCanvas->SetLeftMargin(0.15);
			PlotCanvas->SetRightMargin(0.25);			
			
			gStyle->SetMarkerSize(1.5);
			gStyle->SetPaintTextFormat("4.3f");			
			
			Covariances[WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
			Covariances[WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
			Covariances[WhichPlot]->GetXaxis()->SetTitleSize(TextSize);
			Covariances[WhichPlot]->GetXaxis()->SetLabelSize(TextSize);			
			Covariances[WhichPlot]->GetXaxis()->CenterTitle();
			Covariances[WhichPlot]->GetXaxis()->SetNdivisions(5);
			
			Covariances[WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
			Covariances[WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
			Covariances[WhichPlot]->GetYaxis()->SetTitleSize(TextSize);
			Covariances[WhichPlot]->GetYaxis()->SetLabelSize(TextSize);			
			Covariances[WhichPlot]->GetYaxis()->CenterTitle();
			Covariances[WhichPlot]->GetYaxis()->SetNdivisions(5);
			Covariances[WhichPlot]->GetYaxis()->SetTitleOffset(1.);						

			Covariances[WhichPlot]->SetTitle(Runs[WhichRun] + " Total");	

			double CovMax = 1.05*Covariances[WhichPlot]->GetMaximum();
			double CovMin = TMath::Min(0.,1.05*Covariances[WhichPlot]->GetMinimum());
			Covariances[WhichPlot]->GetZaxis()->SetRangeUser(CovMin,CovMax);
			Covariances[WhichPlot]->GetZaxis()->SetTitle("[x10^{-76} cm^{4}]");
			Covariances[WhichPlot]->GetZaxis()->CenterTitle();
			Covariances[WhichPlot]->GetZaxis()->SetTitleFont(FontStyle);
			Covariances[WhichPlot]->GetZaxis()->SetTitleSize(TextSize);
			Covariances[WhichPlot]->GetZaxis()->SetLabelFont(FontStyle);
			Covariances[WhichPlot]->GetZaxis()->SetLabelSize(TextSize-0.01);
			Covariances[WhichPlot]->GetZaxis()->SetNdivisions(5);

			Covariances[WhichPlot]->SetMarkerColor(kWhite);			
			Covariances[WhichPlot]->SetMarkerSize(1.5);
//			Covariances[WhichPlot]->Draw("text colz e"); 
			Covariances[WhichPlot]->Draw("colz");
			
			PlotCanvas->SaveAs(PlotPath+OverlaySample+"/WienerSVD_Total_CovarianceMatrices_"+PlotNames[WhichPlot]+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
			
			delete PlotCanvas;
		
		} // End of the loop over the plots

		TotalFileCovarianceMatrices->Close();
		cout << endl << "Covariance matrix file " << TotalFileCovarianceName << " has been created" << endl << endl;

		cout << "Merging of covariance matrices for run " << Runs[WhichRun] << " completed!" << endl;

	} // End of the loop over the runs	

} // End of the program 
