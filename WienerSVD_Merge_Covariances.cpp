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

void ReturnUncPlot(TH2D* LocalCovMatrix,TString PlotName, TString Run,TString UncSources,int Color, TLegend* leg) {

	int n = LocalCovMatrix->GetNbinsX();
	TString TitleX =  LocalCovMatrix->GetXaxis()->GetTitle();
	double Nuedges[n+1];
			    
	for (int i = 0; i < n+1; i++) { Nuedges[i] = LocalCovMatrix->GetXaxis()->GetBinLowEdge(i+1); }

	TH1D* unc = new TH1D("unc_"+PlotName+"_"+Run+"_"+UncSources,";"+TitleX+";"+Run+" Uncertainty [%]",n,Nuedges);	

	for (int i = 1; i <= n; i++) { 

		double CovValue = LocalCovMatrix->GetBinContent(i,i);	
		unc->SetBinContent(i,TMath::Sqrt(CovValue)*100.);
		
	}

	unc->GetXaxis()->SetTitleFont(FontStyle);
	unc->GetXaxis()->SetLabelFont(FontStyle);
	unc->GetXaxis()->SetTitleSize(TextSize);
	unc->GetXaxis()->SetLabelSize(TextSize);	
	unc->GetXaxis()->CenterTitle();
	unc->GetXaxis()->SetNdivisions(5);
	
	unc->GetYaxis()->SetLabelFont(FontStyle);
	unc->GetYaxis()->SetTitleFont(FontStyle);
	unc->GetYaxis()->SetTitleSize(TextSize);
	unc->GetYaxis()->SetLabelSize(TextSize);	
	unc->GetYaxis()->CenterTitle();
	unc->GetYaxis()->SetNdivisions(5);
	unc->GetYaxis()->SetTitleOffset(1.);				
	unc->GetYaxis()->SetRangeUser(0.,99.);

	if (PlotName == "DeltaPTPlot") { unc->GetYaxis()->SetRangeUser(0.,32.); }
	if (PlotName == "DeltaAlphaTPlot") { unc->GetYaxis()->SetRangeUser(0.,14.); }
	if (PlotName == "DeltaPhiTPlot") { unc->GetYaxis()->SetRangeUser(0.,47.); }
	if (PlotName == "MuonMomentumPlot") { unc->GetYaxis()->SetRangeUser(0.,37.); }
	if (PlotName == "MuonPhiPlot") { unc->GetYaxis()->SetRangeUser(0.,16.); }
	if (PlotName == "ProtonPhiPlot") { unc->GetYaxis()->SetRangeUser(0.,16.); }
	if (PlotName == "ProtonMomentumPlot") { unc->GetYaxis()->SetRangeUser(0.,19.); }	
				


	unc->SetLineWidth(2);
	if (Color >= 9) { Color = Color - 9; }
	if (Color >= 18) { Color = Color - 18; }		
	unc->SetLineColor(Color+1);
	unc->SetMarkerColor(Color+1);			

	unc->Draw("hist text0 same");

	leg->AddEntry(unc,UncSources);

}

// -----------------------------------------------------------------------------------------------

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

			// Create canvases for 3 categories of systematics: MC event rates, Data event rates, Smearing + Efficiency

			TString DataERCanvasName = "DataERSyst_"+PlotNames[WhichPlot]+OverlaySample+"_"+Runs[WhichRun];
			TCanvas* DataERPlotCanvas = new TCanvas(DataERCanvasName,DataERCanvasName,205,34,1024,768);
			DataERPlotCanvas->SetBottomMargin(0.16);
			DataERPlotCanvas->SetLeftMargin(0.15);
			DataERPlotCanvas->SetRightMargin(0.25);
			DataERPlotCanvas->SetTopMargin(0.15);			

			TLegend* legData = new TLegend(0.02,0.89,0.97,0.99);
			legData->SetBorderSize(0);
			legData->SetTextSize(0.04);
			legData->SetTextFont(FontStyle);
			legData->SetNColumns(5);			

			TString MCERCanvasName = "MCERSyst_"+PlotNames[WhichPlot]+OverlaySample+"_"+Runs[WhichRun];
			TCanvas* MCERPlotCanvas = new TCanvas(MCERCanvasName,MCERCanvasName,205,34,1024,768);
			MCERPlotCanvas->SetBottomMargin(0.16);
			MCERPlotCanvas->SetLeftMargin(0.15);
			MCERPlotCanvas->SetRightMargin(0.25);
			MCERPlotCanvas->SetTopMargin(0.15);			

			TLegend* legMC = new TLegend(0.02,0.89,0.97,0.99);
			legMC->SetBorderSize(0);
			legMC->SetTextSize(0.04);
			legMC->SetTextFont(FontStyle);
			legMC->SetNColumns(5);			

			TString SmEffCanvasName = "SmEffSyst_"+PlotNames[WhichPlot]+OverlaySample+"_"+Runs[WhichRun];
			TCanvas* SmEffPlotCanvas = new TCanvas(SmEffCanvasName,MCERCanvasName,205,34,1024,768);			
			SmEffPlotCanvas->SetBottomMargin(0.16);
			SmEffPlotCanvas->SetLeftMargin(0.15);
			SmEffPlotCanvas->SetRightMargin(0.25);
			SmEffPlotCanvas->SetTopMargin(0.15);			

			TLegend* legSmEff = new TLegend(0.02,0.89,0.97,0.99);
			legSmEff->SetBorderSize(0);
			legSmEff->SetTextSize(0.04);
			legSmEff->SetTextFont(FontStyle);
			legSmEff->SetNColumns(5);

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
						TotalFileCovarianceMatrices->cd();
						//SystCovariances[WhichPlot]->Write(UncSources[WhichSample]+"Covariance_"+PlotNames[WhichPlot]);

					}
					else if (WhichSample > 1) { 
						
						TotalFileCovarianceMatrices->cd();
						//LocalCovMatrix->Write(UncSources[WhichSample]+"Covariance_"+PlotNames[WhichPlot]);
						SystCovariances[WhichPlot]->Add(LocalCovMatrix); 
						
					}
					
				}

				// --------------------------------------------------------------------------------

				TLegend* leg = nullptr;

				if (string(UncSources[WhichSample]).find("SmEff_") != std::string::npos) { SmEffPlotCanvas->cd(); leg = legSmEff;}
				else if (string(UncSources[WhichSample]).find("MC_") != std::string::npos) { MCERPlotCanvas->cd();  leg = legMC; }				
				else { DataERPlotCanvas->cd(); leg = legData; }

				ReturnUncPlot(LocalCovMatrix,PlotNames[WhichPlot],Runs[WhichRun],UncSources[WhichSample],WhichSample,leg);				

			} // End of the loop over the samples

			// ------------------------------------------------------------------

			DataERPlotCanvas->cd();
			legData->Draw();

			MCERPlotCanvas->cd();
			legMC->Draw();

			SmEffPlotCanvas->cd();
			legSmEff->Draw();

			// ------------------------------------------------------------------

			DataERPlotCanvas->SaveAs(PlotPath+OverlaySample+"/"+DataERCanvasName+".pdf");
			MCERPlotCanvas->SaveAs(PlotPath+OverlaySample+"/"+MCERCanvasName+".pdf");
			SmEffPlotCanvas->SaveAs(PlotPath+OverlaySample+"/"+SmEffCanvasName+".pdf");	

			delete 	DataERPlotCanvas;	
			delete 	MCERPlotCanvas;
			delete 	SmEffPlotCanvas;

			// ------------------------------------------------------------------			

			// Storing Total Covariance Matrices

			TotalFileCovarianceMatrices->cd();
			StatCovariances[WhichPlot]->Write("StatCovariance_"+PlotNames[WhichPlot]);
			SystCovariances[WhichPlot]->Write("SystCovariance_"+PlotNames[WhichPlot]);						

			Covariances[WhichPlot] = (TH2D*) (StatCovariances[WhichPlot]);
			Covariances[WhichPlot]->Add(SystCovariances[WhichPlot]);			
			Covariances[WhichPlot]->Write("TotalCovariance_"+PlotNames[WhichPlot]);
			
			// ---------------------------------------------------------------------------------------------
			
			// Plotting Total Fractional Covariances

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

			double CovMax = TMath::Min(1.,1.05*Covariances[WhichPlot]->GetMaximum());
			double CovMin = TMath::Min(0.,1.05*Covariances[WhichPlot]->GetMinimum());
			Covariances[WhichPlot]->GetZaxis()->SetRangeUser(CovMin,CovMax);
//			Covariances[WhichPlot]->GetZaxis()->SetTitle("[x10^{-76} cm^{4}]");
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
