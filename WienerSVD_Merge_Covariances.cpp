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

//	if (string(UncSources).find("POT") != std::string::npos || string(UncSources).find("NTarget") != std::string::npos) {  }
//	else {

	int n = LocalCovMatrix->GetNbinsX();
	TString TitleX =  LocalCovMatrix->GetXaxis()->GetTitle();
	double Nuedges[n+1];
			    
	for (int i = 0; i < n+1; i++) { Nuedges[i] = LocalCovMatrix->GetXaxis()->GetBinLowEdge(i+1); }

	TH1D* unc = new TH1D("unc_"+PlotName+"_"+Run+"_"+UncSources,";"+TitleX+";"+Run+" Uncertainty [%]",n,Nuedges);	

	for (int i = 1; i <= n; i++) { 

		double CovValue = LocalCovMatrix->GetBinContent(i,i);	
		unc->SetBinContent(i,TMath::Sqrt(CovValue)*100.);
		
//	}

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

}

// -----------------------------------------------------------------------------------------------

void WienerSVD_Merge_Covariances(TString OverlaySample = "Overlay9", TString BeamOn9 = "", bool IncludeUnfTech = false) {

	// -------------------------------------------------------------------------------------

	// Store the pdf's only for the OverlaySample = "Overlay9" & BeamOn9 = ""
	// Otherwise, just store the merged fractional covariance matrix	

	bool StorePlots = false;
	if (OverlaySample == "Overlay9" && BeamOn9 == "") { StorePlots = true; }

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

	PlotNames.push_back("CCQEMuonMomentumPlot"); 
	PlotNames.push_back("CCQEMuonCosThetaPlot"); 
	PlotNames.push_back("CCQEProtonMomentumPlot"); 
	PlotNames.push_back("CCQEProtonCosThetaPlot");

	const int NPlots = PlotNames.size();
	//cout << "Number of 1D Plots = " << NPlots << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");				
	Runs.push_back("Combined");				

	int NRuns = (int)(Runs.size());
	//cout << "Number of Runs = " << NRuns << endl;

	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> UncSources;

	if (BeamOn9 != "") { // Fake data study, we need only the Stat, MC_Stat, XSec & SmEff_XSec covariances

	UncSources.push_back("Stat");
	UncSources.push_back("XSec");
	UncSources.push_back("MC_Stat");
 	UncSources.push_back("SmEff_XSec");

	} else {

	UncSources.push_back("Stat");
	UncSources.push_back("POT");     // Potentially add back 
	UncSources.push_back("NTarget"); // Potentially add back
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
//	UncSources.push_back("SmEff_POT");     // Potentially add back
//	UncSources.push_back("SmEff_NTarget"); // Potentially add back
//	UncSources.push_back("SmEff_LY");
//	UncSources.push_back("SmEff_TPC");
//	UncSources.push_back("SmEff_XSec");
//	UncSources.push_back("SmEff_G4");
//	UncSources.push_back("SmEff_Flux");
//	UncSources.push_back("SmEff_Dirt");

	}

	if (IncludeUnfTech) { UncSources.push_back("UnfoldingTechnique"); }

	int NSamples = UncSources.size();

	vector<TFile*> CovFiles;
	CovFiles.resize(NSamples);

	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TH2D*> Covariances;
	Covariances.resize(NPlots);

	vector<TH2D*> StatCovariances;
	StatCovariances.resize(NPlots);

	vector<TH2D*> POTCovariances;
	POTCovariances.resize(NPlots);

	vector<TH2D*> NTargetCovariances;
	NTargetCovariances.resize(NPlots);

	vector<TH2D*> LYCovariances;
	LYCovariances.resize(NPlots);

	vector<TH2D*> TPCCovariances;
	TPCCovariances.resize(NPlots);

	vector<TH2D*> XSecCovariances;
	XSecCovariances.resize(NPlots);

	vector<TH2D*> G4Covariances;
	G4Covariances.resize(NPlots);

	vector<TH2D*> FluxCovariances;
	FluxCovariances.resize(NPlots);

	vector<TH2D*> DirtCovariances;
	DirtCovariances.resize(NPlots);

	vector<TH2D*> SystCovariances;
	SystCovariances.resize(NPlots);	

	vector<TH2D*> ERCovariances;
	ERCovariances.resize(NPlots);

	vector<TH2D*> SmEffCovariances;
	SmEffCovariances.resize(NPlots);	

	// -------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {					

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString TotalFileCovarianceSpecName = BeamOn9 + "WienerSVD_Total_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TString TotalFileCovarianceName = MigrationMatrixPath + TotalFileCovarianceSpecName;
		TFile* TotalFileCovarianceMatrices = new TFile(TotalFileCovarianceName,"recreate");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

//			TCanvas* DataERPlotCanvas = nullptr;
			TCanvas* MCERPlotCanvas = nullptr;
//			TCanvas* SmEffPlotCanvas = nullptr;

//			TLegend* legData = nullptr;
			TLegend* legMC = nullptr;
//			TLegend* legSmEff = nullptr;

//			TString DataERCanvasName = "DataERSyst_"+PlotNames[WhichPlot]+OverlaySample+"_"+Runs[WhichRun];
//			TString SmEffCanvasName = "SmEffSyst_"+PlotNames[WhichPlot]+OverlaySample+"_"+Runs[WhichRun];
			TString MCERCanvasName = "MCERSyst_"+PlotNames[WhichPlot]+OverlaySample+"_"+Runs[WhichRun];

			if (StorePlots) {

				// Create canvases for 3 categories of systematics: MC event rates, Data event rates, Smearing + Efficiency

//				DataERPlotCanvas = new TCanvas(DataERCanvasName,DataERCanvasName,205,34,1024,768);
//				DataERPlotCanvas->SetBottomMargin(0.16);
//				DataERPlotCanvas->SetLeftMargin(0.15);
//				DataERPlotCanvas->SetRightMargin(0.25);
//				DataERPlotCanvas->SetTopMargin(0.15);			

//				legData = new TLegend(0.02,0.89,0.97,0.99);
//				legData->SetBorderSize(0);
//				legData->SetTextSize(0.04);
//				legData->SetTextFont(FontStyle);
//				legData->SetNColumns(5);			

				MCERPlotCanvas = new TCanvas(MCERCanvasName,MCERCanvasName,205,34,1024,768);
				MCERPlotCanvas->SetBottomMargin(0.16);
				MCERPlotCanvas->SetLeftMargin(0.15);
				MCERPlotCanvas->SetRightMargin(0.25);
				MCERPlotCanvas->SetTopMargin(0.15);			

				legMC = new TLegend(0.02,0.89,0.97,0.99);
				legMC->SetBorderSize(0);
				legMC->SetTextSize(0.04);
				legMC->SetTextFont(FontStyle);
				legMC->SetNColumns(5);			

//				SmEffPlotCanvas = new TCanvas(SmEffCanvasName,SmEffCanvasName,205,34,1024,768);			
//				SmEffPlotCanvas->SetBottomMargin(0.16);
//				SmEffPlotCanvas->SetLeftMargin(0.15);
//				SmEffPlotCanvas->SetRightMargin(0.25);
//				SmEffPlotCanvas->SetTopMargin(0.15);			

//				legSmEff = new TLegend(0.02,0.89,0.97,0.99);
//				legSmEff->SetBorderSize(0);
//				legSmEff->SetTextSize(0.04);
//				legSmEff->SetTextFont(FontStyle);
//				legSmEff->SetNColumns(5);

			}

			// ----------------------------------------------------------------------------------------------------

			// Loop over the samples

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				// -----------------------------------------------------------------------------------------------------------------------------------------

				// Opening the file containing the covariance matrices for each one of the systematics

				TString FileCovarianceSpecName = BeamOn9+"WienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";

				// For the detector variation, we follow the PeLEE recipe
				// Only Run3 and propagate across all runs

				if (UncSources[WhichSample] == "LY" || UncSources[WhichSample] == "TPC" || UncSources[WhichSample] == "MC_LY" 
				|| UncSources[WhichSample] == "MC_TPC" || UncSources[WhichSample] == "SmEff_LY" || UncSources[WhichSample] == "SmEff_TPC" ) {

					FileCovarianceSpecName = BeamOn9+"WienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_Run3_"+UBCodeVersion+".root";

				}

				TString FileCovarianceName = MigrationMatrixPath + FileCovarianceSpecName;
				CovFiles[WhichSample] = new TFile(FileCovarianceName,"readonly");

				// -----------------------------------------------------------------------------------------------------------------------------------------

				TString LocalCovMatrixName = UncSources[WhichSample]+"_Covariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun];

				// For the detector variation, we follow the PeLEE recipe
				// Only Run3 and propagate across all runs

				if (UncSources[WhichSample] == "LY" || UncSources[WhichSample] == "TPC" || UncSources[WhichSample] == "MC_LY" 
				|| UncSources[WhichSample] == "MC_TPC" || UncSources[WhichSample] == "SmEff_LY" || UncSources[WhichSample] == "SmEff_TPC" ) {

					LocalCovMatrixName = UncSources[WhichSample]+"_Covariance_"+PlotNames[WhichPlot]+"_Run3";

				}

				TH2D* LocalCovMatrix = (TH2D*)( CovFiles[WhichSample]->Get(LocalCovMatrixName) );

				// -------------------------------------------------------------------------------------------------

				// Descriminate between Smearing / Efficiency covariances (SmEff) and MC Event Rates (ER)

				if (string(UncSources[WhichSample]).find("SmEff") != std::string::npos) { 

					if (SmEffCovariances[WhichPlot] == nullptr) { SmEffCovariances[WhichPlot] = LocalCovMatrix; }
					else { SmEffCovariances[WhichPlot]->Add(LocalCovMatrix); }

				} else {

					if (ERCovariances[WhichPlot] == nullptr) { ERCovariances[WhichPlot] = LocalCovMatrix; }
					else { ERCovariances[WhichPlot]->Add(LocalCovMatrix); }

				}

				// -------------------------------------------------------------------------------------------------


				if (string(UncSources[WhichSample]).find("Stat") != std::string::npos) { 

					if (StatCovariances[WhichPlot] == nullptr) { StatCovariances[WhichPlot] = LocalCovMatrix; }
					else { StatCovariances[WhichPlot]->Add(LocalCovMatrix); }
				
				} else {

					// --------------------------------------------------------------------------------

					if (string(UncSources[WhichSample]).find("POT") != std::string::npos) { 

						if (POTCovariances[WhichPlot] == nullptr) { POTCovariances[WhichPlot] = LocalCovMatrix; }
						else { POTCovariances[WhichPlot]->Add(LocalCovMatrix); }

					}

					if (string(UncSources[WhichSample]).find("NTarget") != std::string::npos) { 

						if (NTargetCovariances[WhichPlot] == nullptr) { NTargetCovariances[WhichPlot] = LocalCovMatrix; }
						else { NTargetCovariances[WhichPlot]->Add(LocalCovMatrix); }

					}

					if (string(UncSources[WhichSample]).find("LY") != std::string::npos) { 

						if (LYCovariances[WhichPlot] == nullptr) { LYCovariances[WhichPlot] = LocalCovMatrix; }
						else { LYCovariances[WhichPlot]->Add(LocalCovMatrix); }

					}

					if (string(UncSources[WhichSample]).find("TPC") != std::string::npos) { 

						if (TPCCovariances[WhichPlot] == nullptr) { TPCCovariances[WhichPlot] = LocalCovMatrix; }
						else { TPCCovariances[WhichPlot]->Add(LocalCovMatrix); }

					}

					if (string(UncSources[WhichSample]).find("XSec") != std::string::npos) { 

						if (XSecCovariances[WhichPlot] == nullptr) { XSecCovariances[WhichPlot] = LocalCovMatrix; }
						else { XSecCovariances[WhichPlot]->Add(LocalCovMatrix); }

					}

					if (string(UncSources[WhichSample]).find("G4") != std::string::npos) { 

						if (G4Covariances[WhichPlot] == nullptr) { G4Covariances[WhichPlot] = LocalCovMatrix; }
						else { G4Covariances[WhichPlot]->Add(LocalCovMatrix); }

					}

					if (string(UncSources[WhichSample]).find("Flux") != std::string::npos) { 

						if (FluxCovariances[WhichPlot] == nullptr) { FluxCovariances[WhichPlot] = LocalCovMatrix; }
						else { FluxCovariances[WhichPlot]->Add(LocalCovMatrix); }

					}

					if (string(UncSources[WhichSample]).find("Dirt") != std::string::npos) { 

						if (DirtCovariances[WhichPlot] == nullptr) { DirtCovariances[WhichPlot] = LocalCovMatrix; }
						else { DirtCovariances[WhichPlot]->Add(LocalCovMatrix); }

					}

					// --------------------------------------------------------------------------------

					if ( !(string(UncSources[WhichSample]).find("Stat") != std::string::npos) ) { 

						if (SystCovariances[WhichPlot] == nullptr) { SystCovariances[WhichPlot] = LocalCovMatrix; }
						else { SystCovariances[WhichPlot]->Add(LocalCovMatrix); }
						//TotalFileCovarianceMatrices->cd();
						//SystCovariances[WhichPlot]->Write(UncSources[WhichSample]+"Covariance_"+PlotNames[WhichPlot]);

					}
//					else if (WhichSample > 2) { 
//						
//						TotalFileCovarianceMatrices->cd();
//						//LocalCovMatrix->Write(UncSources[WhichSample]+"Covariance_"+PlotNames[WhichPlot]);
//						SystCovariances[WhichPlot]->Add(LocalCovMatrix); 
//						
//					}
					
				}

				// --------------------------------------------------------------------------------

				if (StorePlots) {

					TLegend* leg = nullptr;

					MCERPlotCanvas->cd();  leg = legMC;

//					if (string(UncSources[WhichSample]).find("SmEff_") != std::string::npos)
//						{ SmEffPlotCanvas->cd(); leg = legSmEff;}
//					else { MCERPlotCanvas->cd();  leg = legMC; }				
//					else if (string(UncSources[WhichSample]).find("MC_") != std::string::npos) { MCERPlotCanvas->cd();  leg = legMC; }				
//					else { DataERPlotCanvas->cd(); leg = legData; }
//					else { cout << "WARNING !!! No canvas to point to !"; }

					ReturnUncPlot(LocalCovMatrix,PlotNames[WhichPlot],Runs[WhichRun],UncSources[WhichSample],WhichSample,leg);				

				}

			} // End of the loop over the samples

			// ------------------------------------------------------------------

			if (StorePlots) {

//				DataERPlotCanvas->cd();
//				legData->Draw();

				MCERPlotCanvas->cd();
				legMC->Draw();

//				SmEffPlotCanvas->cd();
//				legSmEff->Draw();
		
				// ------------------------------------------------------------------

//				DataERPlotCanvas->SaveAs(PlotPath+OverlaySample+"/"+DataERCanvasName+".pdf");
				MCERPlotCanvas->SaveAs(PlotPath+OverlaySample+"/"+MCERCanvasName+".pdf");
//				SmEffPlotCanvas->SaveAs(PlotPath+OverlaySample+"/"+SmEffCanvasName+".pdf");	

//				delete 	DataERPlotCanvas;	
				delete 	MCERPlotCanvas;
//				delete 	SmEffPlotCanvas;

			}

			// ------------------------------------------------------------------			

			// Storing Total Covariance Matrices in root file

			TotalFileCovarianceMatrices->cd();
			StatCovariances[WhichPlot]->Write("StatCovariance_"+PlotNames[WhichPlot]);
			SystCovariances[WhichPlot]->Write("SystCovariance_"+PlotNames[WhichPlot]);						

			Covariances[WhichPlot] = (TH2D*) (StatCovariances[WhichPlot]);
			Covariances[WhichPlot]->Add(SystCovariances[WhichPlot]);			
			Covariances[WhichPlot]->Write("TotalCovariance_"+PlotNames[WhichPlot]);
			
			// ---------------------------------------------------------------------------------------------
			
			if (StorePlots) {

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

			}
		
		} // End of the loop over the plots

		TotalFileCovarianceMatrices->Close();
		cout << endl << "Covariance matrix file " << TotalFileCovarianceName << " has been created" << endl << endl;

		cout << "Merging of covariance matrices for run " << Runs[WhichRun] << " completed!" << endl;

	} // End of the loop over the runs	

} // End of the program 
