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

void PlotCov(TH2D* h, TString Label, TString PlotNames, TString OverlaySamples, TString Runs) {

	TString CanvasName = "Total_"+PlotNames+OverlaySamples+"_"+Runs;
	TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
	PlotCanvas->cd();
	PlotCanvas->SetBottomMargin(0.16);
	PlotCanvas->SetLeftMargin(0.15);
	PlotCanvas->SetRightMargin(0.25);			
	
	gStyle->SetMarkerSize(1.5);
	gStyle->SetPaintTextFormat("4.3f");			
	
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleSize(TextSize);
	h->GetXaxis()->SetLabelSize(TextSize);			
	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetNdivisions(5);
	
	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetTitleSize(TextSize);
	h->GetYaxis()->SetLabelSize(TextSize);			
	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetNdivisions(5);
	h->GetYaxis()->SetTitleOffset(1.);			

	h->SetTitle(Runs + Label + " Total");	

	double FracCovMax = TMath::Min(1.,1.05 * h->GetMaximum());
	double FracCovMin = TMath::Min(0.,1.05 * h->GetMinimum());
	h->GetZaxis()->SetRangeUser(FracCovMin,FracCovMax);
	h->GetZaxis()->CenterTitle();
	h->GetZaxis()->SetTitleFont(FontStyle);
	h->GetZaxis()->SetTitleSize(TextSize);
	h->GetZaxis()->SetLabelFont(FontStyle);
	h->GetZaxis()->SetLabelSize(TextSize-0.01);
	h->GetZaxis()->SetNdivisions(5);

	h->SetMarkerColor(kWhite);			
	h->SetMarkerSize(1.5);
	//h->Draw("text colz e"); 
	h->Draw("colz");
	
	PlotCanvas->SaveAs(PlotPath+OverlaySamples+"/WienerSVD_Total_"+Label+"CovarianceMatrices_"+PlotNames+OverlaySamples+"_"+Runs+"_"+UBCodeVersion+".pdf");
	
	delete PlotCanvas;

}

// -----------------------------------------------------------------------------------------------

void ReturnUncPlot(TH2D* LocalCovMatrix,TString PlotName, TString Run,TString UncSources,int Color, TLegend* leg) {

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();

	//if (string(UncSources).find("POT") != std::string::npos || string(UncSources).find("NTarget") != std::string::npos) { return; }
	//else {

		std::vector<int> Colors{kBlack, kRed+1,kGreen+2,kOrange+1,kBlue,kMagenta, kViolet+1, kYellow + 2, kCyan+1, kGreen, kRed, kGray-1};

		int n = LocalCovMatrix->GetNbinsX();
		TString TitleX =  LocalCovMatrix->GetXaxis()->GetTitle();
		double Nuedges[n+1];
				    
		for (int i = 0; i < n+1; i++) { Nuedges[i] = LocalCovMatrix->GetXaxis()->GetBinLowEdge(i+1); }

		TH1D* unc = new TH1D("unc_"+PlotName+"_"+Run+"_"+UncSources,";"+TitleX+";Uncertainty [%]",n,Nuedges);	

		double SumBins = 0.;
		double SumWeightBins = 0.;

		for (int i = 1; i <= n; i++) { 

			double CovValue = LocalCovMatrix->GetBinContent(i,i);	
			unc->SetBinContent(i,TMath::Sqrt(CovValue)*100.);

			double BinWidth = unc->GetBinWidth(i);
			double BinCont = unc->GetBinContent(i);

			SumBins += BinWidth;
			SumWeightBins += BinCont*BinWidth;
			
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
		unc->GetYaxis()->SetRangeUser(0.,20.);

		if (PlotName == "DeltaPTPlot") { unc->GetYaxis()->SetRangeUser(0.,49.); }
//		if (PlotName == "DeltaAlphaTPlot") { unc->GetYaxis()->SetRangeUser(0.,14.); }
		if (PlotName == "DeltaPhiTPlot") { unc->GetYaxis()->SetRangeUser(0.,47.); }
		if (PlotName == "MuonMomentumPlot") { unc->GetYaxis()->SetRangeUser(0.,37.); }
		if (PlotName == "MuonPhiPlot") { unc->GetYaxis()->SetRangeUser(0.,30.); }
		if (PlotName == "MuonCosThetaPlot") { unc->GetYaxis()->SetRangeUser(0.,49.); }	
		if (PlotName == "ProtonPhiPlot") { unc->GetYaxis()->SetRangeUser(0.,30.); }
		if (PlotName == "ProtonMomentumPlot") { unc->GetYaxis()->SetRangeUser(0.,39.); }	
		if (PlotName == "ProtonCosThetaPlot") { unc->GetYaxis()->SetRangeUser(0.,49.); }	
					
		unc->SetLineWidth(2);		
		unc->SetLineColor(Colors[Color+1]);
		unc->SetMarkerColor(Colors[Color+1]);			

		unc->Draw("hist text0 same");

		leg->AddEntry(unc,UncSources);

		// ----------------------------------------------------------------

}

// -----------------------------------------------------------------------------------------------

void WienerSVD_Merge_Covariances(TString OverlaySample = "Overlay9", TString BeamOn9 = "") {

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
//	PlotNames.push_back("DeltaPTPlot"); 
//	PlotNames.push_back("DeltaAlphaTPlot"); 
//	PlotNames.push_back("DeltaPhiTPlot");
//	PlotNames.push_back("MuonMomentumPlot");
//	PlotNames.push_back("MuonCosThetaPlot"); 
//	PlotNames.push_back("MuonCosThetaSingleBinPlot"); 
//	PlotNames.push_back("MuonPhiPlot");
//	PlotNames.push_back("ProtonMomentumPlot"); 
//	PlotNames.push_back("ProtonCosThetaPlot");
//	PlotNames.push_back("ProtonPhiPlot");
//	PlotNames.push_back("ECalPlot");
//	PlotNames.push_back("EQEPlot"); 
//	PlotNames.push_back("Q2Plot");

	PlotNames.push_back("CCQEMuonMomentumPlot"); 
	PlotNames.push_back("CCQEMuonCosThetaPlot"); 
	PlotNames.push_back("CCQEProtonMomentumPlot"); 
	PlotNames.push_back("CCQEProtonCosThetaPlot");
	PlotNames.push_back("CCQEECalPlot");
	PlotNames.push_back("CCQEQ2Plot");

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

	if (BeamOn9 != "") { // Fake data study, we need only the Stat & XSec covariances

		UncSources.push_back("Stat");
		UncSources.push_back("XSec");
		UncSources.push_back("MC_Stat");
		UncSources.push_back("NuWro");			

	} else {

		UncSources.push_back("Stat");
		UncSources.push_back("LY");
		UncSources.push_back("TPC");
		UncSources.push_back("SCERecomb2");
		UncSources.push_back("XSec");
		UncSources.push_back("G4");
		UncSources.push_back("Flux");
		UncSources.push_back("Dirt");
		UncSources.push_back("POT"); 
		UncSources.push_back("NTarget");
		UncSources.push_back("MC_Stat");
		UncSources.push_back("NuWro");		

	}

	int NSamples = UncSources.size();

	vector<TFile*> CovFiles;
	CovFiles.resize(NSamples);

	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TH2D*> FracCovariances;
	FracCovariances.resize(NPlots);

	vector<TH2D*> FracStatCovariances;
	FracStatCovariances.resize(NPlots);

	vector<TH2D*> FracSystCovariances;
	FracSystCovariances.resize(NPlots);

	vector<TH2D*> Covariances;
	Covariances.resize(NPlots);

	vector<TH2D*> StatCovariances;
	StatCovariances.resize(NPlots);

	vector<TH2D*> SystCovariances;
	SystCovariances.resize(NPlots);

//	vector<TH2D*> POTCovariances;
//	POTCovariances.resize(NPlots);

//	vector<TH2D*> NTargetCovariances;
//	NTargetCovariances.resize(NPlots);

//	vector<TH2D*> LYCovariances;
//	LYCovariances.resize(NPlots);

//	vector<TH2D*> TPCCovariances;
//	TPCCovariances.resize(NPlots);

//	vector<TH2D*> XSecCovariances;
//	XSecCovariances.resize(NPlots);

//	vector<TH2D*> G4Covariances;
//	G4Covariances.resize(NPlots);

//	vector<TH2D*> FluxCovariances;
//	FluxCovariances.resize(NPlots);

//	vector<TH2D*> DirtCovariances;
//	DirtCovariances.resize(NPlots);	

//	vector<TH2D*> ERCovariances;
//	ERCovariances.resize(NPlots);

//	vector<TH2D*> SmEffCovariances;
//	SmEffCovariances.resize(NPlots);	

	// -------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {					

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString TotalFileCovarianceSpecName = BeamOn9 + "CCQEWienerSVD_Total_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TString TotalFileCovarianceName = MigrationMatrixPath + TotalFileCovarianceSpecName;
		TFile* TotalFileCovarianceMatrices = new TFile(TotalFileCovarianceName,"recreate");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

			TCanvas* MCERPlotCanvas = nullptr;
			TLegend* legMC = nullptr;
			TString MCERCanvasName = "CCQEMCERSyst_"+PlotNames[WhichPlot]+OverlaySample+"_"+Runs[WhichRun];

			if (StorePlots) {

				// Create canvases for MC event rates			

				MCERPlotCanvas = new TCanvas(MCERCanvasName,MCERCanvasName,205,34,1024,768);
				MCERPlotCanvas->SetBottomMargin(0.16);
				MCERPlotCanvas->SetLeftMargin(0.15);
				MCERPlotCanvas->SetRightMargin(0.25);
				MCERPlotCanvas->SetTopMargin(0.15);			

				legMC = new TLegend(0.15,0.89,0.9,0.99);
				legMC->SetBorderSize(0);
				legMC->SetTextSize(0.04);
				legMC->SetTextFont(FontStyle);
				legMC->SetNColumns(5);			

			}

			// ----------------------------------------------------------------------------------------------------

			// Loop over the samples

			int Counter = 0; // Counter for colors

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				// -----------------------------------------------------------------------------------------------------------------------------------------

				// Opening the file containing the covariance matrices for each one of the systematics

				TString FileCovarianceSpecName = BeamOn9+"CCQEWienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";

				// For the detector variation, we follow the PeLEE recipe
				// Only Run3 and propagate across all runs

				if (UncSources[WhichSample] == "LY" || UncSources[WhichSample] == "TPC"  || UncSources[WhichSample] == "SCERecomb2" || UncSources[WhichSample] == "MC_LY" 
				|| UncSources[WhichSample] == "MC_TPC" || UncSources[WhichSample] == "MC_SCERecomb2" || UncSources[WhichSample] == "SmEff_LY" 
				|| UncSources[WhichSample] == "SmEff_TPC" || UncSources[WhichSample] == "SmEff_SCERecomb2" ) {

					FileCovarianceSpecName = BeamOn9+"CCQEWienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_Run3_"+UBCodeVersion+".root";

				}

				TString FileCovarianceName = MigrationMatrixPath + FileCovarianceSpecName;
				CovFiles[WhichSample] = new TFile(FileCovarianceName,"readonly");

				// -----------------------------------------------------------------------------------------------------------------------------------------

				TString LocalCovMatrixName = UncSources[WhichSample]+"_Covariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun];
				TString LocalFracCovMatrixName = UncSources[WhichSample]+"_FracCovariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun];

				// For the detector variation, we follow the PeLEE recipe
				// Only Run3 and propagate across all runs

				if (UncSources[WhichSample] == "LY" || UncSources[WhichSample] == "TPC" || UncSources[WhichSample] == "SCERecomb2" || UncSources[WhichSample] == "MC_LY" 
				|| UncSources[WhichSample] == "MC_TPC" || UncSources[WhichSample] == "MC_SCERecomb2" || UncSources[WhichSample] == "SmEff_LY" 
				|| UncSources[WhichSample] == "SmEff_TPC" || UncSources[WhichSample] == "SmEff_SCERecomb2" ) {

					LocalCovMatrixName = UncSources[WhichSample]+"_Covariance_"+PlotNames[WhichPlot]+"_Run3";
					LocalFracCovMatrixName = UncSources[WhichSample]+"_FracCovariance_"+PlotNames[WhichPlot]+"_Run3";

				}

				TH2D* LocalCovMatrix = (TH2D*)( CovFiles[WhichSample]->Get(LocalCovMatrixName) );
				TH2D* LocalFracCovMatrix = (TH2D*)( CovFiles[WhichSample]->Get(LocalFracCovMatrixName) );

				// -------------------------------------------------------------------------------------------------

				// Merging all the covariance matrices into one

				if (WhichSample == 0) { Covariances[WhichPlot] = LocalCovMatrix; FracCovariances[WhichPlot] = LocalFracCovMatrix; }
				else { Covariances[WhichPlot]->Add(LocalCovMatrix); FracCovariances[WhichPlot]->Add(LocalFracCovMatrix); }

				// -------------------------------------------------------------------------------------------------

				if (string(UncSources[WhichSample]).find("Stat") != std::string::npos) { 

					if (string(UncSources[WhichSample]).find("Stat") != std::string::npos) 
						{ StatCovariances[WhichPlot] = LocalCovMatrix; FracStatCovariances[WhichPlot] = LocalFracCovMatrix; }
					else if (string(UncSources[WhichSample]).find("MC_Stat") != std::string::npos) 
						{ StatCovariances[WhichPlot]->Add(LocalCovMatrix); FracStatCovariances[WhichPlot]->Add(LocalFracCovMatrix); }
					else { cout << "What is this stat covariance matrix ?" << endl; }
				
				} else {

					if (WhichSample == 1) { SystCovariances[WhichPlot] = LocalCovMatrix; FracSystCovariances[WhichPlot] = LocalFracCovMatrix; }
					else { SystCovariances[WhichPlot]->Add(LocalCovMatrix); FracSystCovariances[WhichPlot]->Add(LocalFracCovMatrix); }

				}

				// -------------------------------------------------------------------------------------------------

//				// Descriminate between Smearing / Efficiency covariances (SmEff) and MC Event Rates (ER)

//				if (string(UncSources[WhichSample]).find("SmEff") != std::string::npos) { 

//					if (SmEffCovariances[WhichPlot] == nullptr) { SmEffCovariances[WhichPlot] = LocalCovMatrix; }
//					else { SmEffCovariances[WhichPlot]->Add(LocalCovMatrix); }

//				} else {

//					if (ERCovariances[WhichPlot] == nullptr) { ERCovariances[WhichPlot] = LocalCovMatrix; }
//					else { ERCovariances[WhichPlot]->Add(LocalCovMatrix); }

//				}

				// -------------------------------------------------------------------------------------------------


//				if (string(UncSources[WhichSample]).find("Stat") != std::string::npos) { 

//					if (StatCovariances[WhichPlot] == nullptr) { StatCovariances[WhichPlot] = LocalCovMatrix; }
//					else { StatCovariances[WhichPlot]->Add(LocalCovMatrix); }
//				
//				} else {

//					// --------------------------------------------------------------------------------

//					if (string(UncSources[WhichSample]).find("POT") != std::string::npos) { 

//						if (POTCovariances[WhichPlot] == nullptr) { POTCovariances[WhichPlot] = LocalCovMatrix; }
//						else { POTCovariances[WhichPlot]->Add(LocalCovMatrix); }

//					}

//					if (string(UncSources[WhichSample]).find("NTarget") != std::string::npos) { 

//						if (NTargetCovariances[WhichPlot] == nullptr) { NTargetCovariances[WhichPlot] = LocalCovMatrix; }
//						else { NTargetCovariances[WhichPlot]->Add(LocalCovMatrix); }

//					}

//					if (string(UncSources[WhichSample]).find("LY") != std::string::npos) { 

//						if (LYCovariances[WhichPlot] == nullptr) { LYCovariances[WhichPlot] = LocalCovMatrix; }
//						else { LYCovariances[WhichPlot]->Add(LocalCovMatrix); }

//					}

//					if (string(UncSources[WhichSample]).find("TPC") != std::string::npos) { 

//						if (TPCCovariances[WhichPlot] == nullptr) { TPCCovariances[WhichPlot] = LocalCovMatrix; }
//						else { TPCCovariances[WhichPlot]->Add(LocalCovMatrix); }

//					}

//					if (string(UncSources[WhichSample]).find("XSec") != std::string::npos) { 

//						if (XSecCovariances[WhichPlot] == nullptr) { XSecCovariances[WhichPlot] = LocalCovMatrix; }
//						else { XSecCovariances[WhichPlot]->Add(LocalCovMatrix); }

//					}

//					if (string(UncSources[WhichSample]).find("G4") != std::string::npos) { 

//						if (G4Covariances[WhichPlot] == nullptr) { G4Covariances[WhichPlot] = LocalCovMatrix; }
//						else { G4Covariances[WhichPlot]->Add(LocalCovMatrix); }

//					}

//					if (string(UncSources[WhichSample]).find("Flux") != std::string::npos) { 

//						if (FluxCovariances[WhichPlot] == nullptr) { FluxCovariances[WhichPlot] = LocalCovMatrix; }
//						else { FluxCovariances[WhichPlot]->Add(LocalCovMatrix); }

//					}

//					if (string(UncSources[WhichSample]).find("Dirt") != std::string::npos) { 

//						if (DirtCovariances[WhichPlot] == nullptr) { DirtCovariances[WhichPlot] = LocalCovMatrix; }
//						else { DirtCovariances[WhichPlot]->Add(LocalCovMatrix); }

//					}

//					// --------------------------------------------------------------------------------

//					if ( !(string(UncSources[WhichSample]).find("Stat") != std::string::npos) ) { 

//						if (SystCovariances[WhichPlot] == nullptr) { SystCovariances[WhichPlot] = LocalCovMatrix; }
//						else { SystCovariances[WhichPlot]->Add(LocalCovMatrix); }
//						//TotalFileCovarianceMatrices->cd();
//						//SystCovariances[WhichPlot]->Write(UncSources[WhichSample]+"Covariance_"+PlotNames[WhichPlot]);

//					}
////					else if (WhichSample > 2) { 
////						
////						TotalFileCovarianceMatrices->cd();
////						//LocalCovMatrix->Write(UncSources[WhichSample]+"Covariance_"+PlotNames[WhichPlot]);
////						SystCovariances[WhichPlot]->Add(LocalCovMatrix); 
////						
////					}
					
//				}

				// --------------------------------------------------------------------------------

				if (StorePlots) {

					TLegend* leg = nullptr;

					MCERPlotCanvas->cd();  leg = legMC;

					if ( !( string(UncSources[WhichSample]).find("POT") != std::string::npos || string(UncSources[WhichSample]).find("NTarget") != std::string::npos ) ) {

						ReturnUncPlot(LocalFracCovMatrix,PlotNames[WhichPlot],Runs[WhichRun],UncSources[WhichSample],Counter,leg);
						Counter++;

					}				

				}

			} // End of the loop over the samples

			// ------------------------------------------------------------------			

			// Storing Total Covariance Matrices in root file

			TotalFileCovarianceMatrices->cd();
			Covariances[WhichPlot]->Write("TotalCovariance_"+PlotNames[WhichPlot]);
			FracCovariances[WhichPlot]->Write("FracTotalCovariance_"+PlotNames[WhichPlot]);

			StatCovariances[WhichPlot]->Write("StatCovariance_"+PlotNames[WhichPlot]);
			FracStatCovariances[WhichPlot]->Write("FracStatCovariance_"+PlotNames[WhichPlot]);

			SystCovariances[WhichPlot]->Write("SystCovariance_"+PlotNames[WhichPlot]);
			FracSystCovariances[WhichPlot]->Write("FracSystCovariance_"+PlotNames[WhichPlot]);

			// ------------------------------------------------------------------

			if (StorePlots) {

				// ------------------------------------------------------------------

				MCERPlotCanvas->cd();
				legMC->Draw();

				// ------------------------------------------------------------------

				// Plot the total unc on top of everything else

				ReturnUncPlot(FracCovariances[WhichPlot],PlotNames[WhichPlot],Runs[WhichRun],"Total",-1,legMC);
		
				// ------------------------------------------------------------------

				MCERPlotCanvas->SaveAs(PlotPath+OverlaySample+"/"+MCERCanvasName+".pdf");
				delete 	MCERPlotCanvas;

				// ---------------------------------------------------------------------------------------------

				// Plot the 2D covariance matrices

				PlotCov(Covariances[WhichPlot],"",PlotNames[WhichPlot],OverlaySample,Runs[WhichRun]);

				// ---------------------------------------------------------------------------------------------

				// Plot the 2D fractional covariance matrices

				PlotCov(FracCovariances[WhichPlot],"Frac",PlotNames[WhichPlot],OverlaySample,Runs[WhichRun]);

				// ---------------------------------------------------------------------------------------------

			}
		
		} // End of the loop over the plots

		TotalFileCovarianceMatrices->Close();
		cout << endl << "Covariance matrix file " << TotalFileCovarianceName << " has been created" << endl << endl;

		cout << "Merging of covariance matrices for run " << Runs[WhichRun] << " completed!" << endl;

	} // End of the loop over the runs	

} // End of the program 
