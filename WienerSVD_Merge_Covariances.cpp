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
#include <TGaxis.h>
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
#include "ubana/AnalysisCode/Secondary_Code/myFunctions.cpp"

// -----------------------------------------------------------------------------------------------

void PlotCov(TH2D* h, TString Label, TString PlotNames, TString OverlaySamples, TString Runs, TString Tune = "") {

	TString CanvasName = "Total_"+PlotNames+OverlaySamples+"_"+Runs;
	TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
	PlotCanvas->cd();
	PlotCanvas->SetBottomMargin(0.16);
	PlotCanvas->SetLeftMargin(0.15);
	PlotCanvas->SetRightMargin(0.25);			
	
	gStyle->SetMarkerSize(1.5);
	gStyle->SetPaintTextFormat("4.2f");			
	
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

	double FracCovMax = FindTwoDimHistoMaxValue(h);
	double FracCovMin = FindTwoDimHistoMinValue(h);

//	double FracCovMax = TMath::Min(1.,1.05 * h->GetMaximum());
//	double FracCovMin = TMath::Min(0.,1.05 * h->GetMinimum());	

	TString Title = "Title";
	if (Label == "") { Title = "Cov Matrix"; }
	if (Label == "Frac") { Title = "Frac Cov Matrix"; }
	if (Label == "Corr") { Title = "Corr Matrix"; }		

	h->SetTitle("Total " + Title);	

	h->GetZaxis()->SetRangeUser(FracCovMin,FracCovMax);
	h->GetZaxis()->CenterTitle();
	h->GetZaxis()->SetTitleFont(FontStyle);
	h->GetZaxis()->SetTitleSize(TextSize);
	h->GetZaxis()->SetLabelFont(FontStyle);
	h->GetZaxis()->SetLabelSize(TextSize-0.01);
	h->GetZaxis()->SetNdivisions(5);

	h->SetMarkerColor(kWhite);			
	h->SetMarkerSize(1.);
	//h->Draw("text colz e"); 

	if (string(PlotNames).find("Serial") != std::string::npos) {	

		TString XaxisTitle = h->GetXaxis()->GetTitle();
		XaxisTitle.ReplaceAll("deg","bin #");
		XaxisTitle.ReplaceAll("GeV/c","bin #");
		XaxisTitle.ReplaceAll("GeV","bin #");				
		h->GetXaxis()->SetTitle(XaxisTitle);

		TString YaxisTitle = h->GetYaxis()->GetTitle();
		YaxisTitle.ReplaceAll("deg","bin #");
		YaxisTitle.ReplaceAll("GeV/c","bin #");
		YaxisTitle.ReplaceAll("GeV","bin #");				
		h->GetYaxis()->SetTitle(YaxisTitle);									

	}

	if (Label == "Corr") { h->Draw("colz text"); }
	else { h->Draw("colz"); }
	
	PlotCanvas->SaveAs(PlotPath+OverlaySamples+"/"+Tune+"WienerSVD_Total_"+Label+"CovarianceMatrices_"+PlotNames+OverlaySamples+"_"+Runs+"_"+UBCodeVersion+".pdf");
	
	delete PlotCanvas;

}

// -----------------------------------------------------------------------------------------------

void ReturnUncPlot(TH2D* LocalCovMatrix,TString PlotName, TString Run,TString UncSources,int Color, TLegend* leg) {

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();

	//if (string(UncSources).find("POT") != std::string::npos || string(UncSources).find("NTarget") != std::string::npos) { return; }
	//else {

		std::vector<int> Colors{kBlack, kRed+1,kGreen+2,kOrange+1,kBlue,kMagenta, kViolet+1, kYellow + 2, kCyan+1, kGreen, kRed-10, kGray,kMagenta-10};

		int n = LocalCovMatrix->GetNbinsX();
		TString TitleX =  LocalCovMatrix->GetXaxis()->GetTitle();
		double Nuedges[n+1];
				    
		for (int i = 0; i < n+1; i++) { Nuedges[i] = LocalCovMatrix->GetXaxis()->GetBinLowEdge(i+1); }

		TH1D* unc = new TH1D("unc_"+PlotName+"_"+Run+"_"+UncSources,";"+TitleX+";Uncertainty [%]",n,Nuedges);	

        for (int i = 1; i <= n; i++) {

            double CovValue = LocalCovMatrix->GetBinContent(i,i);
            unc->SetBinContent(i,TMath::Sqrt(CovValue)*100.);

        }

		unc->GetXaxis()->SetTitleFont(FontStyle);
		unc->GetXaxis()->SetLabelFont(FontStyle);
		unc->GetXaxis()->SetTitleSize(TextSize);
		unc->GetXaxis()->SetLabelSize(TextSize);	
		unc->GetXaxis()->CenterTitle();
		unc->GetXaxis()->SetNdivisions(8);
		
		unc->GetYaxis()->SetLabelFont(FontStyle);
		unc->GetYaxis()->SetTitleFont(FontStyle);
		unc->GetYaxis()->SetTitleSize(TextSize);
		unc->GetYaxis()->SetLabelSize(TextSize);	
		unc->GetYaxis()->CenterTitle();
		unc->GetYaxis()->SetNdivisions(5);
		unc->GetYaxis()->SetTitleOffset(0.7);				
		unc->GetYaxis()->SetRangeUser(0.,34.);

		if (PlotName == "DeltaPTPlot") { unc->GetYaxis()->SetRangeUser(0.,46.); }
		if (PlotName == "DeltaAlphaTPlot") { unc->GetYaxis()->SetRangeUser(0.,18.); }
		if (PlotName == "DeltaPhiTPlot") { unc->GetYaxis()->SetRangeUser(0.,34.); }
		if (PlotName == "MuonMomentumPlot") { unc->GetYaxis()->SetRangeUser(0.,34.); }
		if (PlotName == "MuonPhiPlot") { unc->GetYaxis()->SetRangeUser(0.,19.); }
		if (PlotName == "MuonCosThetaPlot") { unc->GetYaxis()->SetRangeUser(0.,39.); }			
		if (PlotName == "ProtonPhiPlot") { unc->GetYaxis()->SetRangeUser(0.,19.); }
		if (PlotName == "ProtonMomentumPlot") { unc->GetYaxis()->SetRangeUser(0.,29.); }	
		if (PlotName == "ProtonCosThetaPlot") { unc->GetYaxis()->SetRangeUser(0.,34.); }	
		if (PlotName == "DeltaPLPlot") { unc->GetYaxis()->SetRangeUser(0.,39.); }
		if (PlotName == "DeltaPnPlot") { unc->GetYaxis()->SetRangeUser(0.,43.); }
		if (PlotName == "DeltaPtxPlot") { unc->GetYaxis()->SetRangeUser(0.,39.); }		
		if (PlotName == "DeltaPtyPlot") { unc->GetYaxis()->SetRangeUser(0.,39.); }	
		if (PlotName == "PMissMinusPlot") { unc->GetYaxis()->SetRangeUser(0.,49.); }
		if (PlotName == "ECalPlot") { unc->GetYaxis()->SetRangeUser(0.,109.); }	

		if (PlotName == "MuonCosThetaSingleBinPlot") { 

			unc->GetXaxis()->SetTitleSize(0.);
			unc->GetXaxis()->SetLabelSize(0.);			

			unc->GetYaxis()->SetRangeUser(0.,12.); 
			
		}

		if ( string(PlotName).find("Serial") != std::string::npos) {

			unc->GetYaxis()->SetRangeUser(0.,119.);

		}							
					
		unc->SetLineWidth(2);		
		unc->SetLineColor(Colors[Color+1]);
		unc->SetMarkerColor(Colors[Color+1]);			

		gStyle->SetPaintTextFormat("4.1f");
		unc->Draw("hist text0 same");

		leg->AddEntry(unc,UncSources);

		// ----------------------------------------------------------------

}

// -----------------------------------------------------------------------------------------------

void WienerSVD_Merge_Covariances(TString OverlaySample = "Overlay9", TString BeamOn9 = "", TString Tune = "") {

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

//	vector<TString> PlotNames;
//	PlotNames.push_back("DeltaPTPlot"); 
//	PlotNames.push_back("DeltaAlphaTPlot"); 

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

	if (BeamOn9 != "") { // Fake data study, we need only the Stat/MC_Stat & XSec covariances

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

	vector<TH2D*> FracCovariances; FracCovariances.resize(NPlots);
	vector<TH2D*> FracStatCovariances; FracStatCovariances.resize(NPlots);
	vector<TH2D*> FracSystCovariances; FracSystCovariances.resize(NPlots);

	vector<TH2D*> Covariances; Covariances.resize(NPlots);
	vector<TH2D*> StatCovariances; StatCovariances.resize(NPlots);
	vector<TH2D*> MCStatCovariances; MCStatCovariances.resize(NPlots);	
	vector<TH2D*> SystCovariances; SystCovariances.resize(NPlots);	
	vector<TH2D*> LYCovariances; LYCovariances.resize(NPlots);	
	vector<TH2D*> TPCCovariances; TPCCovariances.resize(NPlots);
	vector<TH2D*> SCERecomb2Covariances; SCERecomb2Covariances.resize(NPlots);
	vector<TH2D*> XSecCovariances; XSecCovariances.resize(NPlots);
	vector<TH2D*> G4Covariances; G4Covariances.resize(NPlots);
	vector<TH2D*> FluxCovariances; FluxCovariances.resize(NPlots);
	vector<TH2D*> DirtCovariances; DirtCovariances.resize(NPlots);	
	vector<TH2D*> POTCovariances; POTCovariances.resize(NPlots);
	vector<TH2D*> NuWroCovariances; NuWroCovariances.resize(NPlots);
	vector<TH2D*> NTargetCovariances; NTargetCovariances.resize(NPlots);							

	// -------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {					

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString TotalFileCovarianceSpecName = Tune + BeamOn9 + "WienerSVD_Total_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		if (BeamOn9 != "") { TotalFileCovarianceSpecName = Tune + "WienerSVD_Total_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }

		TString TotalFileCovarianceName = MigrationMatrixPath + TotalFileCovarianceSpecName;
		TFile* TotalFileCovarianceMatrices = new TFile(TotalFileCovarianceName,"recreate");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

			TCanvas* MCERPlotCanvas = nullptr;
			TLegend* legMC = nullptr;
			TString MCERCanvasName = Tune + "MCERSyst_"+PlotNames[WhichPlot]+OverlaySample+"_"+Runs[WhichRun];

			if (StorePlots) {

				// Create canvases for MC event rates			

//				MCERPlotCanvas = new TCanvas(MCERCanvasName,MCERCanvasName,205,34,1024,768);
				MCERPlotCanvas = new TCanvas(MCERCanvasName,MCERCanvasName,205,34,2000,1000);
				MCERPlotCanvas->SetBottomMargin(0.16);
				MCERPlotCanvas->SetLeftMargin(0.1);
//				MCERPlotCanvas->SetRightMargin(0.25);
				MCERPlotCanvas->SetRightMargin(0.05);
				MCERPlotCanvas->SetTopMargin(0.15);			

				legMC = new TLegend(0.15,0.89,0.9,0.99);
				legMC->SetBorderSize(0);
				legMC->SetTextSize(0.04);
				legMC->SetTextFont(FontStyle);
				legMC->SetNColumns(5);			

			}

			// ----------------------------------------------------------------------------------------------------

			// Loop over the samples / sources of uncertainty

			int Counter = 0; // Counter for colors

			TH2D* DetFracCovMatrix = nullptr; // matrix for all detector variations

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				// -----------------------------------------------------------------------------------------------------------------------------------------

				// Opening the file containing the covariance matrices for each one of the systematics

//				TString FileCovarianceSpecName = BeamOn9+"WienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
				TString FileCovarianceSpecName = "WienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
				if (BeamOn9 != "" && UncSources[WhichSample] == "Stat") { FileCovarianceSpecName = Tune + "WienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }
				if (BeamOn9 != "" && UncSources[WhichSample] == "NuWro" && ( Tune == "NoTune" || Tune == "TwiceMEC" ) ) 
					{ FileCovarianceSpecName = "WienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }


				// For the detector variation, we follow the PeLEE recipe
				// Only Run3 and propagate across all runs

				if (UncSources[WhichSample] == "LY" || UncSources[WhichSample] == "TPC"  || UncSources[WhichSample] == "SCERecomb2" || UncSources[WhichSample] == "MC_LY" 
				|| UncSources[WhichSample] == "MC_TPC" || UncSources[WhichSample] == "MC_SCERecomb2" || UncSources[WhichSample] == "SmEff_LY" 
				|| UncSources[WhichSample] == "SmEff_TPC" || UncSources[WhichSample] == "SmEff_SCERecomb2" ) {

//					FileCovarianceSpecName = BeamOn9+"WienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_Run3_"+UBCodeVersion+".root";
					FileCovarianceSpecName = Tune + "WienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_Run3_"+UBCodeVersion+".root";

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
				TH2D* LocalCovMatrixClone = (TH2D*)(LocalCovMatrix->Clone() );
				TH2D* LocalFracCovMatrix = (TH2D*)( CovFiles[WhichSample]->Get(LocalFracCovMatrixName) );			

				LocalCovMatrix->SetDirectory(0);
				LocalFracCovMatrix->SetDirectory(0);				
				LocalCovMatrixClone->SetDirectory(0);				
				CovFiles[WhichSample]->Close();

				// -------------------------------------------------------------------------------------------------

				// Stat
				if (WhichSample == 0) { 

					StatCovariances[WhichPlot] = LocalCovMatrix; 
					FracStatCovariances[WhichPlot] = LocalFracCovMatrix;					
				
				} else {

					if (WhichSample == 1) { 

						if ( StatCovariances[WhichPlot]->GetNbinsX() != LocalCovMatrixClone->GetNbinsX() )
							{ cout << PlotNames[WhichPlot] << " covariance matrix with different number of bins" << endl; }
						
						SystCovariances[WhichPlot] = LocalCovMatrixClone; 
						FracSystCovariances[WhichPlot] = LocalFracCovMatrix;

					} else { 

						if ( StatCovariances[WhichPlot]->GetNbinsX() != LocalCovMatrixClone->GetNbinsX() )
							{ cout << PlotNames[WhichPlot] << " covariance matrix with different number of bins" << endl; }						
						
						SystCovariances[WhichPlot]->Add(LocalCovMatrixClone); 
						FracSystCovariances[WhichPlot]->Add(LocalFracCovMatrix);						 
						
					}					

					if (UncSources[WhichSample] == "MC_Stat") { MCStatCovariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "LY") { LYCovariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "TPC") { TPCCovariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "SCERecomb2") { SCERecomb2Covariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "XSec") { XSecCovariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "G4") { G4Covariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "Flux") { FluxCovariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "Dirt") { DirtCovariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "POT") { POTCovariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "NuWro") { NuWroCovariances[WhichPlot] = LocalCovMatrix; }					
					if (UncSources[WhichSample] == "NTarget") { NTargetCovariances[WhichPlot] = LocalCovMatrix; }																				

				}

				// --------------------------------------------------------------------------------

				if (StorePlots) {

					TLegend* leg = nullptr;

					MCERPlotCanvas->cd();  leg = legMC;

/*					if ( !( string(UncSources[WhichSample]).find("POT") != std::string::npos || string(UncSources[WhichSample]).find("NTarget") != std::string::npos ) ) {*/

						//----------------------------------------//
/*
						// Only for the single bin plot

						if (PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot") {

							// Merge the detector variations

							if (UncSources[WhichSample] == "LY" || UncSources[WhichSample] == "TPC" || UncSources[WhichSample] == "SCERecomb2") {

								if (UncSources[WhichSample] == "LY") {

									DetFracCovMatrix = LocalFracCovMatrix;

								} else if (UncSources[WhichSample] == "TPC") {

									//DetFracCovMatrix->Add(LocalFracCovMatrix);

								} else { // SCERecomb2

									//DetFracCovMatrix->Add(LocalFracCovMatrix);
									ReturnUncPlot(DetFracCovMatrix,PlotNames[WhichPlot],Runs[WhichRun],"Det",Counter,leg);
									Counter++;									

								}

							} else { // Explicitly include evenrything else

								ReturnUncPlot(LocalFracCovMatrix,PlotNames[WhichPlot],Runs[WhichRun],UncSources[WhichSample],Counter,leg);
								Counter++;

							}


						} else { // For all the other plots
*/
							ReturnUncPlot(LocalFracCovMatrix,PlotNames[WhichPlot],Runs[WhichRun],UncSources[WhichSample],Counter,leg);
							Counter++;

						/*}*/

						//----------------------------------------//

					/*}*/			

				}

			} // End of the loop over the samples

			// ------------------------------------------------------------------			

			// Storing Stat/Syst/Total Covariance Matrices in root file

			TotalFileCovarianceMatrices->cd();

			MCStatCovariances[WhichPlot]->Write("MCStatCovariance_"+PlotNames[WhichPlot]);

			if (BeamOn9 == "") {

				LYCovariances[WhichPlot]->Write("LYCovariance_"+PlotNames[WhichPlot]);
				TPCCovariances[WhichPlot]->Write("TPCCovariance_"+PlotNames[WhichPlot]);			
				SCERecomb2Covariances[WhichPlot]->Write("SCERecomb2Covariance_"+PlotNames[WhichPlot]);			
				XSecCovariances[WhichPlot]->Write("XSecCovariance_"+PlotNames[WhichPlot]);
				G4Covariances[WhichPlot]->Write("G4Covariance_"+PlotNames[WhichPlot]);			
				FluxCovariances[WhichPlot]->Write("FluxCovariance_"+PlotNames[WhichPlot]);
				DirtCovariances[WhichPlot]->Write("DirtCovariance_"+PlotNames[WhichPlot]);
				POTCovariances[WhichPlot]->Write("POTCovariance_"+PlotNames[WhichPlot]);
				NuWroCovariances[WhichPlot]->Write("NuWroCovariance_"+PlotNames[WhichPlot]);							
				NTargetCovariances[WhichPlot]->Write("NTargetCovariance_"+PlotNames[WhichPlot]);

			}			

			StatCovariances[WhichPlot]->Write("StatCovariance_"+PlotNames[WhichPlot]);
			FracStatCovariances[WhichPlot]->Write("FracStatCovariance_"+PlotNames[WhichPlot]);		

			SystCovariances[WhichPlot]->Write("SystCovariance_"+PlotNames[WhichPlot]);
			FracSystCovariances[WhichPlot]->Write("FracSystCovariance_"+PlotNames[WhichPlot]);		

			TH2D* CloneCovariances = (TH2D*)(StatCovariances[WhichPlot]->Clone());
			if ( SystCovariances[WhichPlot]->GetNbinsX() != CloneCovariances->GetNbinsX() )
				{ cout << PlotNames[WhichPlot] << " covariance matrix with different number of bins" << endl; }			
			CloneCovariances->Add(SystCovariances[WhichPlot]);
			CloneCovariances->Write("TotalCovariance_"+PlotNames[WhichPlot]);

			TH2D* CloneFracCovariances = (TH2D*)(FracStatCovariances[WhichPlot]->Clone());
			if ( FracSystCovariances[WhichPlot]->GetNbinsX() != CloneFracCovariances->GetNbinsX() )
				{ cout << PlotNames[WhichPlot] << " covariance matrix with different number of bins" << endl; }			
			CloneFracCovariances->Add(FracSystCovariances[WhichPlot]);
			CloneFracCovariances->Write("FracTotalCovariance_"+PlotNames[WhichPlot]);		

			// ------------------------------------------------------------------

			if (StorePlots) {

				// ------------------------------------------------------------------

				MCERPlotCanvas->cd();
				legMC->Draw();

				// ------------------------------------------------------------------

				// Plot the total unc on top of everything else

				ReturnUncPlot(CloneFracCovariances,PlotNames[WhichPlot],Runs[WhichRun],"Total",-1,legMC);
		
				// ------------------------------------------------------------------

				MCERPlotCanvas->SaveAs(PlotPath+OverlaySample+"/"+MCERCanvasName+".pdf");
				delete 	MCERPlotCanvas;

				// ---------------------------------------------------------------------------------------------

				// Plot the 2D total covariance matrices

				PlotCov(CloneCovariances,"",PlotNames[WhichPlot],OverlaySample,Runs[WhichRun],Tune);

				// ---------------------------------------------------------------------------------------------

				// Plot the 2D total fractional covariance matrices

				PlotCov(CloneFracCovariances,"Frac",PlotNames[WhichPlot],OverlaySample,Runs[WhichRun],Tune);

				// -------------------------------------------------------------------------------------------	

				// Store correlation matrices

				TH2D* CorrMatrix = (TH2D*)(CloneCovariances->Clone());
				int NBins = CorrMatrix->GetNbinsX();				

				for (int WhichXBin = 0; WhichXBin < NBins; WhichXBin++) { 

					for (int WhichYBin = 0; WhichYBin < NBins; WhichYBin++) {	

						double BinValue = CloneCovariances->GetBinContent(WhichXBin+1,WhichYBin+1);
						double XBinValue = CloneCovariances->GetBinContent(WhichXBin+1,WhichXBin+1);
						double YBinValue = CloneCovariances->GetBinContent(WhichYBin+1,WhichYBin+1);						
						double CorrBinValue = BinValue / ( TMath::Sqrt(XBinValue) * TMath::Sqrt(YBinValue) ); 

						CorrMatrix->SetBinContent(WhichXBin+1,WhichYBin+1,CorrBinValue);

					}

				}			

				PlotCov(CorrMatrix, "Corr", PlotNames[WhichPlot], OverlaySample, Runs[WhichRun],Tune);											

				// ---------------------------------------------------------------------------------------------

			}
		
		} // End of the loop over the plots

		TotalFileCovarianceMatrices->Close();
		cout << endl << "Covariance matrix file " << TotalFileCovarianceName << " has been created" << endl << endl;

		cout << "Merging of covariance matrices for run " << Runs[WhichRun] << " completed!" << endl;

	} // End of the loop over the runs	

} // End of the program 