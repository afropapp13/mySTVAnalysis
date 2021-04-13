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

TString ToStringPOT(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

// -------------------------------------------------------------------------------------------------------------------------------------

void Reweight(TH1D* h, double SF = 1.) {

	int NBins = h->GetXaxis()->GetNbins();

	for (int i = 0; i < NBins; i++) {

		double CurrentEntry = h->GetBinContent(i+1);
		double NewEntry = CurrentEntry * SF / h->GetBinWidth(i+1);

		double CurrentError = h->GetBinError(i+1);
		double NewError = CurrentError * SF / h->GetBinWidth(i+1);

		h->SetBinContent(i+1,NewEntry); 
		h->SetBinError(i+1,NewError); 
//		h->SetBinError(i+1,0.000001); 

	}

}

// -------------------------------------------------------------------------------------------------------------------------------------

void WienerSVD_XSection_Extraction(TString OverlaySample) {

	// -------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;
	gStyle->SetOptStat(0);

	TString Subtract = "";
//	TString Subtract = "_BUnsubtracted";

	int DecimalAccuracy = 2;

	// -------------------------------------------------------------------------------------

	int NEventsPassingSelectionCuts = 0;
	TString CutExtension = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	// v52
	VectorCuts.push_back("");
	VectorCuts.push_back("_PID");
	VectorCuts.push_back("_NuScore");

	int NCuts = (int)(VectorCuts.size());	

	for (int i = 0; i < NCuts; i++) {

		CutExtension = CutExtension + VectorCuts[i];

	}				

	// -------------------------------------------------------------------------------------

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

	const int N1DPlots = PlotNames.size();
	//cout << "Number of 1D Plots = " << N1DPlots << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");

	vector<TString> LabelsOfSamples;
	vector<TString> NameOfSamples;

	NameOfSamples.push_back("Overlay9"); 
	NameOfSamples.push_back("BeamOn9"); 
	NameOfSamples.push_back("ExtBNB9"); 
	NameOfSamples.push_back("OverlayDirt9"); 
	
	int DataIndex = -1.;
	double DataPOT = -99.;

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

	// CV Flux File

	TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));

	// -------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------
	
		if (Runs[WhichRun] == "Run1") { DataPOT = tor860_wcut_Run1 ; }
		if (Runs[WhichRun] == "Run2") { DataPOT = tor860_wcut_Run2 ; }
		if (Runs[WhichRun] == "Run3") { DataPOT = tor860_wcut_Run3 ; }
		if (Runs[WhichRun] == "Run4") { DataPOT = tor860_wcut_Run4 ; }
		if (Runs[WhichRun] == "Run5") { DataPOT = tor860_wcut_Run5 ; }	

		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface) * (SoftFidSurface / Nominal_UB_XY_Surface);					

		// -------------------------------------------------------------------------------------		

		vector<TCanvas*> PlotCanvas; PlotCanvas.clear();

		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
		vector<vector<TH1D*> > PlotsBkgReco; PlotsBkgReco.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();

		vector<TH2D*> ResponseMatrices; ResponseMatrices.clear();
		vector<TH2D*> CovarianceMatrices; CovarianceMatrices.clear();

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString FileResponseName = MigrationMatrixPath+"FileResponseMatrices_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
		TFile* FileResponseMatrices = new TFile(FileResponseName,"readonly");

		TString FileCovarianceName = MigrationMatrixPath+"WienerSVD_Total_CovarianceMatrices_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";

		TFile* FileCovarianceMatrices = new TFile(FileCovarianceName,"readonly");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();
		
		TString PathToFilesUBCodeExtension = PathToFiles+CutExtension;

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Loop over the samples

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {
		
			if (NameOfSamples[WhichSample] == "BeamOn9") { DataIndex = WhichSample; }

			if (
				NameOfSamples[WhichSample] == "BeamOn9" || 
				NameOfSamples[WhichSample] == "ExtBNB9" || 
				NameOfSamples[WhichSample] == "OverlayDirt9"
			) { 
			
				TString FileName = "STVStudies_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+CutExtension+".root";
				FileSample.push_back(TFile::Open(PathToFilesUBCodeExtension+"/"+FileName)); 
			}
			
			if (NameOfSamples[WhichSample] == "Overlay9") { 
			
				TString FileName = "STVStudies_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+CutExtension+".root";
				FileSample.push_back(TFile::Open(PathToFilesUBCodeExtension+"/"+FileName)); 
				
			}
						

			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();
			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();
			vector<TH1D*> CurrentPlotsBkgReco; CurrentPlotsBkgReco.clear();
			vector<TH1D*> CurrentPlotsCC1pReco; CurrentPlotsCC1pReco.clear();

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
				CurrentPlotsReco.push_back(histReco);

				TH1D* histBkgReco = (TH1D*)(FileSample[WhichSample]->Get("NonCC1pReco"+PlotNames[WhichPlot]));
				CurrentPlotsBkgReco.push_back(histBkgReco);

				TH1D* histCC1pReco = (TH1D*)(FileSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]));
				CurrentPlotsCC1pReco.push_back(histCC1pReco);

				TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
				CurrentPlotsTrue.push_back(histTrue);
		
			} // End of the loop over the plots

			PlotsReco.push_back(CurrentPlotsReco);		
			PlotsTrue.push_back(CurrentPlotsTrue);		
			PlotsBkgReco.push_back(CurrentPlotsBkgReco);
			PlotsCC1pReco.push_back(CurrentPlotsCC1pReco);

		} // End of the loop over the samples

		// ----------------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

			// ---------------------------------------------------------------------------------------------------------------------------

			ResponseMatrices.push_back((TH2D*)FileResponseMatrices->Get("POTScaledCC1pReco"+PlotNames[WhichPlot]+"2D"));
			CovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("TotalCovariance_"+PlotNames[WhichPlot]));

			// ---------------------------------------------------------------------------------------------------------------------------

			// CC1p Signal MC

			int n = PlotsCC1pReco[0][WhichPlot]->GetNbinsX();
			double Nuedges[n+1];
			    
			for (int i = 0; i < n+1; i++) { Nuedges[i] = PlotsCC1pReco[0][WhichPlot]->GetBinLowEdge(i+1); }

			// ---------------------------------------------------------------------------------------------------------------------------

			// BeamOn

			PlotsReco[1][WhichPlot]->Add(PlotsReco[2][WhichPlot],-1); // Subtract ExtBNB
			PlotsReco[1][WhichPlot]->Add(PlotsReco[3][WhichPlot],-1); // Subtract Dirt
			PlotsReco[1][WhichPlot]->Add(PlotsBkgReco[0][WhichPlot],-1); // Subtract NonCC1p Beam Related Background

			int m = PlotsReco[1][WhichPlot]->GetNbinsX();

			// ---------------------------------------------------------------------------------------------------------------------------

			// Construct vectors (for 1D histogram) and matrices (for 2D histogram) for input

			TVectorD signal(n);
			TVectorD measure(m);
			TMatrixD response(m, n);
			TMatrixD covariance(m, m);

			// Convert input into mathematical formats, easy and clean to be processed. 
			// Converted defined/implemented in source files, see include/Util.h

			H2V(PlotsCC1pReco[0][WhichPlot], signal);
			H2V(PlotsReco[1][WhichPlot], measure);
			H2M(ResponseMatrices[WhichPlot], response, kTRUE);
			H2M(CovarianceMatrices[WhichPlot], covariance, kTRUE);

			// Construct to record additinal smearing matrix and wiener filter (diagomal matrix) elements. 
    
			TMatrixD AddSmear(n,n);    
			TVectorD WF(n);
			TMatrixD UnfoldCov(n,n);
			TH2D* smear = new TH2D("smear","Additional Smearing Matirx",n,Nuedges,n,Nuedges);
			TH1D* wiener = new TH1D("wiener","Wiener Filter Vector",n,0,n);
			TH2D* unfcov = new TH2D("unfcov","Unfolded spectrum covariance", n, Nuedges, n, Nuedges);

			// Core implementation of Wiener-SVD
			// Input as names read. AddSmear and WF to record the core information in the unfolding.
			TVectorD unfold = WienerSVD(response, signal, measure, covariance, 2, 0, AddSmear, WF, UnfoldCov);
			TH1D* unf = new TH1D("unf","unfolded spectrum",n,Nuedges);
			V2H(unfold, unf);
			Reweight(unf);
			unf->Scale(Units/(IntegratedFlux*NTargets));
unf->Draw("e");

			// ---------------------------------------------------------------------------------------------------------------------------
		
/*
			int NBinsX = PlotsCC1pReco[0][WhichPlot]->GetNbinsX();

			// ---------------------------------------------------------------------------------------------------------------------------

			// Apply the relevant weights

			for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

				double BinWidth = PlotsCC1pReco[0][WhichPlot]->GetBinWidth(WhichXBin+1);

				// -----------------------------------------------------------------------------------------------------------------

				double ScalingFactor = Units / (IntegratedFlux * NTargets * BinWidth);

				// -----------------------------------------------------------------------------------------------------------------
		
				// Effective Efficiencies
				
				double EffectiveEfficiencyXBin = PlotsTEfficiency[0][WhichPlot]->GetBinContent(WhichXBin+1);
				double EffectiveEfficiencyXBinError = PlotsTEfficiency[0][WhichPlot]->GetBinError(WhichXBin+1);

				if (EffectiveEfficiencyXBin == 0) { cout << "Bin with 0 efficiency" << endl; }

				// -----------------------------------------------------------------------------------------------------------------

				// ExtBNB

				double CurrentExtBNBEntry = PlotsReco[2][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentExtBNBError = PlotsReco[2][WhichPlot]->GetBinError(WhichXBin+1);
				double ExtBNBScaledEntry = 0., ExtBNBScaledError = 0.;
				
				if (EffectiveEfficiencyXBin != 0 ) {
				  
					ExtBNBScaledEntry = CurrentExtBNBEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					ExtBNBScaledError = sqrt( 
								TMath::Power(CurrentExtBNBError / EffectiveEfficiencyXBin,2.)  +
								TMath::Power(CurrentExtBNBEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2.)
							     ) * ScalingFactor;
							     
				}

				PlotsReco[2][WhichPlot]->SetBinContent(WhichXBin+1,ExtBNBScaledEntry);
				PlotsReco[2][WhichPlot]->SetBinError(WhichXBin+1,ExtBNBScaledError);

				// -------------------------------------------------------------------------------------------------------

				// Dirt

				double CurrentDirtEntry = PlotsReco[3][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentDirtError = PlotsReco[3][WhichPlot]->GetBinError(WhichXBin+1);
				double DirtScaledEntry = 0., DirtScaledError = 0.;
				
				if (EffectiveEfficiencyXBin != 0 ) {  
				
					DirtScaledEntry = CurrentDirtEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					DirtScaledError = sqrt( 
								TMath::Power(CurrentDirtError / EffectiveEfficiencyXBin,2.)  +
								TMath::Power(CurrentDirtEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2.)
							   ) * ScalingFactor;
							   
				}

				PlotsReco[3][WhichPlot]->SetBinContent(WhichXBin+1,DirtScaledEntry);
				PlotsReco[3][WhichPlot]->SetBinError(WhichXBin+1,DirtScaledError);

				// ----------------------------------------------------------------------------------------------------------------

				// Beam Related Overlay Bkg

				double CurrentBkgEntry = PlotsBkgReco[0][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentBkgError = PlotsBkgReco[0][WhichPlot]->GetBinError(WhichXBin+1);
				double BkgScaledEntry = 0., BkgScaledError = 0.;

				if (EffectiveEfficiencyXBin != 0 ) {
				
					BkgScaledEntry = CurrentBkgEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					BkgScaledError = sqrt( 
								TMath::Power(CurrentBkgError / EffectiveEfficiencyXBin,2.)  +
								TMath::Power(CurrentBkgEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2.)
							  ) * ScalingFactor;

				}

				PlotsBkgReco[0][WhichPlot]->SetBinContent(WhichXBin+1,BkgScaledEntry);
				PlotsBkgReco[0][WhichPlot]->SetBinError(WhichXBin+1,BkgScaledError);

				// ------------------------------------------------------------------------------------------------------------------

				// Overlay

				double CurrentOverlayEntry = PlotsCC1pReco[0][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentOverlayError = PlotsCC1pReco[0][WhichPlot]->GetBinError(WhichXBin+1);
				double OverlayScaledEntry = 0., OverlayScaledError = 0.;

				if (EffectiveEfficiencyXBin != 0 ) {
				
					OverlayScaledEntry = CurrentOverlayEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					OverlayScaledError = sqrt( 
								TMath::Power(CurrentOverlayError / EffectiveEfficiencyXBin,2.)  +
								TMath::Power(CurrentOverlayEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2.)
							      ) * ScalingFactor;

				}

				if (Subtract == "_BUnsubtracted") {

					OverlayScaledEntry += BkgScaledEntry;
					OverlayScaledError = TMath::Sqrt( TMath::Power(OverlayScaledError,2.) + TMath::Power(BkgScaledError,2.) );

				}

				PlotsCC1pReco[0][WhichPlot]->SetBinContent(WhichXBin+1,OverlayScaledEntry);
				PlotsCC1pReco[0][WhichPlot]->SetBinError(WhichXBin+1,OverlayScaledError);
				
				// --------------------------------------------------------------------------------------------------------------
				
				// Genie Overlay for internal overlay closure test

				double CurrentGenieOverlayEntry = PlotsTrue[4][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentGenieOverlayError = PlotsTrue[4][WhichPlot]->GetBinError(WhichXBin+1);

				double GenieOverlayScaledEntry = 0., GenieOverlayScaledError = 0.;
				GenieOverlayScaledEntry = CurrentGenieOverlayEntry * ScalingFactor; 
				GenieOverlayScaledError = CurrentGenieOverlayError * ScalingFactor; 

				PlotsTrue[4][WhichPlot]->SetBinContent(WhichXBin+1,GenieOverlayScaledEntry);
				PlotsTrue[4][WhichPlot]->SetBinError(WhichXBin+1,GenieOverlayScaledError);
				
				// --------------------------------------------------------------------------------------------------------------------

				// BeamOn

				double CurrentDataEntry = PlotsReco[1][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentDataError = PlotsReco[1][WhichPlot]->GetBinError(WhichXBin+1);
				double DataScaledEntry = 0., DataScaledError = 0.;

				if (EffectiveEfficiencyXBin != 0 ) {  
				
					DataScaledEntry = CurrentDataEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					DataScaledError = sqrt( 
								TMath::Power(CurrentDataError / EffectiveEfficiencyXBin,2.)  +
								TMath::Power(CurrentDataEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2.)
							   ) * ScalingFactor;
							   
				}

				double TotalDataScaledEntry = DataScaledEntry - ExtBNBScaledEntry - BkgScaledEntry - DirtScaledEntry; 
				if (Subtract == "_BUnsubtracted") { TotalDataScaledEntry = DataScaledEntry - ExtBNBScaledEntry - DirtScaledEntry; }

				double TotalDataScaledError = TMath::Sqrt( 
										TMath::Power(DataScaledError,2.) +
										TMath::Power(ExtBNBScaledError,2.) +
										TMath::Power(BkgScaledError,2.) +
										TMath::Power(DirtScaledError,2.)						
									   );	
									   
				if (Subtract == "_BUnsubtracted") { TotalDataScaledError = TMath::Sqrt( 
													TMath::Power(DataScaledError,2.) +
													TMath::Power(ExtBNBScaledError,2.) +
													TMath::Power(DirtScaledError,2.)						
									   				); 
								  }
									   
				PlotsReco[1][WhichPlot]->SetBinContent(WhichXBin+1,TotalDataScaledEntry);
				PlotsReco[1][WhichPlot]->SetBinError(WhichXBin+1,TotalDataScaledError);

			} // End of the loop over the bins

			// --------------------------------------------------------------------------------------------------------------------

			// Plotting the xsections

			// Overlay

			PlotsCC1pReco[0][WhichPlot]->SetLineColor(OverlayColor);
			PlotsCC1pReco[0][WhichPlot]->SetFillColor(OverlayColor);

			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->CenterTitle();
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetTitleOffset(1.05);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetTitleSize(0.06);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetLabelSize(0.06);		
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetNdivisions(5);			

			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->CenterTitle();
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitleOffset(1.18);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitle(PlotXAxis[WhichPlot]);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitleSize(0.06);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetLabelSize(0.06);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetNdivisions(5);

			// BeamOn

			PlotsReco[1][WhichPlot]->SetLineWidth(3);
			PlotsReco[1][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsReco[1][WhichPlot]->SetMarkerColor(BeamOnColor);
			PlotsReco[1][WhichPlot]->SetMarkerStyle(20);
			PlotsReco[1][WhichPlot]->SetMarkerSize(1.5);

			PlotsTrue[4][WhichPlot]->SetLineWidth(3);	
			PlotsTrue[4][WhichPlot]->SetLineColor(GenieColor);
			PlotsTrue[4][WhichPlot]->SetFillColor(GenieColor);
			PlotsTrue[4][WhichPlot]->SetMarkerColor(GenieColor);
					
			// -----------------------------------------------------------------------------------------------------------------------

			if (OverlaySample == "") {

				// -----------------------------------------------------------------------------------------------------------------------	

				TString CanvasName = PlotNames[WhichPlot]+"_"+Runs[WhichRun]+OverlaySample;
				TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
				PlotCanvas->cd();

				TPad *midPad = new TPad("midPad", "", 0.005,0.0, 0.995, 0.995);
				midPad->SetBottomMargin(0.14);
				midPad->SetTopMargin(0.17);
				midPad->SetLeftMargin(0.16);
				midPad->Draw();

				TLegend* leg = new TLegend(0.15,0.89,0.85,0.99);
				leg->SetBorderSize(0);
				leg->SetTextSize(0.05);
				leg->SetTextFont(FontStyle);
				leg->SetNColumns(4);
				leg->SetColumnSeparation(0.33);			
				
				TLegend* legData = new TLegend(0.15,0.84,0.9,0.88);
				legData->SetBorderSize(0);
				legData->SetTextSize(0.06);
				legData->SetTextFont(FontStyle);
				legData->SetNColumns(3);	

				// Max on plot so that we can include the integrated xsecs
				
				double max = TMath::Max(PlotsCC1pReco[0][WhichPlot]->GetMaximum(),PlotsReco[1][WhichPlot]->GetMaximum());
				double YmaxScaleFactor = 1.35;
				double min = TMath::Min(0.,1.5*PlotsReco[1][WhichPlot]->GetMinimum());
				PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetRangeUser(min,YmaxScaleFactor * max);

				midPad->cd();

				// Plot MC
				PlotsCC1pReco[0][WhichPlot]->Draw("e2");
				
				// Plot GENIE Overlay // Closure test
				if (Subtract == "") { PlotsTrue[4][WhichPlot]->Draw("e same"); }												

				// Plot BeamOn
				PlotsReco[1][WhichPlot]->Draw("ex0 same");
			
				// ----------------------------------------------------------------------------------------------------------------

				// Legend & POT Normalization

				double tor860_wcut = -99.;
				
				if (Runs[WhichRun] == "Run1") { tor860_wcut = tor860_wcut_Run1; }
				if (Runs[WhichRun] == "Run2") { tor860_wcut = tor860_wcut_Run2; }
				if (Runs[WhichRun] == "Run3") { tor860_wcut = tor860_wcut_Run3; }
				if (Runs[WhichRun] == "Run4") { tor860_wcut = tor860_wcut_Run4; }
				if (Runs[WhichRun] == "Run5") { tor860_wcut = tor860_wcut_Run5; }

				TString Label = ToStringPOT(tor860_wcut)+" POT";

				TLegendEntry* lMC = leg->AddEntry(PlotsCC1pReco[0][WhichPlot],"MC","f");
				lMC->SetTextColor(OverlayColor);

				TLegendEntry* lGenie = leg->AddEntry(PlotsTrue[4][WhichPlot],"GENIE Overlay (uB Tune v2)","l");			
				lGenie->SetTextColor(GenieColor);

				leg->Draw();	

				legData->AddEntry(PlotsReco[1][WhichPlot],"MicroBooNE Data " + Runs[WhichRun] + " " + Label,"ep");
				legData->Draw();
			
				TString CanvasPath = PlotPath+NameOfSamples[0];
				TString FullCanvasName = "/XSections_"+CanvasName+"_"+UBCodeVersion+Subtract+".pdf";
				PlotCanvas->SaveAs(CanvasPath+FullCanvasName);

				delete PlotCanvas;

			} // End of the CV case where we plot & store the canvases
			
			// --------------------------------------------------------------------------------------------------------------------------------
*/
		} // End of the loop over the plots
					
		// --------------------------------------------------------------------------------------------------------------------------------
		
//		FileResponseMatrices->Close();
//		FileCovarianceMatrices->Close();

		// --------------------------------------------------------------------------------------------------------------------------------

		// Store the extracted xsections

		TString NameExtractedXSec = PathToExtractedXSec+"WienerSVD_ExtractedXSec_"+NameOfSamples[0]+"_"\
			+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+Subtract+".root";
		TFile* ExtractedXSec = TFile::Open(NameExtractedXSec,"recreate");

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

			PlotsReco[1][WhichPlot]->Write(); // Data Reco
/*			PlotsCC1pReco[0][WhichPlot]->Write(); // Overlay MC		
			PlotsTrue[4][WhichPlot]->Write(); // Genie Overlay // Closure test			
*/
		}

//		ExtractedXSec->Close();


		std::cout << std::endl << "File " << NameExtractedXSec << " created" << std::endl << std::endl;

//		for (int i = 0; i < NSamples; i++) { FileSample[i]->Close(); }

	} // End of the loop over the runs	

} // End of the program 
