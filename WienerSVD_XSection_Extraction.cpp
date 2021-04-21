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

#include "ubana/AnalysisCode/Secondary_Code/myFunctions.cpp"

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

void ReweightXSec(TH1D* h, double SF = 1.) {

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

int LocateBinWithValue(TH1D* h, double Value) {

	int NBins = h->GetXaxis()->GetNbins();

	for (int i = 1; i <= NBins; i++) {

		double CurrentEntry = h->GetBinContent(i);
		if (CurrentEntry == Value) { return i; } 

	}

	return -99;

}

// -------------------------------------------------------------------------------------------------------------------------------------

void WienerSVD_XSection_Extraction(TString OverlaySample = "", bool ClosureTest = false) {

	// -------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();	
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
	NameOfSamples.push_back("GenieOverlay");	 
	
	int DataIndex = -1.;

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

	// CV Flux File

	TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));

	// -------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------
	
		double DataPOT = ReturnBeamOnRunPOT(Runs[WhichRun]);						
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface) * (SoftFidSurface / Nominal_UB_XY_Surface);					

		// -------------------------------------------------------------------------------------		

		vector<TCanvas*> PlotCanvas; PlotCanvas.clear();

		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
		vector<vector<TH1D*> > PlotsBkgReco; PlotsBkgReco.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();

		vector<TH2D*> ResponseMatrices; ResponseMatrices.clear();
		vector<TH2D*> CovarianceMatrices; CovarianceMatrices.clear();

		vector<TH2D*> StatCovarianceMatrices; StatCovarianceMatrices.clear();	
		vector<TH2D*> SystCovarianceMatrices; SystCovarianceMatrices.clear();
		vector<TH2D*> POTCovarianceMatrices; POTCovarianceMatrices.clear();
		vector<TH2D*> NTargetCovarianceMatrices; NTargetCovarianceMatrices.clear();
		vector<TH2D*> LYCovarianceMatrices; LYCovarianceMatrices.clear();
		vector<TH2D*> TPCCovarianceMatrices; TPCCovarianceMatrices.clear();
		vector<TH2D*> XSecCovarianceMatrices; XSecCovarianceMatrices.clear();									
		vector<TH2D*> G4CovarianceMatrices; G4CovarianceMatrices.clear();
		vector<TH2D*> FluxCovarianceMatrices; FluxCovarianceMatrices.clear();

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString FileResponseName = MigrationMatrixPath+"FileResponseMatrices_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
		TFile* FileResponseMatrices = new TFile(FileResponseName,"readonly");

		TString FileCovarianceName = MigrationMatrixPath+"WienerSVD_Total_CovarianceMatrices_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";

		TFile* FileCovarianceMatrices = new TFile(FileCovarianceName,"readonly");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();
		
		TString PathToFilesUBCodeExtension = PathToFiles+CutExtension;

		// Store the extracted xsections & associated files in dedicated file

		TFile* ExtractedXSec = nullptr;
		TString NameExtractedXSec = "";

		if (ClosureTest == false) {

			NameExtractedXSec = PathToExtractedXSec+"WienerSVD_ExtractedXSec_"+NameOfSamples[0]+"_"\
				+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+Subtract+".root";
			ExtractedXSec = TFile::Open(NameExtractedXSec,"recreate");

		}		

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

			if (NameOfSamples[WhichSample] == "GenieOverlay") { 
			
				TString FileName = "TruthSTVAnalysis_Overlay9_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
				FileSample.push_back(TFile::Open(PathToFiles+FileName));  
				
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

			StatCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("StatCovariance_"+PlotNames[WhichPlot]));
			SystCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("SystCovariance_"+PlotNames[WhichPlot]));

			// ---------------------------------------------------------------------------------------------------------------------------

			// True CC1p Signal MC // No detector effects

			int n = PlotsTrue[4][WhichPlot]->GetNbinsX();
			double Nuedges[n+1];
			    
			for (int i = 0; i < n+1; i++) { Nuedges[i] = PlotsTrue[4][WhichPlot]->GetBinLowEdge(i+1); }

			// ---------------------------------------------------------------------------------------------------------------------------

			// BeamOn

			TH1D* DataPlot = (TH1D*)(PlotsReco[1][WhichPlot]->Clone());

			DataPlot->Add(PlotsReco[2][WhichPlot],-1); // Subtract ExtBNB
			DataPlot->Add(PlotsReco[3][WhichPlot],-1); // Subtract Dirt
			DataPlot->Add(PlotsBkgReco[0][WhichPlot],-1); // Subtract NonCC1p Beam Related Background

			if (ClosureTest == true) { DataPlot = PlotsCC1pReco[0][WhichPlot]; }

			int m = DataPlot->GetNbinsX();
			TString XTitle = DataPlot->GetXaxis()->GetTitle();
			TString YTitle = DataPlot->GetYaxis()->GetTitle();			 

			// ---------------------------------------------------------------------------------------------------------------------------

			// Construct vectors (for 1D histogram) and matrices (for 2D histogram) for input

			TVectorD signal(n);
			TVectorD measure(m);
			TMatrixD response(m, n);
			TMatrixD covariance(m, m);

			TMatrixD covarianceStat(m, m);
			TMatrixD covarianceSyst(m, m);																											

			// Convert input into mathematical formats, easy and clean to be processed. 
			// Converted defined/implemented in source files, see include/Util.h

			H2V(PlotsTrue[4][WhichPlot], signal);
			H2V(DataPlot, measure);
			H2M(ResponseMatrices[WhichPlot], response, kFALSE); // X axis: Reco, Y axis: True
			H2M(CovarianceMatrices[WhichPlot], covariance, kTRUE);

			H2M(StatCovarianceMatrices[WhichPlot], covarianceStat, kTRUE);
			H2M(SystCovarianceMatrices[WhichPlot], covarianceSyst, kTRUE);

			// ------------------------------------------------------------------------------------------

			// Construct to record additinal smearing matrix and wiener filter (diagomal matrix) elements. 
    
			TMatrixD AddSmear(n,n);

			TMatrixD StatAddSmear(n,n);
			TMatrixD SystAddSmear(n,n);

			// ------------------------------------------------------------------------------------------			

			TVectorD WF(n);

			TVectorD StatWF(n);
			TVectorD SystWF(n);

			// ------------------------------------------------------------------------------------------			

			TMatrixD UnfoldCov(n,n);

			TMatrixD StatUnfoldCov(n,n);
			TMatrixD SystUnfoldCov(n,n);	

			// ------------------------------------------------------------------------------------------	

			TH2D* smear = new TH2D("smear_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"Additional Smearing Matirx",n,Nuedges,n,Nuedges);
			TH1D* wiener = new TH1D("wiener_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"Wiener Filter Vector",n,0,n);
			TH2D* unfcov = new TH2D("unfcov_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"Unfolded spectrum covariance", n, Nuedges, n, Nuedges);

			// --------------------------------------------------------------------------------------------------

			// Use corresponding matrices for each one of the sources of systematic uncertainty

			// Core implementation of Wiener-SVD
			// AddSmear and WF to record the core information in the unfolding.

			TVectorD unfold = WienerSVD(response, signal, measure, covariance, 2, 0, AddSmear, WF, UnfoldCov);

			TVectorD unfoldStat = WienerSVD(response, signal, measure, covarianceStat, 2, 0, StatAddSmear, StatWF, StatUnfoldCov);			
			TVectorD unfoldSyst = WienerSVD(response, signal, measure, covarianceSyst, 2, 0, SystAddSmear, SystWF, SystUnfoldCov);

			// --------------------------------------------------------------------------------------------------

			// Start plotting
		
			TString CanvasName = PlotNames[WhichPlot]+"_"+Runs[WhichRun];
			TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.16);
			PlotCanvas->SetTopMargin(0.13);			
			PlotCanvas->SetLeftMargin(0.2);			
		
			TH1D* unf = new TH1D("unf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);

			TH1D* unfStat = new TH1D("unfStat_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* unfSyst = new TH1D("unfSyst_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);			

			TH1D* unfStatUnc = new TH1D("unfStatUnc_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* unfSystUnc = new TH1D("unfSystUnc_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);

			// --------------------------------------------------------------------------------------------------

			V2H(unfold, unf);

			V2H(unfoldStat, unfStat);
			V2H(unfoldSyst, unfSyst);	

			// --------------------------------------------------------------------------------------------------																												

			// Scaling by correct units & bin width division

			unf->Scale(Units/(IntegratedFlux*NTargets));

			unfStat->Scale(Units/(IntegratedFlux*NTargets));
			unfSyst->Scale(Units/(IntegratedFlux*NTargets));

			for (int i = 1; i <= n;i++ ) { 

				double CVInBin = unf->GetBinContent(i);
				double CovUnc = TMath::Sqrt(UnfoldCov(i-1,i-1) );
				double StatCovUnc = TMath::Sqrt(StatUnfoldCov(i-1,i-1) );
				double SystCovUnc = TMath::Sqrt(SystUnfoldCov(i-1,i-1) );

				unf->SetBinError(i, CovUnc );

				//unfStat->SetBinContent(i, CVInBin );
				//unfSyst->SetBinContent(i, CVInBin );

				unfStat->SetBinError(i, StatCovUnc );
				unfSyst->SetBinError(i, SystCovUnc );

				unfStatUnc->SetBinContent(i, StatCovUnc / CovUnc * 100.);
				unfSystUnc->SetBinContent(i, SystCovUnc / CovUnc * 100.);

			}

			ReweightXSec(unf);

			ReweightXSec(unfStat);
			ReweightXSec(unfSyst);																									

			// ------------------------------------------------------------------------------------

			// Start plotting 					

			unf->GetXaxis()->CenterTitle();
			unf->GetXaxis()->SetLabelSize(TextSize);
			unf->GetXaxis()->SetLabelFont(FontStyle);
			unf->GetXaxis()->SetTitleSize(TextSize);
			unf->GetXaxis()->SetTitleFont(FontStyle);			
			unf->GetXaxis()->SetNdivisions(6);	

			unf->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);
			unf->GetYaxis()->CenterTitle();
			unf->GetYaxis()->SetLabelSize(TextSize);
			unf->GetYaxis()->SetLabelFont(FontStyle);
			unf->GetYaxis()->SetTitleSize(TextSize);
			unf->GetYaxis()->SetTitleFont(FontStyle);			
			unf->GetYaxis()->SetNdivisions(6);

			double MaxValue = unf->GetMaximum();
			int MaxValueBin = LocateBinWithValue(unf,MaxValue);
			double MaxValueError = unf->GetBinError(MaxValueBin);

			double MinValue = unf->GetMinimum();
			int MinValueBin = LocateBinWithValue(unf,MinValue);
			double MinValueError = unf->GetBinError(MinValueBin);			

			double min = TMath::Min(0., 0.8*(MinValue-MinValueError));
			double max = TMath::Max(0., 1.2*(MaxValue+MaxValueError));											
			unf->GetYaxis()->SetRangeUser(min,max);
			unf->SetLineColor(BeamOnColor);
			unf->SetMarkerColor(BeamOnColor);
			unf->SetMarkerStyle(20);
			unf->SetMarkerSize(2.);	

			PlotCanvas->cd();
			if (ClosureTest == true) { unf->Draw("p0 hist"); }
			else { unf->Draw("ex0"); }			

			// The MC CC1p prediction has to be multiplied by the smearing matrix Ac

			TH1D* TrueUnf = new TH1D("TrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TVectorD AcTrueUnfold = AddSmear * signal;
			V2H(AcTrueUnfold, TrueUnf);

			ReweightXSec(TrueUnf);
			TrueUnf->Scale(Units/(IntegratedFlux*NTargets));
			TrueUnf->SetLineColor(OverlayColor);
			//TrueUnf->SetFillColor(OverlayColor);
			TrueUnf->SetMarkerColor(OverlayColor);	
			PlotCanvas->cd();					
			TrueUnf->Draw("hist same");

			PlotCanvas->cd();
			if (ClosureTest == true) { unf->Draw("p0 hist same"); }
			else { 
				
				unf->Draw("ex0 same"); 

				// unfStat->SetLineColor(kRed);
				// unfStat->Draw("ex0 same");

				// unfSyst->SetLineColor(kBlue);
				// unfSyst->Draw("ex0 same");				
				
			}

			// ------------------------------------------------------------------------------

			// Legend & POT Normalization

			double tor860_wcut = ReturnBeamOnRunPOT(Runs[WhichRun]);
			TString Label = ToStringPOT(tor860_wcut)+" POT";

			TLegend* legData = new TLegend(0.15,0.89,0.95,0.98);
			legData->SetBorderSize(0);
			legData->SetTextSize(0.06);
			legData->SetTextFont(FontStyle);
			legData->SetNColumns(1);
			legData->SetMargin(0.1);			

			legData->AddEntry(unf,"MicroBooNE Data " + Runs[WhichRun] + " " + Label,"ep");

			TLegendEntry* lMC = legData->AddEntry(TrueUnf,"MC (uB Tune v2)","l");
			lMC->SetTextColor(OverlayColor);	

			legData->Draw();
			
			TString CanvasPath = PlotPath+NameOfSamples[0];
			TString FullCanvasName = "/WienerSVD_XSections_"+CanvasName+"_"+UBCodeVersion+Subtract+".pdf";
			if (ClosureTest == true) { FullCanvasName = "/ClosureTest_WienerSVD_XSections_"+CanvasName+"_"+UBCodeVersion+Subtract+".pdf"; }
			PlotCanvas->SaveAs(CanvasPath+FullCanvasName);	
			delete PlotCanvas;			

			// ---------------------------------------------------------------------------------------

			// Systematics contributions on a bin-by-bin basis if NOT closure test only	

			if (!ClosureTest && OverlaySample == "") {

				// Unfolding procedure using the covariance matrix for each one of the sources of systematic uncertainty

				POTCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("POTCovariance_"+PlotNames[WhichPlot]));			
				NTargetCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("NTargetCovariance_"+PlotNames[WhichPlot]));
				LYCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("LYCovariance_"+PlotNames[WhichPlot]));			
				TPCCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("TPCCovariance_"+PlotNames[WhichPlot]));
				XSecCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("XSecCovariance_"+PlotNames[WhichPlot]));
				G4CovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("G4Covariance_"+PlotNames[WhichPlot]));
				FluxCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("FluxCovariance_"+PlotNames[WhichPlot]));

				TMatrixD covariancePOT(m, m);
				TMatrixD covarianceNTarget(m, m);
				TMatrixD covarianceLY(m, m);
				TMatrixD covarianceTPC(m, m);
				TMatrixD covarianceXSec(m, m);
				TMatrixD covarianceG4(m, m);
				TMatrixD covarianceFlux(m, m);

				H2M(POTCovarianceMatrices[WhichPlot], covariancePOT, kTRUE);
				H2M(NTargetCovarianceMatrices[WhichPlot], covarianceNTarget, kTRUE);												
				H2M(LYCovarianceMatrices[WhichPlot], covarianceLY, kTRUE);
				H2M(TPCCovarianceMatrices[WhichPlot], covarianceTPC, kTRUE);
				H2M(XSecCovarianceMatrices[WhichPlot], covarianceXSec, kTRUE);
				H2M(G4CovarianceMatrices[WhichPlot], covarianceG4, kTRUE);
				H2M(FluxCovarianceMatrices[WhichPlot], covarianceFlux, kTRUE);

				TMatrixD POTAddSmear(n,n);						
				TMatrixD NTargetAddSmear(n,n);
				TMatrixD LYAddSmear(n,n);
				TMatrixD TPCAddSmear(n,n);
				TMatrixD XSecAddSmear(n,n);
				TMatrixD G4AddSmear(n,n);
				TMatrixD FluxAddSmear(n,n);

				TVectorD POTWF(n);
				TVectorD NTargetWF(n);
				TVectorD LYWF(n);
				TVectorD TPCWF(n);
				TVectorD XSecWF(n);
				TVectorD G4WF(n);
				TVectorD FluxWF(n);					

				TMatrixD POTUnfoldCov(n,n);	
				TMatrixD NTargetUnfoldCov(n,n);
				TMatrixD LYUnfoldCov(n,n);
				TMatrixD TPCUnfoldCov(n,n);
				TMatrixD XSecUnfoldCov(n,n);
				TMatrixD G4UnfoldCov(n,n);
				TMatrixD FluxUnfoldCov(n,n);				

				TVectorD unfoldPOT = WienerSVD(response, signal, measure, covariancePOT, 2, 0, POTAddSmear, POTWF, POTUnfoldCov);						
				TVectorD unfoldNTarget = WienerSVD(response, signal, measure, covarianceNTarget, 2, 0, NTargetAddSmear, NTargetWF, NTargetUnfoldCov);		
				TVectorD unfoldLY = WienerSVD(response, signal, measure, covarianceLY, 2, 0, LYAddSmear, LYWF, LYUnfoldCov);
				TVectorD unfoldTPC = WienerSVD(response, signal, measure, covarianceTPC, 2, 0, TPCAddSmear, TPCWF, TPCUnfoldCov);			
				TVectorD unfoldXSec = WienerSVD(response, signal, measure, covarianceXSec, 2, 0, XSecAddSmear, XSecWF, XSecUnfoldCov);
				TVectorD unfoldG4 = WienerSVD(response, signal, measure, covarianceG4, 2, 0, G4AddSmear, G4WF, G4UnfoldCov);			
				TVectorD unfoldFlux = WienerSVD(response, signal, measure, covarianceFlux, 2, 0, FluxAddSmear, FluxWF, FluxUnfoldCov);
		
				// ---------------------------------------------------------------------------------------------------

				// Declaring histos

				TH1D* unfPOT = new TH1D("unfPOT_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);	
				TH1D* unfNTarget = new TH1D("unfNTarget_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
				TH1D* unfLY = new TH1D("unfLY_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
				TH1D* unfTPC = new TH1D("unfTPC_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
				TH1D* unfXSec = new TH1D("unfXSec_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
				TH1D* unfG4 = new TH1D("unfG4_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
				TH1D* unfFlux = new TH1D("unfFlux_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
		
				TH1D* unfPOTUnc = new TH1D("unfPOTUnc_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);	
				TH1D* unfNTargetUnc = new TH1D("unfNTargetUnc_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
				TH1D* unfLYUnc = new TH1D("unfLYUnc_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
				TH1D* unfTPCUnc = new TH1D("unfTPCUnc_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
				TH1D* unfXSecUnc = new TH1D("unfXSecUnc_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
				TH1D* unfG4Unc = new TH1D("unfG4Unc_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
				TH1D* unfFluxUnc = new TH1D("unfFluxUnc_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);

				TH1D* unfPOTUncOverCV = new TH1D("unfPOTUncOverCV_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);	
				TH1D* unfNTargetUncOverCV = new TH1D("unfNTargetUncOverCV_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
				TH1D* unfLYUncOverCV = new TH1D("unfLYUncOverCV_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
				TH1D* unfTPCUncOverCV = new TH1D("unfTPCUncOverCV_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
				TH1D* unfXSecUncOverCV = new TH1D("unfXSecUncOverCV_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
				TH1D* unfG4UncOverCV = new TH1D("unfG4UncOverCV_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
				TH1D* unfFluxUncOverCV = new TH1D("unfFluxUncOverCV_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);				

				V2H(unfoldPOT, unfPOT);
				V2H(unfoldNTarget, unfNTarget);	
				V2H(unfoldLY, unfLY);
				V2H(unfoldTPC, unfTPC);
				V2H(unfoldXSec, unfXSec);
				V2H(unfoldG4, unfG4);
				V2H(unfoldFlux, unfFlux);		

				// ----------------------------------------------------------------------------------------------------

				// Scaling factors 

				unfPOT->Scale(Units/(IntegratedFlux*NTargets));						
				unfNTarget->Scale(Units/(IntegratedFlux*NTargets));
				unfLY->Scale(Units/(IntegratedFlux*NTargets));
				unfTPC->Scale(Units/(IntegratedFlux*NTargets));
				unfXSec->Scale(Units/(IntegratedFlux*NTargets));
				unfG4->Scale(Units/(IntegratedFlux*NTargets));
				unfFlux->Scale(Units/(IntegratedFlux*NTargets));	

				for (int i = 1; i <= n;i++ ) { 				

					double CVBinEntry = unf->GetBinContent(i);

					double CovUnc = TMath::Sqrt(UnfoldCov(i-1,i-1) );
					double POTCovUnc = TMath::Sqrt(POTUnfoldCov(i-1,i-1) );												
					double NTargetCovUnc = TMath::Sqrt(NTargetUnfoldCov(i-1,i-1) );
					double LYCovUnc = TMath::Sqrt(LYUnfoldCov(i-1,i-1) );
					double TPCCovUnc = TMath::Sqrt(TPCUnfoldCov(i-1,i-1) );
					double XSecCovUnc = TMath::Sqrt(XSecUnfoldCov(i-1,i-1) );				
					double G4CovUnc = TMath::Sqrt(G4UnfoldCov(i-1,i-1) );
					double FluxCovUnc = TMath::Sqrt(FluxUnfoldCov(i-1,i-1) );

					unfPOT->SetBinError(i, POTCovUnc );												 
					unfNTarget->SetBinError(i, NTargetCovUnc );				
					unfLY->SetBinError(i, LYCovUnc );
					unfTPC->SetBinError(i, TPCCovUnc );
					unfXSec->SetBinError(i, XSecCovUnc );
					unfG4->SetBinError(i, G4CovUnc );
					unfFlux->SetBinError(i, FluxCovUnc );
							 
					unfNTargetUnc->SetBinContent(i, NTargetCovUnc / CovUnc * 100.);		
					unfLYUnc->SetBinContent(i, LYCovUnc / CovUnc * 100.);
					unfTPCUnc->SetBinContent(i, TPCCovUnc / CovUnc * 100.);
					unfXSecUnc->SetBinContent(i, XSecCovUnc / CovUnc * 100.);
					unfG4Unc->SetBinContent(i, G4CovUnc / CovUnc * 100.);
					unfFluxUnc->SetBinContent(i, FluxCovUnc / CovUnc * 100.);	

				}				


				ReweightXSec(unfPOT);
				ReweightXSec(unfNTarget);
				ReweightXSec(unfLY);
				ReweightXSec(unfTPC);
				ReweightXSec(unfXSec);
				ReweightXSec(unfG4);
				ReweightXSec(unfFlux);				

				TString CanvasNameUnc = "Unc_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun];
				TCanvas* PlotCanvasUnc = new TCanvas(CanvasNameUnc,CanvasNameUnc,205,34,1024,768);
				PlotCanvasUnc->cd();
				PlotCanvasUnc->SetBottomMargin(0.16);
				PlotCanvasUnc->SetTopMargin(0.13);			
				PlotCanvasUnc->SetLeftMargin(0.2);

				TLegend* legUnc = new TLegend(0.15,0.88,0.95,0.98);
				legUnc->SetBorderSize(0);
				legUnc->SetTextSize(0.06);
				legUnc->SetTextFont(FontStyle);
				legUnc->SetNColumns(5);
				legUnc->SetMargin(0.1);								

				unfStatUnc->GetXaxis()->CenterTitle();
				unfStatUnc->GetXaxis()->SetLabelSize(TextSize);
				unfStatUnc->GetXaxis()->SetLabelFont(FontStyle);
				unfStatUnc->GetXaxis()->SetTitleSize(TextSize);
				unfStatUnc->GetXaxis()->SetTitleFont(FontStyle);			
				unfStatUnc->GetXaxis()->SetNdivisions(6);	

				unfStatUnc->GetYaxis()->SetTitle("Uncertainty Source / Total [%]");
				unfStatUnc->GetYaxis()->CenterTitle();
				unfStatUnc->GetYaxis()->SetLabelSize(TextSize);
				unfStatUnc->GetYaxis()->SetLabelFont(FontStyle);
				unfStatUnc->GetYaxis()->SetTitleSize(TextSize);
				unfStatUnc->GetYaxis()->SetTitleFont(FontStyle);			
				unfStatUnc->GetYaxis()->SetNdivisions(6);
				unfStatUnc->GetYaxis()->SetRangeUser(0.,109.);					

				unfStatUnc->SetLineWidth(3);
				unfPOTUnc->SetLineWidth(3);
				unfNTargetUnc->SetLineWidth(3);
				unfLYUnc->SetLineWidth(3);
				unfTPCUnc->SetLineWidth(3);
				unfXSecUnc->SetLineWidth(3);												
				unfG4Unc->SetLineWidth(3);
				unfFluxUnc->SetLineWidth(3);

				unfStatUnc->SetLineColor(BeamOnColor);
				unfPOTUnc->SetLineColor(OverlayColor);
				unfNTargetUnc->SetLineColor(GiBUUColor);
				unfLYUnc->SetLineColor(NuWroColor);
				unfTPCUnc->SetLineColor(NEUTColor);
				unfXSecUnc->SetLineColor(SuSav2Color);
				unfG4Unc->SetLineColor(GENIEv2Color);
				unfFluxUnc->SetLineColor(GenieColor);																												

				unfStatUnc->Draw("hist text0 same");
				unfPOTUnc->Draw("hist text0 same");				
				unfNTargetUnc->Draw("hist text0 same");
				unfLYUnc->Draw("hist text0 same");				
				unfTPCUnc->Draw("hist text0 same");
				unfXSecUnc->Draw("hist text0 same");
				unfG4Unc->Draw("hist text0 same");
				unfFluxUnc->Draw("hist text0 same");											

				legUnc->AddEntry(unfStatUnc,"Stat","l");
				legUnc->AddEntry(unfPOTUnc,"POT","l");
				legUnc->AddEntry(unfNTargetUnc,"N_{targets}","l");
				legUnc->AddEntry(unfLYUnc,"LY","l");
				legUnc->AddEntry(unfTPCUnc,"TPC","l");
				legUnc->AddEntry(unfXSecUnc,"XSec","l");
				legUnc->AddEntry(unfG4Unc,"G4","l");
				legUnc->AddEntry(unfFluxUnc,"Flux","l");																				

				legUnc->Draw();											

				TString FullCanvasNameUnc = "/Syst_WienerSVD_XSections_"+CanvasName+"_"+UBCodeVersion+Subtract+".pdf";
				PlotCanvasUnc->SaveAs(CanvasPath+FullCanvasNameUnc);
				delete PlotCanvasUnc;

			}							

			// ----------------------------------------------------------------------------------------------------------------		

			TH1D* diff = new TH1D("diff_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"Fractional difference of unf and signal model",n, Nuedges);
			
			for(int i=1; i <= n; i++) {
			
				double s2s = 0;
				double u = unf->GetBinContent(i);
				double t = signal(i-1);
				if (t != 0) { s2s = u/t - 1; }
				else { s2s = 1.; } // attention to this t = 0
				diff->SetBinContent(i, s2s); // in percentage 
			}	

			// ---------------------------------------------------------------------------------------------------------------------------

            // intrinsic bias (Ac-I) * s_bar formula

            TH1D* bias = new TH1D("bias_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"intrinsic bias w.r.t. model",n, Nuedges);
            TH1D* bias2 = new TH1D("bias2_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"intrinsic bias2 w.r.t. unfolded result",n, Nuedges);
            TMatrixD unit(n,n);
            unit.UnitMatrix();
            TVectorD intrinsicbias = (AddSmear - unit)*signal;
            TVectorD intrinsicbias2 = (AddSmear - unit)*unfold;

            for (int i = 0; i < n; i++) {

                if (signal(i) != 0) { intrinsicbias(i) = intrinsicbias(i)/signal(i); }
                else { intrinsicbias(i) = 0.; }
                if (unfold(i) != 0) { intrinsicbias2(i) = intrinsicbias2(i)/unfold(i); }
                else { intrinsicbias2(i) = 0.; }
            }	

		    V2H(intrinsicbias, bias);
            V2H(intrinsicbias2, bias2);					

			// ---------------------------------------------------------------------------------------------------------------------------

            // Diagonal Uncertainty

            TH1D* fracError = new TH1D("fracError_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun], "Fractional uncertainty", n, Nuedges);
            TH1D* absError = new TH1D("absError_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun], "absolute uncertainty", n, Nuedges);

            for (int i = 1; i <= n; i++) {

                fracError->SetBinContent(i, TMath::Sqrt(UnfoldCov(i-1, i-1))/unfold(i-1));
                absError->SetBinContent(i, TMath::Sqrt(UnfoldCov(i-1, i-1)));
            }
    
            /// MSE
            TH1D* MSE = new TH1D("MSE_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun], "Mean Square Error: variance+bias^2", n, Nuedges);
            TH1D* MSE2 = new TH1D("MSE2_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun], "Mean Square Error: variance", n, Nuedges);
    
	        for (int i = 0; i < n; i++) {

                MSE->SetBinContent(i+1, TMath::Power(intrinsicbias2(i)*unfold(i),2)+UnfoldCov(i,i));
                MSE2->SetBinContent(i+1, UnfoldCov(i,i));

            }

            // convert matrix/vector to histogram and save

            M2H(AddSmear, smear);
            V2H(WF, wiener);
            M2H(UnfoldCov, unfcov);

			// ---------------------------------------------------------------------------------------------------------------------------
    
			if (ClosureTest == false ) {

				ExtractedXSec->cd();

				unf->Write("Reco"+PlotNames[WhichPlot]);
				TrueUnf->Write("True"+PlotNames[WhichPlot]);			
				smear->Write("Ac"+PlotNames[WhichPlot]);
				unfcov->Write("UnfCov"+PlotNames[WhichPlot]);	
				CovarianceMatrices[WhichPlot]->Write("Cov"+PlotNames[WhichPlot]);					
				wiener->Write("Wiener"+PlotNames[WhichPlot]);
				diff->Write("Diff"+PlotNames[WhichPlot]);
				bias->Write("Bias"+PlotNames[WhichPlot]);
				bias2->Write("Bias2"+PlotNames[WhichPlot]);
				fracError->Write("FracErr"+PlotNames[WhichPlot]);
				absError->Write("AbsErr"+PlotNames[WhichPlot]);
				MSE->Write("MSE"+PlotNames[WhichPlot]);
				MSE2->Write("MSE2"+PlotNames[WhichPlot]);
				ResponseMatrices[WhichPlot]->Write("Response"+PlotNames[WhichPlot]);	

			}

			// ---------------------------------------------------------------------------------------------------------------------------

		} // End of the loop over the plots

		if (ClosureTest == false ) {

			ExtractedXSec->Close();
			std::cout << std::endl << "File " << NameExtractedXSec << " created" << std::endl << std::endl;

		}

//		for (int i = 0; i < NSamples; i++) { FileSample[i]->Close(); }

	} // End of the loop over the runs	

} // End of the program 
