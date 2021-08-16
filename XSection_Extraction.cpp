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

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "ubana/AnalysisCode/Secondary_Code/myFunctions.cpp"

TString ToStringPOT(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

void XSection_Extraction(TString OverlaySample,int Universe = -1, bool DetVar = false) { // Universe != -1 ONLY for the flux systematics to deal with the flux of the different universes

	// -------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
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

//	vector<TString> PlotNames;
//	PlotNames.push_back("DeltaPTPlot"); 
//	PlotNames.push_back("DeltaAlphaTPlot"); 
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
//	PlotNames.push_back("kMissPlot");
//	PlotNames.push_back("PMissPlot");
//	PlotNames.push_back("PMissMinusPlot");

//	PlotNames.push_back("CCQEMuonMomentumPlot"); 
//	PlotNames.push_back("CCQEMuonCosThetaPlot"); 
//	PlotNames.push_back("CCQEProtonMomentumPlot"); 
//	PlotNames.push_back("CCQEProtonCosThetaPlot");

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
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");				
	Runs.push_back("Combined");

	if (DetVar) {

		Runs.clear();
		Runs.push_back("Run3");

	}				

	int NRuns = (int)(Runs.size());
	//cout << "Number of Runs = " << NRuns << endl;

	// -------------------------------------------------------------------------------------------------------------------------------------
	
	// CV Flux File

	TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));
	
	// If flux universe, glad the relevant flux universe histo

	if ( Universe != -1 && string(OverlaySample).find("fluxes") != std::string::npos ) {
	
//		TString DublicateOverlaySample = OverlaySample;
//		TString ReducedOverlaySample = DublicateOverlaySample.ReplaceAll("m_","m");
//		if ( !(string(OverlaySample).find("expskin") != std::string::npos) ) { ReducedOverlaySample = ReducedOverlaySample.ReplaceAll("n_","n"); }
//		ReducedOverlaySample = ReducedOverlaySample.ReplaceAll("g_","g");				
//			
//		for (int i = 0; i < 10;i++) { ReducedOverlaySample.ReplaceAll(TString(std::to_string(i)),""); }
//			
//		TString FluxHistoName = "numu_ms"+ReducedOverlaySample+"/hEnumu"+ReducedOverlaySample+"_ms_"+TString(std::to_string(Universe));
//		HistoFlux = (TH1D*)(FluxFile->Get(FluxHistoName));

		TString FluxHistoName = "numu_ms_total/hEnumu_ms_"+TString(std::to_string(Universe));
		HistoFlux = (TH1D*)(FluxFile->Get(FluxHistoName));

	}

	// -------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------
	
		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);													
		
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface) * (SoftFidSurface / Nominal_UB_XY_Surface);	
		//cout << Runs[WhichRun] << " integrated flux = " << IntegratedFlux << endl;		

		// -------------------------------------------------------------------------------------		

		vector<TCanvas*> PlotCanvas; PlotCanvas.clear();

		vector<TH1D*> PlotsRecoUnfOnly; PlotsRecoUnfOnly.resize(N1DPlots);
		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
		vector<vector<TH1D*> > PlotsBkgReco; PlotsBkgReco.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();
		vector<vector<TH1D*> > PlotsTEfficiency; PlotsTEfficiency.clear();

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString FileEfficienciesName = "FileEfficiences_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
		TFile* FileEfficiences = new TFile(FileEfficienciesPath+FileEfficienciesName,"readonly");

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();
		
		TString PathToFilesUBCodeExtension = PathToFiles+CutExtension;

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
			vector<TH1D*> CurrentPlotsTEfficiency; CurrentPlotsTEfficiency.clear();

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

				TH1D* histTEfficiency = (TH1D*)(FileEfficiences->Get("CC1pReco"+PlotNames[WhichPlot]));
				CurrentPlotsTEfficiency.push_back(histTEfficiency);
		
			} // End of the loop over the plots

			PlotsReco.push_back(CurrentPlotsReco);		
			PlotsTrue.push_back(CurrentPlotsTrue);		
			PlotsBkgReco.push_back(CurrentPlotsBkgReco);
			PlotsCC1pReco.push_back(CurrentPlotsCC1pReco);
			PlotsTEfficiency.push_back(CurrentPlotsTEfficiency);


		} // End of the loop over the samples

		// ----------------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {			

			int NBinsX = PlotsCC1pReco[0][WhichPlot]->GetNbinsX();
				
			PlotsRecoUnfOnly[WhichPlot] = (TH1D*)(PlotsReco[1][WhichPlot]->Clone());

			// ---------------------------------------------------------------------------------------------------------------------------

			// Apply the relevant weights

			for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

				double BinWidth = PlotsCC1pReco[0][WhichPlot]->GetBinWidth(WhichXBin+1);

				// We want the number of events, as if the bin width is 1
				if (PlotNames[WhichPlot] == "MuonCosThetaPlot") { BinWidth = 1.; }

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
				double ExtBNBScaledEntry = 0., ExtBNBScaledError = 0., ExtBNBScaledErrorUnfOnly = 0.;
				
				if (EffectiveEfficiencyXBin != 0 ) {
				  
					ExtBNBScaledEntry = CurrentExtBNBEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					ExtBNBScaledError = sqrt( 
								TMath::Power(CurrentExtBNBError / EffectiveEfficiencyXBin,2.)  +
								TMath::Power(CurrentExtBNBEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2.)
							     ) * ScalingFactor;
					ExtBNBScaledErrorUnfOnly = CurrentExtBNBEntry * EffectiveEfficiencyXBinError / (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin) * ScalingFactor;
							     
				}

				PlotsReco[2][WhichPlot]->SetBinContent(WhichXBin+1,ExtBNBScaledEntry);
				PlotsReco[2][WhichPlot]->SetBinError(WhichXBin+1,ExtBNBScaledError);

				// -------------------------------------------------------------------------------------------------------

				// Dirt

				double CurrentDirtEntry = PlotsReco[3][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentDirtError = PlotsReco[3][WhichPlot]->GetBinError(WhichXBin+1);
				double DirtScaledEntry = 0., DirtScaledError = 0., DirtScaledErrorUnfOnly = 0.;
				
				if (EffectiveEfficiencyXBin != 0 ) {  
				
					DirtScaledEntry = CurrentDirtEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					DirtScaledError = sqrt( 
								TMath::Power(CurrentDirtError / EffectiveEfficiencyXBin,2.)  +
								TMath::Power(CurrentDirtEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2.)
							   ) * ScalingFactor;
					DirtScaledErrorUnfOnly = CurrentDirtEntry * EffectiveEfficiencyXBinError / (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin) * ScalingFactor;
							   
				}

				PlotsReco[3][WhichPlot]->SetBinContent(WhichXBin+1,DirtScaledEntry);
				PlotsReco[3][WhichPlot]->SetBinError(WhichXBin+1,DirtScaledError);

				// ----------------------------------------------------------------------------------------------------------------

				// Beam Related Overlay Bkg

				double CurrentBkgEntry = PlotsBkgReco[0][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentBkgError = PlotsBkgReco[0][WhichPlot]->GetBinError(WhichXBin+1);
				double BkgScaledEntry = 0., BkgScaledError = 0., BkgScaledErrorUnfOnly = 0.;

				if (EffectiveEfficiencyXBin != 0 ) {
				
					BkgScaledEntry = CurrentBkgEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					BkgScaledError = sqrt( 
								TMath::Power(CurrentBkgError / EffectiveEfficiencyXBin,2.)  +
								TMath::Power(CurrentBkgEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2.)
							  ) * ScalingFactor;
					BkgScaledErrorUnfOnly = CurrentBkgEntry * EffectiveEfficiencyXBinError / (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin) * ScalingFactor;

				}

				PlotsBkgReco[0][WhichPlot]->SetBinContent(WhichXBin+1,BkgScaledEntry);
				PlotsBkgReco[0][WhichPlot]->SetBinError(WhichXBin+1,BkgScaledError);

				// ------------------------------------------------------------------------------------------------------------------

				// Overlay

				double CurrentOverlayEntry = PlotsCC1pReco[0][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentOverlayError = PlotsCC1pReco[0][WhichPlot]->GetBinError(WhichXBin+1);
				double OverlayScaledEntry = 0., OverlayScaledError = 0., OverlayScaledErrorUnfOnly = 0.;

				if (EffectiveEfficiencyXBin != 0 ) {
				
					OverlayScaledEntry = CurrentOverlayEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					OverlayScaledError = sqrt( 
								TMath::Power(CurrentOverlayError / EffectiveEfficiencyXBin,2.)  +
								TMath::Power(CurrentOverlayEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2.)
							      ) * ScalingFactor;
					OverlayScaledErrorUnfOnly = CurrentOverlayEntry * EffectiveEfficiencyXBinError / (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin);

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
				double DataScaledEntry = 0., DataScaledError = 0., DataScaledErrorUnfOnly = 0.;

				if (EffectiveEfficiencyXBin != 0 ) {  
				
					DataScaledEntry = CurrentDataEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					DataScaledError = sqrt( 
								TMath::Power(CurrentDataError / EffectiveEfficiencyXBin,2.)  +
								TMath::Power(CurrentDataEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2.)
							   ) * ScalingFactor;
					DataScaledErrorUnfOnly = CurrentDataEntry * EffectiveEfficiencyXBinError / (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin) * ScalingFactor;
							   
				}

				double TotalDataScaledEntry = DataScaledEntry - ExtBNBScaledEntry - BkgScaledEntry - DirtScaledEntry; 
				if (Subtract == "_BUnsubtracted") { TotalDataScaledEntry = DataScaledEntry - ExtBNBScaledEntry - DirtScaledEntry; }

				double TotalDataScaledError = TMath::Sqrt( 
										TMath::Power(DataScaledError,2.) +
										TMath::Power(ExtBNBScaledError,2.) +
										TMath::Power(BkgScaledError,2.) +
										TMath::Power(DirtScaledError,2.)						
									   );	

				double TotalDataScaledErrorUnfOnly = TMath::Sqrt( 
										TMath::Power(DataScaledErrorUnfOnly,2.) +
										TMath::Power(ExtBNBScaledErrorUnfOnly,2.) +
										TMath::Power(BkgScaledErrorUnfOnly,2.) +
										TMath::Power(DirtScaledErrorUnfOnly,2.)						
									   );
									   
				if (Subtract == "_BUnsubtracted") { TotalDataScaledError = TMath::Sqrt( 
													TMath::Power(DataScaledError,2.) +
													TMath::Power(ExtBNBScaledError,2.) +
													TMath::Power(DirtScaledError,2.)						
									   				); 
								  }
									   
				PlotsReco[1][WhichPlot]->SetBinContent(WhichXBin+1,TotalDataScaledEntry);
				PlotsReco[1][WhichPlot]->SetBinError(WhichXBin+1,TotalDataScaledError);

				PlotsRecoUnfOnly[WhichPlot]->SetBinContent(WhichXBin+1,TotalDataScaledEntry);
				PlotsRecoUnfOnly[WhichPlot]->SetBinError(WhichXBin+1,TotalDataScaledErrorUnfOnly);

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
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);
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

				TString Label = ToStringPOT(DataPOT)+" POT";

				TLegendEntry* lMC = leg->AddEntry(PlotsCC1pReco[0][WhichPlot],"Unfolded MC","f");
				lMC->SetTextColor(OverlayColor);

				TLegendEntry* lGenie = leg->AddEntry(PlotsTrue[4][WhichPlot],"True MC","l");			
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

		} // End of the loop over the plots
					
		// --------------------------------------------------------------------------------------------------------------------------------
		
		FileEfficiences->Close();

		// --------------------------------------------------------------------------------------------------------------------------------

		// Store the extracted xsections

		TString NameExtractedXSec = PathToExtractedXSec+"ExtractedXSec_"+NameOfSamples[0]+"_"\
			+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+Subtract+".root";
		TFile* ExtractedXSec = TFile::Open(NameExtractedXSec,"recreate");

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

			PlotsReco[1][WhichPlot]->Write(); // Data Reco
			PlotsRecoUnfOnly[WhichPlot]->Write("UnfOnlyReco"+PlotNames[WhichPlot]); // Only with uncertainties due to unfolding procedure
			PlotsCC1pReco[0][WhichPlot]->Write(); // Overlay MC		
			PlotsTrue[4][WhichPlot]->Write(); // Genie Overlay // Closure test			

		}

		ExtractedXSec->Close();

		std::cout << std::endl << "File " << NameExtractedXSec << " created" << std::endl << std::endl;

		for (int i = 0; i < NSamples; i++) { FileSample[i]->Close(); }

	} // End of the loop over the runs	

	FluxFile->Close();

} // End of the program 
