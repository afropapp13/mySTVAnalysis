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

//#include "/home/afroditi/Dropbox/PhD/Secondary_Code/myFunctions.cpp"

#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

TString ToString(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

void XSection_Extraction(TString OverlaySample,int Universe = -1) {

	// -------------------------------------------------------------------------------------

	TString PathToFiles = "/uboone/data/users/"+UserID+"/myEvents/OutputFiles/"+UBCodeVersion+"/";
	TString PathToExtractedXSec = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myXSec/"+UBCodeVersion+"/"; 
	TString FileEfficienciesPath = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myEfficiencies/"+UBCodeVersion+"/";
	TString PlotPath = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myPlots/"+UBCodeVersion+"/"; 

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
	VectorCuts.push_back("");
	VectorCuts.push_back("_NuScore");
	VectorCuts.push_back("_ThreePlaneLogChi2");
	VectorCuts.push_back("_Collinearity");

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
	PlotNames.push_back("ECalPlot");
	PlotNames.push_back("EQEPlot"); 
	PlotNames.push_back("Q2Plot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");

	vector<TString> LabelsOfSamples;
	vector<TString> NameOfSamples;

	NameOfSamples.push_back("Overlay9"); 
	NameOfSamples.push_back("BeamOn9"); 
	NameOfSamples.push_back("ExtBNB9"); 
	NameOfSamples.push_back("OverlayDirt9"); 
//	NameOfSamples.push_back("Genie");
	NameOfSamples.push_back("GenieOverlay");	
//	NameOfSamples.push_back("GiBUU");	
//	NameOfSamples.push_back("NuWro");	
//	NameOfSamples.push_back("NEUT");	
//	NameOfSamples.push_back("GENIEv2");	
//	NameOfSamples.push_back("SuSav2");	
	
	int DataIndex = -1.;

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
	
	// CV Flux File

	TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));
	double DataPOT = -99.;
	
	if ( Universe != -1 ) {
	
		TString DublicateOverlaySample = OverlaySample;
		TString ReducedOverlaySample = DublicateOverlaySample.ReplaceAll("m_","m");
		if ( !(string(OverlaySample).find("expskin") != std::string::npos) ) { ReducedOverlaySample = ReducedOverlaySample.ReplaceAll("n_","n"); }
		ReducedOverlaySample = ReducedOverlaySample.ReplaceAll("g_","g");				
		
		for (int i = 0; i < 10;i++) { ReducedOverlaySample.ReplaceAll(TString(std::to_string(i)),""); }
		
		TString FluxHistoName = "numu_ms"+ReducedOverlaySample+"/hEnumu"+ReducedOverlaySample+"_ms_"+TString(std::to_string(Universe));
		HistoFlux = (TH1D*)(FluxFile->Get(FluxHistoName));
	
	}

	// -------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {
	
		// -------------------------------------------------------------------------------------			
	
		if (Runs[WhichRun] == "Run1") { DataPOT = tor860_wcut_Run1 ; }
//		if (Runs[WhichRun] == "Run2") { DataPOT = tor860_wcut_Run2 ; }
		if (Runs[WhichRun] == "Run3") { DataPOT = tor860_wcut_Run3 ; }
//		if (Runs[WhichRun] == "Run4") { DataPOT = tor860_wcut_Run4 ; }
//		if (Runs[WhichRun] == "Run5") { DataPOT = tor860_wcut_Run5 ; }								
		
			
		// Zarko Pavlovic, Jun 22 2020
		double IntegratedFlux = HistoFlux->Integral() * DataPOT / (4997.*5e8) / (256.35*233.);			

		// -------------------------------------------------------------------------------------		

		vector<TCanvas*> PlotCanvas; PlotCanvas.clear();

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
			
//			if (NameOfSamples[WhichSample] == "Genie" || NameOfSamples[WhichSample] == "GiBUU" 
//				|| NameOfSamples[WhichSample] == "NuWro" || NameOfSamples[WhichSample] == "NEUT" 
//				|| NameOfSamples[WhichSample] == "GENIEv2" || NameOfSamples[WhichSample] == "SuSav2") { 
//			
//				TString FileName = "myFiles/"+UBCodeVersion+"/STVAnalysis_"+NameOfSamples[WhichSample]+"_"+UBCodeVersion+".root";	
//				FileSample.push_back(TFile::Open(FileName)); 
//				
//			}
			
			if (NameOfSamples[WhichSample] == "GenieOverlay") { 
			
				TString FileName = "TruthSTVAnalysis_Overlay9_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
				FileSample.push_back(TFile::Open(PathToFiles+FileName));  
				
			}			

			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();
			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();
			vector<TH1D*> CurrentPlotsBkgReco; CurrentPlotsBkgReco.clear();
			vector<TH1D*> CurrentPlotsCC1pReco; CurrentPlotsCC1pReco.clear();
			vector<TH1D*> CurrentPlotsTEfficiency; CurrentPlotsTEfficiency.clear();

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
		
			}

			PlotsReco.push_back(CurrentPlotsReco);		
			PlotsTrue.push_back(CurrentPlotsTrue);		
			PlotsBkgReco.push_back(CurrentPlotsBkgReco);
			PlotsCC1pReco.push_back(CurrentPlotsCC1pReco);
			PlotsTEfficiency.push_back(CurrentPlotsTEfficiency);


		}

		// ----------------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {
	
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
//			leg->SetNColumns(3);			
//			leg->SetNColumns(2);
			leg->SetColumnSeparation(0.33);			
			
			TLegend* legData = new TLegend(0.15,0.84,0.9,0.88);
			legData->SetBorderSize(0);
			legData->SetTextSize(0.06);
			legData->SetTextFont(FontStyle);
			legData->SetNColumns(3);			

			int NBinsX = PlotsCC1pReco[0][WhichPlot]->GetNbinsX();

			// ---------------------------------------------------------------------------------------------------------------------------

			// Apply the relevant weights

			for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

				double BinWidth = PlotsCC1pReco[0][WhichPlot]->GetBinWidth(WhichXBin+1);

				// ------------------------------------------------------------------------------------------------------------------

				double FluxWeight = 1.;

				// New Integrated Flux			

// Fix it !!!!
//				if (string(fWhichSample).find("FluxUnisim") != std::string::npos || string(fWhichSample).find("Primary") != std::string::npos) {
//
//					double POT = -99.;
//					if (string(OverlaySample).find("Run1") != std::string::npos) { POT = tor860_wcut_Run1 ; }
//
//					TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root");
//					TH1D* FluxHisto = (TH1D*)(FluxFile->Get("numu_ms_"+fWhichSample+"/hEnumu_"+fWhichSample+"_"+fUniverse)); // have to specify the universe
//					double FluxIntegral = FluxHisto->Integral() * POT / (2.43e11 * 256.35 * 233.);
//
//					FluxWeight = Flux / FluxIntegral;
//				}

				// -----------------------------------------------------------------------------------------------------------------

				double ScalingFactor = (Units / (IntegratedFlux * NTargets * BinWidth)) * FluxWeight;

				// -----------------------------------------------------------------------------------------------------------------
		
				// Effective Efficiencies
				
				double EffectiveEfficiencyXBin = PlotsTEfficiency[0][WhichPlot]->GetBinContent(WhichXBin+1);
				double EffectiveEfficiencyXBinError = PlotsTEfficiency[0][WhichPlot]->GetBinError(WhichXBin+1);

				// -----------------------------------------------------------------------------------------------------------------

				// ExtBNB

				double CurrentExtBNBEntry = PlotsReco[2][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentExtBNBError = PlotsReco[2][WhichPlot]->GetBinError(WhichXBin+1);
				double ExtBNBScaledEntry = 0., ExtBNBScaledError = 0.;
				
				if (EffectiveEfficiencyXBin != 0 ) {
				  
					ExtBNBScaledEntry = CurrentExtBNBEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					ExtBNBScaledError = sqrt( 
								TMath::Power(CurrentExtBNBError / EffectiveEfficiencyXBin,2)  +
								TMath::Power(CurrentExtBNBEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2)
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
								TMath::Power(CurrentDirtError / EffectiveEfficiencyXBin,2)  +
								TMath::Power(CurrentDirtEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2)
							   ) * ScalingFactor;
							   
				}

				PlotsReco[3][WhichPlot]->SetBinContent(WhichXBin+1,DirtScaledEntry);
				PlotsReco[3][WhichPlot]->SetBinError(WhichXBin+1,DirtScaledError);

				// ------------------------------------------------------------------------------------------------------------------

				// Overlay

				double CurrentOverlayEntry = PlotsCC1pReco[0][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentOverlayError = PlotsCC1pReco[0][WhichPlot]->GetBinError(WhichXBin+1);
				double OverlayScaledEntry = 0., OverlayScaledError = 0.;
				if (EffectiveEfficiencyXBin != 0 ) {
				
					OverlayScaledEntry = CurrentOverlayEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					OverlayScaledError = sqrt( 
								TMath::Power(CurrentOverlayError / EffectiveEfficiencyXBin,2)  +
								TMath::Power(CurrentOverlayEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2)
							      ) * ScalingFactor;

				}

				PlotsCC1pReco[0][WhichPlot]->SetBinContent(WhichXBin+1,OverlayScaledEntry);
				PlotsCC1pReco[0][WhichPlot]->SetBinError(WhichXBin+1,OverlayScaledError);

				// ----------------------------------------------------------------------------------------------------------------

				// Bkg

				double CurrentBkgEntry = PlotsBkgReco[0][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentBkgError = PlotsBkgReco[0][WhichPlot]->GetBinError(WhichXBin+1);
				double BkgScaledEntry = 0., BkgScaledError = 0.;
				if (EffectiveEfficiencyXBin != 0 ) {
				
					BkgScaledEntry = CurrentBkgEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					BkgScaledError = sqrt( 
								TMath::Power(CurrentBkgError / EffectiveEfficiencyXBin,2)  +
								TMath::Power(CurrentBkgEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2)
							  ) * ScalingFactor;

				}

				PlotsBkgReco[0][WhichPlot]->SetBinContent(WhichXBin+1,BkgScaledEntry);
				PlotsBkgReco[0][WhichPlot]->SetBinError(WhichXBin+1,BkgScaledError);

				// --------------------------------------------------------------------------------------------------------

//				// Genie v3.0.6 Out of the box

//				double CurrentGenieEntry = PlotsTrue[4][WhichPlot]->GetBinContent(WhichXBin+1);
//				double CurrentGenieError = PlotsTrue[4][WhichPlot]->GetBinError(WhichXBin+1);
//				double GenieScaledEntry = 0., GenieScaledError = 0.;

//				GenieScaledEntry = FluxIntegratedXSection * CurrentGenieEntry / BinWidth ; 
//				GenieScaledError = FluxIntegratedXSection * CurrentGenieError / BinWidth ; 

//				PlotsTrue[4][WhichPlot]->SetBinContent(WhichXBin+1,GenieScaledEntry);
//				PlotsTrue[4][WhichPlot]->SetBinError(WhichXBin+1,GenieScaledError);
				
				// --------------------------------------------------------------------------------------------------------------
				
				// Genie Overlay
			
//				double CurrentGenieOverlayEntry = PlotsTrue[5][WhichPlot]->GetBinContent(WhichXBin+1);
//				double CurrentGenieOverlayError = PlotsTrue[5][WhichPlot]->GetBinError(WhichXBin+1);

				double CurrentGenieOverlayEntry = PlotsTrue[4][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentGenieOverlayError = PlotsTrue[4][WhichPlot]->GetBinError(WhichXBin+1);

				double GenieOverlayScaledEntry = 0., GenieOverlayScaledError = 0.;
				GenieOverlayScaledEntry = CurrentGenieOverlayEntry * ScalingFactor; 
				GenieOverlayScaledError = CurrentGenieOverlayError * ScalingFactor; 

//				PlotsTrue[5][WhichPlot]->SetBinContent(WhichXBin+1,GenieOverlayScaledEntry);
//				PlotsTrue[5][WhichPlot]->SetBinError(WhichXBin+1,GenieOverlayScaledError);

				PlotsTrue[4][WhichPlot]->SetBinContent(WhichXBin+1,GenieOverlayScaledEntry);
				PlotsTrue[4][WhichPlot]->SetBinError(WhichXBin+1,GenieOverlayScaledError);
				
				// --------------------------------------------------------------------------------------------------------

//				// GiBUU

//				double CurrentGiBUUEntry = PlotsTrue[6][WhichPlot]->GetBinContent(WhichXBin+1);
//				double CurrentGiBUUError = PlotsTrue[6][WhichPlot]->GetBinError(WhichXBin+1);
//				double GiBUUScaledEntry = 0., GiBUUScaledError = 0.;

//				// All the scaling factors have been included
//				
//				GiBUUScaledEntry = CurrentGiBUUEntry / BinWidth ; 
////				GiBUUScaledError = CurrentGiBUUError / BinWidth ; 
//				GiBUUScaledError = 0.00001 ; 

//				PlotsTrue[6][WhichPlot]->SetBinContent(WhichXBin+1,GiBUUScaledEntry);
//				PlotsTrue[6][WhichPlot]->SetBinError(WhichXBin+1,GiBUUScaledError);				

//				// --------------------------------------------------------------------------------------------------------

//				// NuWro

//				double CurrentNuWroEntry = PlotsTrue[7][WhichPlot]->GetBinContent(WhichXBin+1);
//				double CurrentNuWroError = PlotsTrue[7][WhichPlot]->GetBinError(WhichXBin+1);
//				double NuWroScaledEntry = 0., NuWroScaledError = 0.;

//				// All the scaling factors have been included
//				
//				NuWroScaledEntry = CurrentNuWroEntry / BinWidth ; 
////				NuWroScaledError = CurrentNuWroError / BinWidth ; 
//				NuWroScaledError = 0.00001 ; 

//				PlotsTrue[7][WhichPlot]->SetBinContent(WhichXBin+1,NuWroScaledEntry);
//				PlotsTrue[7][WhichPlot]->SetBinError(WhichXBin+1,NuWroScaledError);	
//				
//				// --------------------------------------------------------------------------------------------------------

//				// NEUT

//				double CurrentNEUTEntry = PlotsTrue[8][WhichPlot]->GetBinContent(WhichXBin+1);
//				double CurrentNEUTError = PlotsTrue[8][WhichPlot]->GetBinError(WhichXBin+1);
//				double NEUTScaledEntry = 0., NEUTScaledError = 0.;

//				// All the scaling factors have been included
//				
//				NEUTScaledEntry = CurrentNEUTEntry / BinWidth ; 
////				NEUTScaledError = CurrentNEUTError / BinWidth ; 
//				NEUTScaledError = 0.00001 ; 

//				PlotsTrue[8][WhichPlot]->SetBinContent(WhichXBin+1,NEUTScaledEntry);
//				PlotsTrue[8][WhichPlot]->SetBinError(WhichXBin+1,NEUTScaledError);	
//				
//				// --------------------------------------------------------------------------------------------------------

//				// GENIEv2

//				double CurrentGENIEv2Entry = PlotsTrue[9][WhichPlot]->GetBinContent(WhichXBin+1);
//				double CurrentGENIEv2Error = PlotsTrue[9][WhichPlot]->GetBinError(WhichXBin+1);
//				double GENIEv2ScaledEntry = 0., GENIEv2ScaledError = 0.;

//				// All the scaling factors have been included
//				
//				GENIEv2ScaledEntry = CurrentGENIEv2Entry / BinWidth ; 
//if (PlotNames[WhichPlot] == "DeltaPTPlot"  && WhichXBin == 0) { GENIEv2ScaledEntry = 0.8*NEUTScaledEntry; }		
////				GENIEv2ScaledError = CurrentGENIEv2Error / BinWidth ; 
//				GENIEv2ScaledError = 0.00001 ; 

//				PlotsTrue[9][WhichPlot]->SetBinContent(WhichXBin+1,GENIEv2ScaledEntry);
//				PlotsTrue[9][WhichPlot]->SetBinError(WhichXBin+1,GENIEv2ScaledError);
//				
//				// --------------------------------------------------------------------------------------------------------

//				// SuSav2

//				double CurrentSuSav2Entry = PlotsTrue[10][WhichPlot]->GetBinContent(WhichXBin+1);
//				double CurrentSuSav2Error = PlotsTrue[10][WhichPlot]->GetBinError(WhichXBin+1);
//				double SuSav2ScaledEntry = 0., SuSav2ScaledError = 0.;

//				// All the scaling factors have been included
//				
//				SuSav2ScaledEntry = SuSav2FluxIntegratedXSection * CurrentSuSav2Entry / BinWidth ; 
////				SuSav2ScaledError = SuSav2FluxIntegratedXSection * CurrentSuSav2Error / BinWidth ; 
//				SuSav2ScaledError = 0.00001 ; 

//				PlotsTrue[10][WhichPlot]->SetBinContent(WhichXBin+1,SuSav2ScaledEntry);
//				PlotsTrue[10][WhichPlot]->SetBinError(WhichXBin+1,SuSav2ScaledError);							
				
				// --------------------------------------------------------------------------------------------------------------------

				// BeamOn

				double CurrentDataEntry = PlotsReco[1][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentDataError = PlotsReco[1][WhichPlot]->GetBinError(WhichXBin+1);
				double DataScaledEntry = 0., DataScaledError = 0.;
				if (EffectiveEfficiencyXBin != 0 ) {  
				
					DataScaledEntry = CurrentDataEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					DataScaledError = sqrt( 
								TMath::Power(CurrentDataError / EffectiveEfficiencyXBin,2)  +
								TMath::Power(CurrentDataEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2)
							   ) * ScalingFactor;
							   
				}

//				PlotsReco[1][WhichPlot]->SetBinContent(WhichXBin+1,DataScaledEntry);
//				PlotsReco[1][WhichPlot]->SetBinError(WhichXBin+1,DataScaledError);


				double TotalDataScaledEntry = DataScaledEntry - ExtBNBScaledEntry - BkgScaledEntry - DirtScaledEntry; 
				double TotalDataScaledError = TMath::Sqrt( 
										TMath::Power(DataScaledError,2.) +
										TMath::Power(ExtBNBScaledError,2.) +
										TMath::Power(BkgScaledError,2.) +
										TMath::Power(DirtScaledError,2.)						
									   );	
									   
									   
				PlotsReco[1][WhichPlot]->SetBinContent(WhichXBin+1,TotalDataScaledEntry);
				PlotsReco[1][WhichPlot]->SetBinError(WhichXBin+1,TotalDataScaledError);

				// -------------------------------------------------------------------------------------				

				// Samples for systematics

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

			midPad->cd();

			// BeamOn

			PlotsReco[1][WhichPlot]->SetLineWidth(3);
			PlotsReco[1][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsReco[1][WhichPlot]->SetMarkerColor(BeamOnColor);
			PlotsReco[1][WhichPlot]->SetMarkerStyle(20);
			PlotsReco[1][WhichPlot]->SetMarkerSize(1.5);
//			// Subtract ExtBNB
////			PlotsReco[1][WhichPlot]->Add(PlotsReco[2][WhichPlot],-1); 
//			// Subtract Dirt
////			PlotsReco[1][WhichPlot]->Add(PlotsReco[3][WhichPlot],-1); 
//			// Subtract MC Bkg Subtraction
////			if (Subtract == "") { PlotsReco[1][WhichPlot]->Add(PlotsBkgReco[0][WhichPlot],-1); }

			// GENIE v3.0.6

//			PlotsTrue[4][WhichPlot]->SetLineWidth(3);
//			PlotsTrue[4][WhichPlot]->SetLineColor(GenieColor);
//			PlotsTrue[4][WhichPlot]->SetFillColor(GenieColor);
			
			// GENIE Overlay

//			PlotsTrue[5][WhichPlot]->SetLineWidth(3);
////			PlotsTrue[5][WhichPlot]->SetLineColor(GenieOverlayColor);
////			PlotsTrue[5][WhichPlot]->SetFillColor(GenieOverlayColor);	
//			PlotsTrue[5][WhichPlot]->SetLineColor(GenieColor);
//			PlotsTrue[5][WhichPlot]->SetFillColor(GenieColor);

			PlotsTrue[4][WhichPlot]->SetLineWidth(3);
//			PlotsTrue[4][WhichPlot]->SetLineColor(GenieOverlayColor);
//			PlotsTrue[4][WhichPlot]->SetFillColor(GenieOverlayColor);	
			PlotsTrue[4][WhichPlot]->SetLineColor(GenieColor);
			PlotsTrue[4][WhichPlot]->SetFillColor(GenieColor);
			
//			// GiBUU

//			PlotsTrue[6][WhichPlot]->SetLineWidth(3);
//			PlotsTrue[6][WhichPlot]->SetLineColor(GiBUUColor);
//			PlotsTrue[6][WhichPlot]->SetFillColor(GiBUUColor);
//			
//			// NuWro

//			PlotsTrue[7][WhichPlot]->SetLineWidth(3);
//			PlotsTrue[7][WhichPlot]->SetLineColor(NuWroColor);
//			PlotsTrue[7][WhichPlot]->SetFillColor(NuWroColor);	
//			
//			// NEUT

//			PlotsTrue[8][WhichPlot]->SetLineWidth(3);
//			PlotsTrue[8][WhichPlot]->SetLineColor(NEUTColor);
//			PlotsTrue[8][WhichPlot]->SetFillColor(NEUTColor);
//			
//			// GENIEv2

//			PlotsTrue[9][WhichPlot]->SetLineWidth(3);
//			PlotsTrue[9][WhichPlot]->SetLineColor(GENIEv2Color);
//			PlotsTrue[9][WhichPlot]->SetFillColor(GENIEv2Color);
//			
//			// SuSav2

//			PlotsTrue[10][WhichPlot]->SetLineWidth(3);
//			PlotsTrue[10][WhichPlot]->SetLineColor(SuSav2Color);
//			PlotsTrue[10][WhichPlot]->SetFillColor(SuSav2Color);									
			
			// -----------------------------------------------------------------------------------------------------------------------		

			// Max on plot so that we can include the integrated xsecs
			
			double max = TMath::Max(PlotsCC1pReco[0][WhichPlot]->GetMaximum(),PlotsReco[1][WhichPlot]->GetMaximum());
			double YmaxScaleFactor = 1.35;
//			if (PlotNames[WhichPlot] == "MuonMomentumPlot" ) { YmaxScaleFactor = 1.4; }
			double min = TMath::Min(0.,1.5*PlotsReco[1][WhichPlot]->GetMinimum());
//			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetRangeUser(min,1.8 * max);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetRangeUser(min,YmaxScaleFactor * max);

			// Plot MC
			PlotsCC1pReco[0][WhichPlot]->Draw("e2");

			// Plot GENIE v3.0.6
//			PlotsTrue[4][WhichPlot]->Draw("e same");


// Put it back for the closure test			
			// Plot GENIE Overlay
//			PlotsTrue[5][WhichPlot]->Draw("e same");			
			PlotsTrue[4][WhichPlot]->Draw("e same");			
			
//			// Plot GiBUU
//			PlotsTrue[6][WhichPlot]->Draw("e same");
//			
//			// Plot NuWro
//			PlotsTrue[7][WhichPlot]->Draw("e same");
//			
//			// Plot NEUT
//			PlotsTrue[8][WhichPlot]->Draw("e same");
//			
//			// Plot GENIEv2
//			PlotsTrue[9][WhichPlot]->Draw("e same");
//			
//			// Plot SuSav2
//			PlotsTrue[10][WhichPlot]->Draw("e same");											

			// Plot BeamOn
			PlotsReco[1][WhichPlot]->Draw("ex0 same");

			// -------------------------------------------------------------------------------------------------------------------------

			// Integrated cross-sections & chi2 w/o correlations

//			double IntegratedGenieXSection = round(IntegratedXSec(PlotsTrue[4][WhichPlot]),DecimalAccuracy);
//			double IntegratedGenieXSectionError = round(IntegratedXSecError(PlotsTrue[4][WhichPlot]),DecimalAccuracy);

//			double IntegratedOverlayXSection = round(IntegratedXSec(PlotsCC1pReco[0][WhichPlot]),DecimalAccuracy);
//			double IntegratedOverlayXSectionError = round(IntegratedXSecError(PlotsCC1pReco[0][WhichPlot]),DecimalAccuracy);

//			double IntegratedDataXSection = round(IntegratedXSec(PlotsReco[1][WhichPlot]),DecimalAccuracy);
//			double IntegratedDataXSectionError = round(IntegratedXSecError(PlotsReco[1][WhichPlot]),DecimalAccuracy);

//			double chi2 = Chi2(PlotsReco[1][WhichPlot],PlotsCC1pReco[0][WhichPlot]);

//			TString LabelData = "#splitline{#splitline{#color["+TString(std::to_string(BeamOnColor))+"]{#sigma_{Data} = (" 
//					    +TString(std::to_string(IntegratedDataXSection))+" #pm "
//					    +TString(std::to_string(IntegratedDataXSectionError))+") #upoint 10^{-38} cm^{2}}}{#color["
//					    +TString(std::to_string(OverlayColor))+"]{#sigma_{MC} = ("
//					    +TString(std::to_string(IntegratedOverlayXSection))+" #pm "
//					    +TString(std::to_string(IntegratedOverlayXSectionError))+") #upoint 10^{-38} cm^{2}}}}{#color["
//					    +TString(std::to_string(GenieColor))+"]{#sigma_{GENIE} = ("
//					    +TString(std::to_string(IntegratedGenieXSection))+" #pm "
//					    +TString(std::to_string(IntegratedGenieXSectionError))+") #upoint 10^{-38} cm^{2}}}";
					    
//			TString LabelData = "#splitline{#color["+TString(std::to_string(BeamOnColor))+"]{#sigma_{Data} = (" 
//					    +TString(std::to_string(IntegratedDataXSection))+" #pm "
//					    +TString(std::to_string(IntegratedDataXSectionError))+") #upoint 10^{-38} cm^{2}}}{#color["
//					    +TString(std::to_string(OverlayColor))+"]{#sigma_{MC} = ("
//					    +TString(std::to_string(IntegratedOverlayXSection))+" #pm "
//					    +TString(std::to_string(IntegratedOverlayXSectionError))+") #upoint 10^{-38} cm^{2}}}";				    

			TLatex latexSigma;
			latexSigma.SetTextFont(FontStyle);
			latexSigma.SetTextSize(0.06);
////			latexSigma.DrawLatexNDC(0.27,0.67, LabelData);
//			if (PlotNames[WhichPlot] == "MuonMomentumPlot" ) { latexSigma.DrawLatexNDC(0.27,0.74, LabelData); }			
//			if (PlotNames[WhichPlot] == "MuonMomentumPlot" ) { latexSigma.DrawLatexNDC(0.36,0.74, LabelData); }			
			// ----------------------------------------------------------------------------------------------------------------

			// Legend & POT Normalization

			double tor860_wcut = -99.;
			
			if (Runs[WhichRun] == "Run1") { tor860_wcut = tor860_wcut_Run1; }

			TString Label = ToString(tor860_wcut)+" POT";
//			latex.DrawLatexNDC(0.47,0.9, Label);

			TLegendEntry* lMC = leg->AddEntry(PlotsCC1pReco[0][WhichPlot],"MC","f");
			lMC->SetTextColor(OverlayColor);

//			TLegendEntry* lGenie = leg->AddEntry(PlotsTrue[4][WhichPlot],"GENIE v3.0.6","l");
			TLegendEntry* lGenie = leg->AddEntry(PlotsTrue[4][WhichPlot],"GENIE v3","l");			
			lGenie->SetTextColor(GenieColor);

// Put it back for the closure test
/*		
//			TLegendEntry* lGenieOverlay = leg->AddEntry(PlotsTrue[5][WhichPlot],"GENIE Overlay","l");
//			lGenieOverlay->SetTextColor(GenieOverlayColor);			
//			TLegendEntry* lGenieOverlay = leg->AddEntry(PlotsTrue[5][WhichPlot],"GENIE v3.0.6","l");
			TLegendEntry* lGenieOverlay = leg->AddEntry(PlotsTrue[5][WhichPlot],"GENIE v3","l");			
			lGenieOverlay->SetTextColor(GenieColor);
*/
			
//			TLegendEntry* lGENIEv2 = leg->AddEntry(PlotsTrue[9][WhichPlot],"GENIE v2","l");
//			lGENIEv2->SetTextColor(GENIEv2Color);
//			
//			TLegendEntry* lSuSav2 = leg->AddEntry(PlotsTrue[10][WhichPlot],"SuSav2","l");
//			lSuSav2->SetTextColor(SuSav2Color);						
//			
//			TLegendEntry* lGiBUU = leg->AddEntry(PlotsTrue[6][WhichPlot],"GiBUU","l");
//			lGiBUU->SetTextColor(GiBUUColor);
//			
//			TLegendEntry* lNuWro = leg->AddEntry(PlotsTrue[7][WhichPlot],"NuWro","l");
//			lNuWro->SetTextColor(NuWroColor);
//			
//			TLegendEntry* lNEUT = leg->AddEntry(PlotsTrue[8][WhichPlot],"NEUT","l");
//			lNEUT->SetTextColor(NEUTColor);									

			leg->Draw();	

			legData->AddEntry(PlotsReco[1][WhichPlot],"MicroBooNE Data " + Runs[WhichRun] + " " + Label,"ep");
			legData->Draw();
			
			// -----------------------------------------------------------------------------------------------------------------------

			if (OverlaySample == "") {
			
				TString CanvasPath = PlotPath+"/"+NameOfSamples[0];
				TString CanvasName = "/XSections_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+Subtract+".pdf";
				PlotCanvas->SaveAs(CanvasPath+CanvasName);
			}

			if (OverlaySample != "") { delete PlotCanvas; }

			// --------------------------------------------------------------------------------------------------------------------------------
			
			// chi2 calculation
			
			// Only if the data sample is included 
			
//			if (DataIndex != -1 && OverlaySample == "") {
//			 
//			 	// Loop over the samples
//			 	
//				for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {	
//				
//					// Only with respect to the simulation lines
//					
//					if (NameOfSamples[WhichSample] == "Genie" || NameOfSamples[WhichSample] == "GiBUU" 
//					|| NameOfSamples[WhichSample] == "NuWro" || NameOfSamples[WhichSample] == "NEUT" 
//					|| NameOfSamples[WhichSample] == "GENIEv2" || NameOfSamples[WhichSample] == "SuSav2") {
//					
//						int NBinsTrue = PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->GetNbins();
//						double chi2 = Chi2(PlotsReco[DataIndex][WhichPlot],PlotsTrue[WhichSample][WhichPlot]);
//							
//						cout << PlotNames[WhichPlot] << "  " << NameOfSamples[WhichSample] << "  chi2/d.o.f = ";
//						cout << round(chi2,DecimalAccuracy) << "/" << NBinsTrue << endl; 		
//		
//					} // End of the simulation predictions 
//				
//				} // End of the loop over the samples
//				
//				cout << endl << endl;

//			} // End of the case that we have a data sample
			
			// --------------------------------------------------------------------------------------------------------------------------------

		} // End of the loop over the plots
					
		// --------------------------------------------------------------------------------------------------------------------------------
		
		FileEfficiences->Close();

		// --------------------------------------------------------------------------------------------------------------------------------

		// Store the extracted xsections

		TString NameExtractedXSec = PathToExtractedXSec+"ExtractedXSec_"+NameOfSamples[0]+"_"\
			+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
		TFile* ExtractedXSec = TFile::Open(NameExtractedXSec,"recreate");

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

			PlotsReco[1][WhichPlot]->Write();
			PlotsCC1pReco[0][WhichPlot]->Write();
			//PlotsTrue[4][WhichPlot]->Write(); // Genie v3.0.6
//			PlotsTrue[5][WhichPlot]->Write(); // Genie Overlay			
			PlotsTrue[4][WhichPlot]->Write(); // Genie Overlay			

		}

		ExtractedXSec->Close();
		FluxFile->Close();

		std::cout << std::endl << "File " << NameExtractedXSec << " created" << std::endl << std::endl;

	} // End of the loop over the runs	

} // End of the program 
