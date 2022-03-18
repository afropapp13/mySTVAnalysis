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

//----------------------------------------//

TString ToStringPOT(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

//----------------------------------------//

void FakeData_XSection_Extraction(TString OverlaySample = "Overlay9",TString FakeDataSample = "Overlay9NuWro") {

	//----------------------------------------//

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	int DecimalAccuracy = 2;

	//----------------------------------------//

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

	//----------------------------------------//

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

	//----------------------------------------//

	gStyle->SetPalette(55); const Int_t NCont = 999; 
	gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");

	vector<TString> LabelsOfSamples;
	vector<TString> NameOfSamples;

	NameOfSamples.push_back(OverlaySample); int GenieRecoCC1pIndex = 0;	
	NameOfSamples.push_back(FakeDataSample); int AltGenRecoCC1pIndex = 1;	
	NameOfSamples.push_back("GenieOverlay"); int GenieTruthIndex = 2;
	NameOfSamples.push_back("AltEventGen");	int AltGenTruthIndex = 3;	
	
	int DataIndex = -1.;

	//----------------------------------------//

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");				
	Runs.push_back("Combined");			

	int NRuns = (int)(Runs.size());
	//cout << "Number of Runs = " << NRuns << endl;

	//----------------------------------------//
	
	// CV Flux File

	TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));

	//----------------------------------------//

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		//----------------------------------------//
	
		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);													
		
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface);	
		//cout << Runs[WhichRun] << " integrated flux = " << IntegratedFlux << endl;		

		//----------------------------------------//	

		vector<TCanvas*> PlotCanvas; PlotCanvas.clear();

		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();
		vector<vector<TH1D*> > PlotsTEfficiency; PlotsTEfficiency.clear();

		//----------------------------------------//

		TString FileEfficienciesName = "FileEfficiences_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TFile* FileEfficiences = new TFile(FileEfficienciesPath+FileEfficienciesName,"readonly");

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();
		
		TString PathToFilesUBCodeExtension = PathToFiles+CutExtension;

		//----------------------------------------//		

		// Loop over the samples

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			//----------------------------------------//
			
			if (NameOfSamples[WhichSample] == "Overlay9" || NameOfSamples[WhichSample] == "Overlay9NuWro") { 
			
				TString FileName = "STVStudies_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+CutExtension+".root";
				FileSample.push_back(TFile::Open(PathToFilesUBCodeExtension+"/"+FileName)); 
				
			}
			
			if (NameOfSamples[WhichSample] == "GenieOverlay") { 
			
				TString FileName = "TruthSTVAnalysis_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
				FileSample.push_back(TFile::Open(PathToFiles+FileName));  
				
			}	

			if (NameOfSamples[WhichSample] == "AltEventGen") { 
			
				TString FileName = "TruthSTVAnalysis_"+FakeDataSample+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";	
				if (FakeDataSample == "Overlay9NuWro") { FileName = "TruthSTVAnalysis_Overlay9NuWro_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }
				FileSample.push_back(TFile::Open(PathToFiles+FileName));  
			
			}									

			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();
			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();
			vector<TH1D*> CurrentPlotsCC1pReco; CurrentPlotsCC1pReco.clear();
			vector<TH1D*> CurrentPlotsTEfficiency; CurrentPlotsTEfficiency.clear();

			//----------------------------------------//			

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
				CurrentPlotsReco.push_back(histReco);

				TH1D* histCC1pReco = (TH1D*)(FileSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]));
				CurrentPlotsCC1pReco.push_back(histCC1pReco);

				TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
				CurrentPlotsTrue.push_back(histTrue);

				TH1D* histTEfficiency = (TH1D*)(FileEfficiences->Get("CC1pReco"+PlotNames[WhichPlot]));
				CurrentPlotsTEfficiency.push_back(histTEfficiency);
		
			} // End of the loop over the plots

			PlotsReco.push_back(CurrentPlotsReco);		
			PlotsTrue.push_back(CurrentPlotsTrue);		
			PlotsCC1pReco.push_back(CurrentPlotsCC1pReco);
			PlotsTEfficiency.push_back(CurrentPlotsTEfficiency);

		} // End of the loop over the samples

		//----------------------------------------//

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {			

			int NBinsX = PlotsCC1pReco[0][WhichPlot]->GetNbinsX();

			//----------------------------------------//

			// Apply the relevant weights

			for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

				double BinWidth = PlotsCC1pReco[0][WhichPlot]->GetBinWidth(WhichXBin+1);

				// We want the number of events, as if the bin width is 1
				if (PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot") { BinWidth = 1.; }

				//----------------------------------------//

				double ScalingFactor = Units / (IntegratedFlux * NTargets * BinWidth);

				//----------------------------------------//
		
				// Effective Efficiencies
				
				double EffectiveEfficiencyXBin = PlotsTEfficiency[0][WhichPlot]->GetBinContent(WhichXBin+1);
				double EffectiveEfficiencyXBinError = PlotsTEfficiency[0][WhichPlot]->GetBinError(WhichXBin+1);

				if (EffectiveEfficiencyXBin == 0) { cout << "Bin with 0 efficiency" << endl; }

				//----------------------------------------//

				// Genie CC1p Reco			

				double CurrentOverlayEntry = PlotsCC1pReco[GenieRecoCC1pIndex][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentOverlayError = PlotsCC1pReco[GenieRecoCC1pIndex][WhichPlot]->GetBinError(WhichXBin+1);
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

				PlotsCC1pReco[GenieRecoCC1pIndex][WhichPlot]->SetBinContent(WhichXBin+1,OverlayScaledEntry);
				PlotsCC1pReco[GenieRecoCC1pIndex][WhichPlot]->SetBinError(WhichXBin+1,OverlayScaledError);
				
				//----------------------------------------//
				
				// Genie True CC1p

				double CurrentGenieOverlayEntry = PlotsTrue[GenieTruthIndex][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentGenieOverlayError = PlotsTrue[GenieTruthIndex][WhichPlot]->GetBinError(WhichXBin+1);

				double GenieOverlayScaledEntry = 0., GenieOverlayScaledError = 0.;
				GenieOverlayScaledEntry = CurrentGenieOverlayEntry * ScalingFactor; 
				GenieOverlayScaledError = CurrentGenieOverlayError * ScalingFactor; 

				PlotsTrue[GenieTruthIndex][WhichPlot]->SetBinContent(WhichXBin+1,GenieOverlayScaledEntry);
				PlotsTrue[GenieTruthIndex][WhichPlot]->SetBinError(WhichXBin+1,GenieOverlayScaledError);
				
				//----------------------------------------//

				// Reco CC1p Alternative Generator // Fake Data				

				double CurrentDataEntry = PlotsCC1pReco[AltGenRecoCC1pIndex][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentDataError = PlotsCC1pReco[AltGenRecoCC1pIndex][WhichPlot]->GetBinError(WhichXBin+1);
				double DataScaledEntry = 0., DataScaledError = 0., DataScaledErrorUnfOnly = 0.;

				if (EffectiveEfficiencyXBin != 0 ) {  
				
					DataScaledEntry = CurrentDataEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					DataScaledError = sqrt( 
								TMath::Power(CurrentDataError / EffectiveEfficiencyXBin,2.)  +
								TMath::Power(CurrentDataEntry * EffectiveEfficiencyXBinError 
								/ (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2.)
							   ) * ScalingFactor;
							   
				}
									   
				PlotsCC1pReco[AltGenRecoCC1pIndex][WhichPlot]->SetBinContent(WhichXBin+1,DataScaledEntry);
				PlotsCC1pReco[AltGenRecoCC1pIndex][WhichPlot]->SetBinError(WhichXBin+1,DataScaledError);

				//----------------------------------------//
				
				// Alternative Generator truth level prediction

				double CurrentAltGenEntry = PlotsTrue[AltGenTruthIndex][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentAltGenError = PlotsTrue[AltGenTruthIndex][WhichPlot]->GetBinError(WhichXBin+1);

				double AltGenScaledEntry = CurrentAltGenEntry * ScalingFactor; 
				double AltGenScaledError = CurrentAltGenError * ScalingFactor; 

				PlotsTrue[AltGenTruthIndex][WhichPlot]->SetBinContent(WhichXBin+1,AltGenScaledEntry);
				PlotsTrue[AltGenTruthIndex][WhichPlot]->SetBinError(WhichXBin+1,AltGenScaledError);				

			} // End of the loop over the bins
			
			//----------------------------------------//

			TString CanvasName = PlotNames[WhichPlot]+"_"+Runs[WhichRun]+OverlaySample;
			TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.14);
			PlotCanvas->SetLeftMargin(0.17);			

			TLegend* leg = new TLegend(0.2,0.91,0.95,0.99);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.05);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(4);
			leg->SetColumnSeparation(0.33);	

			//----------------------------------------//								

			// Plotting the xsections

			// Truth AltGen

			PlotsTrue[AltGenTruthIndex][WhichPlot]->SetLineColor(OverlayColor);
			PlotsTrue[AltGenTruthIndex][WhichPlot]->SetLineWidth(3);			

			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetXaxis()->CenterTitle();
			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetXaxis()->SetTitleOffset(1.05);
			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetXaxis()->SetTitleSize(0.06);
			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetXaxis()->SetLabelSize(0.06);		
			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetXaxis()->SetNdivisions(5);			

			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetYaxis()->CenterTitle();
			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetYaxis()->SetTitleOffset(1.18);
			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);
			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetYaxis()->SetTitleSize(0.06);
			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetYaxis()->SetLabelSize(0.06);
			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetYaxis()->SetNdivisions(5);

			double max = 1.05*TMath::Max(PlotsTrue[AltGenTruthIndex][WhichPlot]->GetMaximum(),PlotsCC1pReco[AltGenRecoCC1pIndex][WhichPlot]->GetMaximum());
			PlotsTrue[AltGenTruthIndex][WhichPlot]->GetYaxis()->SetRangeUser(0.,max);

			PlotsTrue[AltGenTruthIndex][WhichPlot]->Draw("hist same");
			leg->AddEntry(PlotsTrue[AltGenTruthIndex][WhichPlot],"True NuWro","l");

			//----------------------------------------//	

			PlotsCC1pReco[AltGenRecoCC1pIndex][WhichPlot]->SetLineWidth(3);
			PlotsCC1pReco[AltGenRecoCC1pIndex][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsCC1pReco[AltGenRecoCC1pIndex][WhichPlot]->SetMarkerColor(BeamOnColor);
			PlotsCC1pReco[AltGenRecoCC1pIndex][WhichPlot]->SetMarkerStyle(20);
			PlotsCC1pReco[AltGenRecoCC1pIndex][WhichPlot]->SetMarkerSize(1.5);

			PlotsCC1pReco[AltGenRecoCC1pIndex][WhichPlot]->Draw("p hist same");		
			leg->AddEntry(PlotsCC1pReco[AltGenRecoCC1pIndex][WhichPlot],"Unfolded CC1p NuWro","p");			

			//----------------------------------------//

			leg->Draw();	

			TLatex *text = new TLatex();
			text->SetTextFont(FontStyle);
			text->SetTextSize(0.05);
			text->DrawTextNDC(0.21, 0.86,"Effective Efficiencies");			
			
			TString CanvasPath = PlotPath+NameOfSamples[0];
			TString FullCanvasName = "/"+FakeDataSample+"FakeData_XSections_"+CanvasName+"_"+UBCodeVersion+".pdf";
			PlotCanvas->SaveAs(CanvasPath+FullCanvasName);

			delete PlotCanvas;
			
			// --------------------------------------------------------------------------------------------------------------------------------

		} // End of the loop over the plots
					
		// --------------------------------------------------------------------------------------------------------------------------------
		
		FileEfficiences->Close();

	} // End of the loop over the runs	

	FluxFile->Close();

} // End of the program 
