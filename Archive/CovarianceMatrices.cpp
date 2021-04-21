#include <TFile.h>
#include <TF1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TLatex.h>

#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

void CovarianceMatrices(TString OverlaySample,int Universe = -1,TString BeamOnSample = "BeamOn9",TString BeamOffSample = "ExtBNB9",TString DirtSample = "OverlayDirt9") {

	// -------------------------------------------------------------------------------------

	TH2D::SetDefaultSumw2();
	
	double TextSize = 0.07;
	
	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(TextSize,"t"); 
	gStyle->SetTitleFont(FontStyle,"t");
	gStyle->SetOptStat(0);

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

	vector<TString> PlotNames;
	PlotNames.push_back("DeltaPTPlot"); 
	PlotNames.push_back("DeltaAlphaTPlot"); 
	PlotNames.push_back("DeltaPhiTPlot"); 
	PlotNames.push_back("MuonMomentumPlot"); 
	PlotNames.push_back("MuonPhiPlot"); 
	PlotNames.push_back("MuonCosThetaPlot");
	PlotNames.push_back("ProtonMomentumPlot"); 
	PlotNames.push_back("ProtonPhiPlot"); 
	PlotNames.push_back("ProtonCosThetaPlot");
//	PlotNames.push_back("ECalPlot"); 
//	PlotNames.push_back("EQEPlot"); 
//	PlotNames.push_back("Q2Plot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// -------------------------------------------------------------------------------------------------------------------------------------
	
	vector<TString> NameOfSamples;
	NameOfSamples.push_back("Overlay9");
	const int NSamples = NameOfSamples.size();
		
	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");
//	Runs.push_back("Run2");
	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");				

	const int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;	

	// -------------------------------------------------------------------------------------

	vector< vector<TFile*>> MCFileSample;
	MCFileSample.resize(NSamples, vector<TFile*>(NRuns));

	vector< vector<TFile*>> BeamOnFileSample;
	BeamOnFileSample.resize(NSamples, vector<TFile*>(NRuns));

	vector< vector<TFile*>> BeamOffFileSample;
	BeamOffFileSample.resize(NSamples, vector<TFile*>(NRuns));

	vector< vector<TFile*>> DirtFileSample;
	DirtFileSample.resize(NSamples, vector<TFile*>(NRuns));

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	vector< vector <TH1D*> > CC1pPlots;
	CC1pPlots.resize(NSamples, vector<TH1D*>(N1DPlots));

	vector< vector <TH1D*> > NonCC1pPlots;
	NonCC1pPlots.resize(NSamples, vector<TH1D*>(N1DPlots));

	vector< vector <TH1D*> > BeamOnPlots;
	BeamOnPlots.resize(NSamples, vector<TH1D*>(N1DPlots));

	vector< vector <TH1D*> > BeamOffPlots;
	BeamOffPlots.resize(NSamples, vector<TH1D*>(N1DPlots));

	vector< vector <TH1D*> > DirtPlots;
	DirtPlots.resize(NSamples, vector<TH1D*>(N1DPlots));

	vector< vector <TH2D*> > Covariances;
	Covariances.resize(NSamples, vector<TH2D*>(N1DPlots));

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

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------
	
		if (Runs[WhichRun] == "Run1") { DataPOT = tor860_wcut_Run1 ; }
		if (Runs[WhichRun] == "Run2") { DataPOT = tor860_wcut_Run2 ; }
		if (Runs[WhichRun] == "Run3") { DataPOT = tor860_wcut_Run3 ; }
		if (Runs[WhichRun] == "Run4") { DataPOT = tor860_wcut_Run4 ; }
		if (Runs[WhichRun] == "Run5") { DataPOT = tor860_wcut_Run5 ; }								
		
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface) * (SoftFidSurface / Nominal_UB_XY_Surface);

		// -------------------------------------------------------------------------------------

		TString FileName = MigrationMatrixPath+"WienerSVDCovarianceMatrices_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
		TFile* FileCovarianceMatrices = new TFile(FileName,"recreate");

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {
		
			TString ExactFileLocation = PathToFiles+CutExtension;

			MCFileSample[WhichSample][WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+NameOfSamples[WhichSample]+OverlaySample+"_"+Runs[WhichRun]+CutExtension+".root");

			BeamOnFileSample[WhichSample][WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOnSample+"_"+Runs[WhichRun]+CutExtension+".root");

			BeamOffFileSample[WhichSample][WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOffSample+"_"+Runs[WhichRun]+CutExtension+".root");

			DirtFileSample[WhichSample][WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+DirtSample+"_"+Runs[WhichRun]+CutExtension+".root");

		}

		// -------------------------------------------------------------------------------------

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				CC1pPlots[WhichSample][WhichPlot] = (TH1D*)(MCFileSample[WhichSample][WhichRun]->Get("CC1pReco"+PlotNames[WhichPlot]));
				NonCC1pPlots[WhichSample][WhichPlot] = (TH1D*)(MCFileSample[WhichSample][WhichRun]->Get("NonCC1pReco"+PlotNames[WhichPlot]));
				BeamOnPlots[WhichSample][WhichPlot] = (TH1D*)(BeamOnFileSample[WhichSample][WhichRun]->Get("Reco"+PlotNames[WhichPlot]));
				BeamOffPlots[WhichSample][WhichPlot] = (TH1D*)(BeamOffFileSample[WhichSample][WhichRun]->Get("Reco"+PlotNames[WhichPlot]));
				DirtPlots[WhichSample][WhichPlot] = (TH1D*)(DirtFileSample[WhichSample][WhichRun]->Get("Reco"+PlotNames[WhichPlot]));

				// -------------------------------------------------------------------------------------------------------

				TH1D* DataPlot = (TH1D*)BeamOnPlots[WhichSample][WhichPlot]->Clone();

				if (BeamOnSample == "BeamOn9") { // Beam On Sample, need to subtract ALL backgrounds

					DataPlot->Add(NonCC1pPlots[WhichSample][WhichPlot],-1.);
					DataPlot->Add(BeamOffPlots[WhichSample][WhichPlot],-1.);
					DataPlot->Add(DirtPlots[WhichSample][WhichPlot],-1.);

				} else { // Fake Data Studies, working with MC CC1p signal events as BeamOn

					DataPlot = CC1pPlots[WhichSample][WhichPlot];

				}

				// -------------------------------------------------------------------------------------------------------
				
				int NBinsX = CC1pPlots[WhichSample][WhichPlot]->GetXaxis()->GetNbins();
				int NBinsY = BeamOnPlots[WhichSample][WhichPlot]->GetXaxis()->GetNbins();	

				double XLow = CC1pPlots[WhichSample][WhichPlot]->GetXaxis()->GetBinLowEdge(1);
				double YLow = BeamOnPlots[WhichSample][WhichPlot]->GetXaxis()->GetBinLowEdge(1);

				double XHigh = CC1pPlots[WhichSample][WhichPlot]->GetXaxis()->GetBinLowEdge(NBinsX) + CC1pPlots[WhichSample][WhichPlot]->GetXaxis()->GetBinWidth(NBinsX);
				double YHigh = BeamOnPlots[WhichSample][WhichPlot]->GetXaxis()->GetBinLowEdge(NBinsY) + BeamOnPlots[WhichSample][WhichPlot]->GetXaxis()->GetBinWidth(NBinsY);

				TString XTitle = CC1pPlots[WhichSample][WhichPlot]->GetXaxis()->GetTitle();
				TString YTitle = BeamOnPlots[WhichSample][WhichPlot]->GetXaxis()->GetTitle();

				Covariances[WhichSample][WhichPlot] = new TH2D("Covariance_"+PlotNames[WhichPlot]+OverlaySample,";i bin "+XTitle+";j bin "+YTitle,NBinsX,XLow,XHigh,NBinsY,YLow,YHigh);

				// -------------------------------------------------------------------------------------------------------

				// Normalizing columns to Nj: Number of true events in the true bin j
				// Effectively creating a 2D efficiency map

				for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) { // MC bins

					// X Bin entry / error

					double DataEntryX = DataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
					double DataErrorX = DataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

					double MCEntryX = CC1pPlots[WhichSample][WhichPlot]->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
					double MCErrorX = CC1pPlots[WhichSample][WhichPlot]->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units; 

					for (int WhichYBin = 0; WhichYBin < NBinsY; WhichYBin++) { // Data bins

						// Y Bin entry / error

						double DataEntryY = DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
						double DataErrorY = DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

						double MCEntryY = CC1pPlots[WhichSample][WhichPlot]->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
						double MCErrorY = CC1pPlots[WhichSample][WhichPlot]->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units; 

						// CV

						double CovEntry = (DataEntryX - MCEntryX) * (DataEntryY - MCEntryY);
						Covariances[WhichSample][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovEntry);

						// Error

						double CovError = TMath::Sqrt( 
									TMath::Power(DataEntryY - MCEntryY,2.) * ( TMath::Power(DataErrorX,2.) + TMath::Power(MCErrorX,2.) ) +
									TMath::Power(DataEntryX - MCEntryX,2.) * ( TMath::Power(DataErrorY,2.) + TMath::Power(MCErrorY,2.) )
								);

						Covariances[WhichSample][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovError);

					}

				}
				
				FileCovarianceMatrices->cd();
				Covariances[WhichSample][WhichPlot]->Write();
				
				// ---------------------------------------------------------------------------------------				
	
				if (OverlaySample == "") {
		
					TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+NameOfSamples[WhichSample],
							    PlotNames[WhichPlot]+NameOfSamples[WhichSample],205,34,1024,768);
					PlotCanvas->cd();
					PlotCanvas->SetBottomMargin(0.16);
					PlotCanvas->SetLeftMargin(0.15);
					PlotCanvas->SetRightMargin(0.2);				
					
					gStyle->SetMarkerSize(1.5);
					gStyle->SetPaintTextFormat("4.3f");				
					
					Covariances[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
					Covariances[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
					Covariances[WhichSample][WhichPlot]->GetXaxis()->SetTitleSize(TextSize);
					Covariances[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(TextSize);				
					Covariances[WhichSample][WhichPlot]->GetXaxis()->CenterTitle();
					Covariances[WhichSample][WhichPlot]->GetXaxis()->SetNdivisions(5);
					
					Covariances[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
					Covariances[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
					Covariances[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(TextSize);
					Covariances[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(TextSize);				
					Covariances[WhichSample][WhichPlot]->GetYaxis()->CenterTitle();
					Covariances[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(5);
					Covariances[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(1.);				
									
					Covariances[WhichSample][WhichPlot]->GetZaxis()->SetLabelFont(FontStyle);
					Covariances[WhichSample][WhichPlot]->GetZaxis()->SetLabelSize(TextSize);
					Covariances[WhichSample][WhichPlot]->GetZaxis()->SetNdivisions(5);				

					Covariances[WhichSample][WhichPlot]->SetTitle(Runs[WhichRun]);	

//					Covariances[WhichSample][WhichPlot]->GetZaxis()->SetRangeUser(-0.015,0.015);
					Covariances[WhichSample][WhichPlot]->GetZaxis()->SetRangeUser(-0.05,0.05);
					Covariances[WhichSample][WhichPlot]->SetMarkerColor(kWhite);				
					Covariances[WhichSample][WhichPlot]->SetMarkerSize(1.5);
					Covariances[WhichSample][WhichPlot]->Draw("text colz e"); 

					TLatex* lat = new TLatex();
					lat->SetTextFont(FontStyle);
					lat->SetTextSize(TextSize);
					lat->DrawLatexNDC(0.8,0.93,"x10^{-76} cm^{4}");
					
					PlotCanvas->SaveAs(PlotPath+NameOfSamples[0]+"/WienerSVDCovarianceMatrices_"+PlotNames[WhichPlot]
						+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".pdf");
					
					delete PlotCanvas;				
				 
				}

			} // End of the loop over the plots
			
			MCFileSample[WhichSample][WhichRun]->Close();
			BeamOnFileSample[WhichSample][WhichRun]->Close();
			BeamOffFileSample[WhichSample][WhichRun]->Close();
			DirtFileSample[WhichSample][WhichRun]->Close();

		} // End of the loop over the samples

		FileCovarianceMatrices->Close();

		std::cout << "File with Covariance matrices " << FileName << " has been created" << std::endl;

	} // End of the loop over the runs	

} // End of the program
