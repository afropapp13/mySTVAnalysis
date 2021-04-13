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

void WienerSVD_POT_CovarianceMatrices(TString OverlaySample,int Universe = -1,TString BeamOnSample = "BeamOn9",TString BeamOffSample = "ExtBNB9",TString DirtSample = "OverlayDirt9") {

	// -------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
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

		TString FileName = MigrationMatrixPath+"WienerSVD_POT_CovarianceMatrices_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
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
				
				int NBins = BeamOnPlots[WhichSample][WhichPlot]->GetXaxis()->GetNbins();
				const double* ArrayBins = BeamOnPlots[WhichSample][WhichPlot]->GetXaxis()->GetXbins()->GetArray();
				TString XTitle = BeamOnPlots[WhichSample][WhichPlot]->GetXaxis()->GetTitle();

				Covariances[WhichSample][WhichPlot] = new TH2D("Covariance_"+PlotNames[WhichPlot]+OverlaySample,";i bin "+XTitle+";j bin "+XTitle,NBins,ArrayBins,NBins,ArrayBins);

				// -------------------------------------------------------------------------------------------------------

				for (int WhichXBin = 0; WhichXBin < NBins; WhichXBin++) { // MC bins

					// X Bin entry / error

					double DataEntryX = DataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
					double DataErrorX = DataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

					double POTDataEntryX = (1+POTUncertainty) * DataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
					double POTDataErrorX = (1+POTUncertainty) * DataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

					for (int WhichYBin = 0; WhichYBin < NBins; WhichYBin++) { // Data bins

						// Y Bin entry / error

						double DataEntryY = DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
						double DataErrorY = DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

						double POTDataEntryY = (1+POTUncertainty) * DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
						double POTDataErrorY = (1+POTUncertainty) * DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

						// CV

						double CovEntry = TMath::Max((POTDataEntryX - DataEntryX) * (POTDataEntryY - DataEntryY),1E-8);
						Covariances[WhichSample][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovEntry);

						// Error

//						double CovError = TMath::Sqrt( 
//									TMath::Power(DataEntryY - POTDataEntryY,2.) * ( TMath::Power(DataErrorX,2.) + TMath::Power(POTDataErrorX,2.) ) +
//									TMath::Power(DataEntryX - POTDataEntryX,2.) * ( TMath::Power(DataErrorY,2.) + TMath::Power(POTDataErrorY,2.) )
//								);

						double CovError = 1E-10;

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
					Covariances[WhichSample][WhichPlot]->GetZaxis()->SetRangeUser(0.,0.000199);
					Covariances[WhichSample][WhichPlot]->SetMarkerColor(kWhite);				
					Covariances[WhichSample][WhichPlot]->SetMarkerSize(1.5);
//					Covariances[WhichSample][WhichPlot]->Draw("text colz e"); 
					Covariances[WhichSample][WhichPlot]->Draw("colz"); 

					TLatex* lat = new TLatex();
					lat->SetTextFont(FontStyle);
					lat->SetTextSize(TextSize);
					lat->DrawLatexNDC(0.6,0.94,"x10^{-76} cm^{4}");
					
					PlotCanvas->SaveAs(PlotPath+NameOfSamples[0]+"/WienerSVD_POT_CovarianceMatrices_"+PlotNames[WhichPlot]
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

		std::cout << "File with POT Covariance matrices " << FileName << " has been created" << std::endl;

	} // End of the loop over the runs	

} // End of the program
