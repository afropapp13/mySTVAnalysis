#include <TFile.h>
#include <TF1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>

#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

void MigrationMatrices(TString OverlaySample) {

	// -------------------------------------------------------------------------------------

//	TString PathToFiles = "/uboone/data/users/"+UserID+"/myEvents/OutputFiles/"+UBCodeVersion+"/";
//	TString PlotPath = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myPlots/"+UBCodeVersion+"/"; 
//	TString MigrationMatrixPath = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myMigrationMatrices/"+UBCodeVersion+"/"; 

	// -------------------------------------------------------------------------------------

	TH2D::SetDefaultSumw2();
	
	double TextSize = 0.07;
	const Int_t NCont = 999;	
	gStyle->SetPalette(55); 
	gStyle->SetNumberContours(NCont); 
	gStyle->SetTitleSize(TextSize,"t"); 
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

//	vector<TString> PlotNames;
//	PlotNames.push_back("DeltaPTPlot"); 
//	PlotNames.push_back("DeltaAlphaTPlot"); 
//	PlotNames.push_back("DeltaPhiTPlot"); 
//	PlotNames.push_back("MuonMomentumPlot"); 
//	PlotNames.push_back("MuonPhiPlot"); 
//	PlotNames.push_back("MuonCosThetaPlot");
//	PlotNames.push_back("ProtonMomentumPlot"); 
//	PlotNames.push_back("ProtonPhiPlot"); 
//	PlotNames.push_back("ProtonCosThetaPlot");
//	PlotNames.push_back("ECalPlot"); 
//	PlotNames.push_back("EQEPlot"); 
//	PlotNames.push_back("Q2Plot");

//	PlotNames.push_back("CCQEMuonMomentumPlot"); 
//	PlotNames.push_back("CCQEMuonCosThetaPlot"); 
//	PlotNames.push_back("CCQEProtonMomentumPlot"); 
//	PlotNames.push_back("CCQEProtonCosThetaPlot");

	const int N2DPlots = PlotNames.size();
	cout << "Number of 2D Plots = " << N2DPlots << endl;

	// -------------------------------------------------------------------------------------------------------------------------------------
	
	vector<TString> NameOfSamples;
	NameOfSamples.push_back("Overlay9");
	const int NSamples = NameOfSamples.size();
		
	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");				
	Runs.push_back("Combined");				

	const int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;	

	// -------------------------------------------------------------------------------------

	vector< vector<TFile*>> FileSample;
	FileSample.resize(NSamples, vector<TFile*>(NRuns));
	vector< vector <TH2D*> > Plots;
	Plots.resize(NSamples, vector<TH2D*>(N2DPlots));

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// -------------------------------------------------------------------------------------

		TString FileName = MigrationMatrixPath+"FileMigrationMatrices_"+\
				NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
		TFile* FileMigrationMatrices = new TFile(FileName,"recreate");

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {
		
			TString ExactFileLocation = PathToFiles+CutExtension;

			FileSample[WhichSample][WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+NameOfSamples[WhichSample]+"_"+\
							  Runs[WhichRun]+CutExtension+".root");

		}

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			for (int WhichPlot = 0; WhichPlot < N2DPlots; WhichPlot ++) {


				Plots[WhichSample][WhichPlot] = (TH2D*)(FileSample[WhichSample][WhichRun]->Get("CC1pReco"+PlotNames[WhichPlot]+"2D"));
				
				int NBinsX = Plots[WhichSample][WhichPlot]->GetXaxis()->GetNbins();
				int NBinsY = Plots[WhichSample][WhichPlot]->GetYaxis()->GetNbins();								

				// Normalizing columns to 1

				for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

					int NEventsInColumn = 0;

					for (int WhichYBin = 0; WhichYBin < NBinsY; WhichYBin++) {

						NEventsInColumn += Plots[WhichSample][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1);

					}

					for (int WhichYBin = 0; WhichYBin < NBinsY; WhichYBin++) {
	
						if (NEventsInColumn > 0) {

							double FracErr = Plots[WhichSample][WhichPlot]->GetBinError(WhichXBin+1,WhichYBin+1) /\
									 Plots[WhichSample][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1);

							// CV

							double CV = double(Plots[WhichSample][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1))/\
								    double(NEventsInColumn);

							Plots[WhichSample][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CV);

							// Error

							double error = CV * TMath::Sqrt( TMath::Power(FracErr,2.) + 1./double(NEventsInColumn) );
							
							Plots[WhichSample][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,error) ; 

						} else { 

							Plots[WhichSample][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,0.); 
							Plots[WhichSample][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,0.); 
						}

					}

				}
				
				FileMigrationMatrices->cd();
				Plots[WhichSample][WhichPlot]->Write();
				
				// ---------------------------------------------------------------------------------------				
	
				if (OverlaySample == "") {
		
					TString PlotCanvasName = Runs[WhichRun]+"_"+PlotNames[WhichPlot]+NameOfSamples[WhichSample];
					TCanvas* PlotCanvas = new TCanvas(PlotCanvasName,PlotCanvasName,205,34,1024,768);
					PlotCanvas->cd();
					PlotCanvas->SetBottomMargin(0.16);
					PlotCanvas->SetLeftMargin(0.15);
					PlotCanvas->SetRightMargin(0.12);				
					
					gStyle->SetMarkerSize(1.5);
					gStyle->SetPaintTextFormat("4.2f");				
					
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleSize(TextSize);
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(TextSize);				
					Plots[WhichSample][WhichPlot]->GetXaxis()->CenterTitle();
					Plots[WhichSample][WhichPlot]->GetXaxis()->SetNdivisions(5);
					
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(TextSize);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(TextSize);				
					Plots[WhichSample][WhichPlot]->GetYaxis()->CenterTitle();
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(5);
					Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(1.);				
									
					Plots[WhichSample][WhichPlot]->GetZaxis()->SetLabelFont(FontStyle);
					Plots[WhichSample][WhichPlot]->GetZaxis()->SetLabelSize(TextSize);
					Plots[WhichSample][WhichPlot]->GetZaxis()->SetNdivisions(8);				
					Plots[WhichSample][WhichPlot]->GetZaxis()->SetRangeUser(-0.1,1.);

					Plots[WhichSample][WhichPlot]->SetMarkerColor(kWhite);				
					Plots[WhichSample][WhichPlot]->SetMarkerSize(0.9);
					Plots[WhichSample][WhichPlot]->SetTitle(Runs[WhichRun]);					
					Plots[WhichSample][WhichPlot]->Draw("text colz e"); 
					
					PlotCanvas->SaveAs(PlotPath+NameOfSamples[0]+"/MigrationMatrices_"+PlotNames[WhichPlot]
						+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".pdf");
					
					delete PlotCanvas;				
				 
				}

			} // End of the loop over the plots
			
			FileSample[WhichSample][WhichRun]->Close();

		} // End of the loop over the samples

		FileMigrationMatrices->Close();

		std::cout << "File with migration matrices " << FileName << " has been created" << std::endl;

	} // End of the loop over the runs	

} // End of the program
