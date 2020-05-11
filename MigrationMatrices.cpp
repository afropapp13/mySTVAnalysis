#include <TFile.h>
#include <TF1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>

#include <iostream>
#include <vector>

#include  "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/ToString.cpp"
#include "./Constants.h"

using namespace std;
using namespace Constants;

void MigrationMatrices(TString OverlaySample) {

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

	TString PathToFiles = "../myEvents/OutputFiles";

	// -------------------------------------------------------------------------------------


	TH2D::SetDefaultSumw2();

	vector<TString> PlotNames;
	PlotNames.push_back("DeltaPTPlot"); PlotNames.push_back("DeltaAlphaTPlot"); PlotNames.push_back("DeltaPhiTPlot"); 
	PlotNames.push_back("MuonMomentumPlot"); PlotNames.push_back("MuonPhiPlot"); PlotNames.push_back("MuonCosThetaPlot");
	PlotNames.push_back("ProtonMomentumPlot"); PlotNames.push_back("ProtonPhiPlot"); PlotNames.push_back("ProtonCosThetaPlot");
	PlotNames.push_back("ECalPlot"); PlotNames.push_back("EQEPlot"); PlotNames.push_back("Q2Plot");

	const int N2DPlots = PlotNames.size();
	cout << "Number of 2D Plots = " << N2DPlots << endl;

	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); gStyle->SetTitleFont(FontStyle,"t"); SetOffsetAndSize();

	vector<TString> NameOfSamples;

	NameOfSamples.push_back("Overlay9");

	// -------------------------------------------------------------------------------------

	const int NSamples = NameOfSamples.size();
	TCanvas* PlotCanvas[NSamples][N2DPlots] = {}; TFile* FileSample[NSamples] = {};
	TH2D* Plots[NSamples][N2DPlots] = {};

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// -------------------------------------------------------------------------------------

		TString FileName = "myMigrationMatrices/"+UBCodeVersion+"/FileMigrationMatrices_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
		TFile* FileMigrationMatrices = new TFile(FileName,"recreate");

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			FileSample[WhichSample] = TFile::Open(PathToFiles+"/"+UBCodeVersion+"/"+CutExtension+"/STVStudies_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+CutExtension+".root");

		}

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			for (int WhichPlot = 0; WhichPlot < N2DPlots; WhichPlot ++) {

				PlotCanvas[WhichSample][WhichPlot] = new TCanvas(PlotNames[WhichPlot]+NameOfSamples[WhichSample],
						    PlotNames[WhichPlot]+NameOfSamples[WhichSample],205,34,1024,768);
				PlotCanvas[WhichSample][WhichPlot]->cd();
				gStyle->SetMarkerSize(1.5);
				gStyle->SetPaintTextFormat("4.2f");
				Plots[WhichSample][WhichPlot] = (TH2D*)(FileSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]+"2D"));
				Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
				Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
				Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
				Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
				Plots[WhichSample][WhichPlot]->GetZaxis()->SetLabelFont(FontStyle);
				Plots[WhichSample][WhichPlot]->GetZaxis()->SetLabelSize(0.03);

				CenterAxisTitle(Plots[WhichSample][WhichPlot]);

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

							// CV
							double CV = double(Plots[WhichSample][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1))/double(NEventsInColumn);
							Plots[WhichSample][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CV);

							// Error
							double error = sqrt(
							TMath::Power(Plots[WhichSample][WhichPlot]->GetBinError(WhichXBin+1,WhichYBin+1)/double(NEventsInColumn),2.) +
							TMath::Power(Plots[WhichSample][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1) * sqrt(NEventsInColumn)/double(NEventsInColumn*NEventsInColumn),2.)
							);
							Plots[WhichSample][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,error) ; 

						} else { 

							Plots[WhichSample][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,0.); 
							Plots[WhichSample][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,0.); 
						}

					}

				}

				// ------------------------------------------------------------------------------------------------------------------------------

				// Normalizing rows to 1

	//			for (int WhichYBin = 0; WhichYBin < NBinsY; WhichYBin++) {

	//				int NEventsInRow = 0;

	//				for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

	//					NEventsInRow += Plots[WhichSample][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1);

	//				}

	//				for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {
	//	
	//					if (NEventsInRow > 0) {
	//						// CV
	//						double CV = Plots[WhichSample][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1)/double(NEventsInRow);
	//						Plots[WhichSample][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CV);
	//	
	//						// Error
	//						double error = sqrt(
	//						TMath::Power(Plots[WhichSample][WhichPlot]->GetBinError(WhichXBin+1,WhichYBin+1)/double(NEventsInRow),2.) +
	//						TMath::Power(Plots[WhichSample][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1) * sqrt(NEventsInRow)/double(NEventsInRow*NEventsInRow),2.)
	//						);
	//						Plots[WhichSample][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,error) ; 

	//					} else { 

	//						Plots[WhichSample][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,0.); 
	//						Plots[WhichSample][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,0.); 
	//					}

	//				}

	//			}

				Plots[WhichSample][WhichPlot]->GetZaxis()->SetRangeUser(0,1.);
				FileMigrationMatrices->cd();
				Plots[WhichSample][WhichPlot]->Write();
				Plots[WhichSample][WhichPlot]->Draw("text colz");

				if (OverlaySample == "") {
					PlotCanvas[WhichSample][WhichPlot]->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/"+NameOfSamples[0]+"/MigrationMatrices_"+PlotNames[WhichPlot]
					+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".pdf");
				}

				delete PlotCanvas[WhichSample][WhichPlot];

			} // End of the loop over the plots

		} // End of the loop over the samples

	//	FileMigrationMatrices->cd();
	//	FileMigrationMatrices->Write();

		FileMigrationMatrices->Close();

		std::cout << "File with migration matrices " << FileName << " has been created" << std::endl;

		//delete []PlotNames;

	} // End of the loop over the runs	

} // End of the program
