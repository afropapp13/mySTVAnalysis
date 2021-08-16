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

void PrettyPlot(TH2D* h, TString Label, TString Runs, TString PlotNames, TString NameOfSamples, TString OverlaySample, bool Store = false) {

	TString PlotCanvasName = Label + "_" + Runs+"_"+PlotNames+NameOfSamples;
	TCanvas* PlotCanvas = new TCanvas(PlotCanvasName,PlotCanvasName,205,34,1024,768);
	PlotCanvas->cd();
	PlotCanvas->SetBottomMargin(0.16);
	PlotCanvas->SetLeftMargin(0.15);
	PlotCanvas->SetRightMargin(0.15);				
	
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
					
	h->GetZaxis()->SetLabelFont(FontStyle);
	h->GetZaxis()->SetLabelSize(TextSize);
	h->GetZaxis()->SetNdivisions(5);				

	h->SetTitle(Runs);	

	//h->GetZaxis()->SetRangeUser(0,1.);
	h->SetMarkerColor(kWhite);				
	h->SetMarkerSize(0.8);
	h->Draw("text colz e"); 
	
	if (Store) { PlotCanvas->SaveAs(PlotPath+NameOfSamples+"/"+Label+"BinInOutMigrationMatrices_"+PlotNames+NameOfSamples+"_"+Runs+OverlaySample+"_"+UBCodeVersion+".pdf"); }
	
	//delete PlotCanvas;

}

// ---------------------------------------------------------------------------------------------------------

void BinInBinOutMigration(TString OverlaySample) {

	// -------------------------------------------------------------------------------------

//	TString PathToFiles = "/uboone/data/users/"+UserID+"/myEvents/OutputFiles/"+UBCodeVersion+"/";
//	TString PlotPath = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myPlots/"+UBCodeVersion+"/"; 
//	TString MigrationMatrixPath = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myMigrationMatrices/"+UBCodeVersion+"/"; 

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

//		TString FileName = MigrationMatrixPath+"FileMigrationMatrices_"+\
//				NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
//		TFile* FileMigrationMatrices = new TFile(FileName,"recreate");

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {
		
			TString ExactFileLocation = PathToFiles+CutExtension;

			FileSample[WhichSample][WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+NameOfSamples[WhichSample]+"_"+\
							  Runs[WhichRun]+OverlaySample+CutExtension+".root");

		}

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			for (int WhichPlot = 0; WhichPlot < N2DPlots; WhichPlot ++) {

				// -----------------------------------------------------------------------------------------------------------------------

				// Grab the base 2D plot with the number of events

				Plots[WhichSample][WhichPlot] = (TH2D*)(FileSample[WhichSample][WhichRun]->Get("CC1pReco"+PlotNames[WhichPlot]+"2D"));

				int NBinsX = Plots[WhichSample][WhichPlot]->GetXaxis()->GetNbins();
				int NBinsY = Plots[WhichSample][WhichPlot]->GetYaxis()->GetNbins();

				// -----------------------------------------------------------------------------------------------------------------------

				// Clone the base plot to modify it

				TH2D* ClonedPlots = (TH2D*)(Plots[WhichSample][WhichPlot]->Clone());

				double NonDiagEventsInRow[NBinsX];
				memset( NonDiagEventsInRow, 0, NBinsX*sizeof(double) );

				double NonDiagEventsInColumn[NBinsY];
				memset( NonDiagEventsInColumn, 0, NBinsY*sizeof(double) );

				// -----------------------------------------------------------------------------------------------------------------------

				// Zero out everything that is not a diagonal element in the clone
				// Calculate the sum of the non-diagonal elements in a row

				for (int WhichYBin = 1; WhichYBin < NBinsY+1; WhichYBin++) {

					NonDiagEventsInRow[WhichYBin-1] = 0;

					for (int WhichXBin = 1; WhichXBin < NBinsX+1; WhichXBin++) {

						if (WhichXBin != WhichYBin) {

							NonDiagEventsInRow[WhichYBin-1] += Plots[WhichSample][WhichPlot]->GetBinContent(WhichXBin,WhichYBin);
							ClonedPlots->SetBinContent(WhichXBin,WhichYBin,0.);
							ClonedPlots->SetBinError(WhichXBin,WhichYBin,0.);

						}

					}

				}

				// -----------------------------------------------------------------------------------------------------------------------
											
				// Calculate the sum of the non-diagonal elements in a column

				for (int WhichXBin = 1; WhichXBin < NBinsX+1; WhichXBin++) {

					NonDiagEventsInColumn[WhichXBin-1] = 0;

					for (int WhichYBin = 1; WhichYBin < NBinsY+1; WhichYBin++) {

						if (WhichXBin != WhichYBin) {

							NonDiagEventsInColumn[WhichXBin-1] += Plots[WhichSample][WhichPlot]->GetBinContent(WhichXBin,WhichYBin);

						}

					}

				}

				// -----------------------------------------------------------------------------------------------------------------------

				// Grab the clone with the diagonal elements
				// Add the events that migrate in, subtract the events that migrate out

	 
				for (int WhichXBin = 1; WhichXBin < NBinsX+1; WhichXBin++) {

					double CurrentContent = ClonedPlots->GetBinContent(WhichXBin,WhichXBin);
					double NewContent = CurrentContent + NonDiagEventsInRow[WhichXBin-1] - NonDiagEventsInColumn[WhichXBin-1];
					double BinInOut = NonDiagEventsInRow[WhichXBin-1] - NonDiagEventsInColumn[WhichXBin-1];

//					ClonedPlots->SetBinContent(WhichXBin,WhichXBin,NewContent);
					ClonedPlots->SetBinContent(WhichXBin,WhichXBin,BinInOut/CurrentContent*100.);

				}

				// -----------------------------------------------------------------------------------------------------------------------
				
//				FileMigrationMatrices->cd();
//				Plots[WhichSample][WhichPlot]->Write();
				
				// ---------------------------------------------------------------------------------------				
	
				if (OverlaySample == "") {

					// ---------------------------------------------------------------------------------------				

					PrettyPlot(Plots[WhichSample][WhichPlot],"Default",Runs[WhichRun],PlotNames[WhichPlot],NameOfSamples[WhichSample],OverlaySample,true);
		
					// ---------------------------------------------------------------------------------------

					ClonedPlots->GetZaxis()->SetTitle("Bin InOut / Diag [%]");
					ClonedPlots->GetZaxis()->SetTitleOffset(1.3);
					PrettyPlot(ClonedPlots,"BinInOut",Runs[WhichRun],PlotNames[WhichPlot],NameOfSamples[WhichSample],OverlaySample,true);

					// ---------------------------------------------------------------------------------------								
				 
				}

			} // End of the loop over the plots
			
			//FileSample[WhichSample][WhichRun]->Close();

		} // End of the loop over the samples

//		FileMigrationMatrices->Close();

//		std::cout << "File with migration matrices " << FileName << " has been created" << std::endl;

	} // End of the loop over the runs	

} // End of the program
