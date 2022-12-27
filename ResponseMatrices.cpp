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

void ResponseMatrices(TString OverlaySample, bool DetVar = false, TString Tune = "", bool NuWro = false) {

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

	const int N2DPlots = PlotNames.size();
	cout << "Number of 2D Plots = " << N2DPlots << endl;

	// -------------------------------------------------------------------------------------------------------------------------------------
	
	vector<TString> NameOfSamples;
	NameOfSamples.push_back("Overlay9");
	if (NuWro) { NameOfSamples[0] = "Overlay9NuWro"; }
	const int NSamples = NameOfSamples.size();
		
	// -------------------------------------------------------------------------------------------------------------------------------------

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

	const int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;	

	// -------------------------------------------------------------------------------------

	vector< vector<TFile*>> FileSample;
	FileSample.resize(NSamples, vector<TFile*>(NRuns));

	vector< vector<TFile*>> TrueFileSample;
	TrueFileSample.resize(NSamples, vector<TFile*>(NRuns));

	vector< vector <TH2D*> > Plots;
	Plots.resize(NSamples, vector<TH2D*>(N2DPlots));

	vector< vector <TH1D*> > TruePlots;
	TruePlots.resize(NSamples, vector<TH1D*>(N2DPlots));

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// -------------------------------------------------------------------------------------

		TString FileName = MigrationMatrixPath+Tune+"FileResponseMatrices_"+\
				NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
		TFile* FileResponseMatrices = new TFile(FileName,"recreate");

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {
		
			TString ExactFileLocation = PathToFiles+CutExtension;

			FileSample[WhichSample][WhichRun] = TFile::Open(ExactFileLocation+"/"+Tune+"STVStudies_"+NameOfSamples[WhichSample]+"_"+\
							  Runs[WhichRun]+OverlaySample+CutExtension+".root","readonly");

			TrueFileSample[WhichSample][WhichRun] = TFile::Open(PathToFiles+"/"+Tune+"TruthSTVAnalysis_"+NameOfSamples[WhichSample]+"_"+\
							  Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root","readonly");

			// Jul 8 2021: after discussion with Xin, if flux variations, the truth level should always be the CV

			if (string(OverlaySample).find("fluxes") != std::string::npos) {

				TrueFileSample[WhichSample][WhichRun] = TFile::Open(PathToFiles+"/"+Tune+"TruthSTVAnalysis_"+NameOfSamples[WhichSample]+"_"+\
							  Runs[WhichRun]+"_"+UBCodeVersion+".root","readonly");

			}

		}

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			for (int WhichPlot = 0; WhichPlot < N2DPlots; WhichPlot ++) {

				Plots[WhichSample][WhichPlot] = (TH2D*)(FileSample[WhichSample][WhichRun]->Get("POTScaledCC1pReco"+PlotNames[WhichPlot]+"2D"));
				TruePlots[WhichSample][WhichPlot] = (TH1D*)(TrueFileSample[WhichSample][WhichRun]->Get("True"+PlotNames[WhichPlot]));
				
				int NBinsX = Plots[WhichSample][WhichPlot]->GetXaxis()->GetNbins();
				int NBinsY = Plots[WhichSample][WhichPlot]->GetYaxis()->GetNbins();	

				// -------------------------------------------------------------------------------------------------------

				// Normalizing columns to Nj: Number of true events in that bin
				// Effectively creating 2D efficiancy map

				for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) { // True bins

					// True bin entry / error

					double TrueEntry = TruePlots[WhichSample][WhichPlot]->GetBinContent(WhichXBin+1);
					double TrueError = TruePlots[WhichSample][WhichPlot]->GetBinError(WhichXBin+1); 

					for (int WhichYBin = 0; WhichYBin < NBinsY; WhichYBin++) { // Reco bins

						// Local 2D Bin

						double Entry = Plots[WhichSample][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1);
						double Error = Plots[WhichSample][WhichPlot]->GetBinError(WhichXBin+1,WhichYBin+1); 

						// CV

						double Eff = Entry / TrueEntry;
						Plots[WhichSample][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,Eff);

						// Error

						double EffError = Eff * TMath::Sqrt( TMath::Power(Error/Entry,2.) + TMath::Power(TrueError/TrueEntry,2.) );
						Plots[WhichSample][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,EffError);

					}

				}
				
				FileResponseMatrices->cd();
				Plots[WhichSample][WhichPlot]->Write();
				
				// ---------------------------------------------------------------------------------------				
	
				if (OverlaySample == "") {
		
					TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+NameOfSamples[WhichSample],
							    PlotNames[WhichPlot]+NameOfSamples[WhichSample],205,34,1024,768);
					PlotCanvas->cd();
					PlotCanvas->SetBottomMargin(0.17);
					PlotCanvas->SetLeftMargin(0.15);
					PlotCanvas->SetRightMargin(0.15);				
					
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
					Plots[WhichSample][WhichPlot]->GetZaxis()->SetNdivisions(5);				

					Plots[WhichSample][WhichPlot]->SetTitle(Runs[WhichRun] + " Response Matrix, " + LatexLabel[PlotNames[WhichPlot]]);	

					Plots[WhichSample][WhichPlot]->GetZaxis()->SetRangeUser(0,1.1*Plots[WhichSample][WhichPlot]->GetMaximum());
					Plots[WhichSample][WhichPlot]->SetMarkerColor(kWhite);				
					Plots[WhichSample][WhichPlot]->SetMarkerSize(0.9);

					//------------------------------//

					// The N-dimensional analysis has been developed based on the bin number, not the actual range

					if (string(PlotNames[WhichPlot]).find("Serial") != std::string::npos) {	

						TString XaxisTitle = Plots[WhichSample][WhichPlot]->GetXaxis()->GetTitle();
						XaxisTitle.ReplaceAll("deg","bin #");
						XaxisTitle.ReplaceAll("GeV/c","bin #");
						XaxisTitle.ReplaceAll("GeV","bin #");				
						Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitle(XaxisTitle);

						TString YaxisTitle = Plots[WhichSample][WhichPlot]->GetYaxis()->GetTitle();
						YaxisTitle.ReplaceAll("deg","bin #");
						YaxisTitle.ReplaceAll("GeV/c","bin #");
						YaxisTitle.ReplaceAll("GeV","bin #");				
						Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitle(YaxisTitle);						

						Plots[WhichSample][WhichPlot]->Draw("colz");						

					} else {

						Plots[WhichSample][WhichPlot]->Draw("text colz e");

					} 
					
					PlotCanvas->SaveAs(PlotPath+NameOfSamples[0]+"/"+Tune+"ResponseMatrices_"+PlotNames[WhichPlot]
						+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".pdf");
					
					delete PlotCanvas;				
				 
				}

			} // End of the loop over the plots
			
			FileSample[WhichSample][WhichRun]->Close();
			TrueFileSample[WhichSample][WhichRun]->Close();

		} // End of the loop over the samples

		FileResponseMatrices->Close();

		std::cout << "File with Response matrices " << FileName << " has been created" << std::endl;

	} // End of the loop over the runs	

} // End of the program
