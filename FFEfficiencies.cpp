#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TLatex.h>
#include <TGaxis.h>

#include <iostream>
#include <vector>

#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

void FFEfficiencies(TString OverlaySample) {

	// -------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;
	gStyle->SetOptStat(0);	
	TGaxis::SetMaxDigits(3);
	TGaxis::SetExponentOffset(-0.05, 0., "y");	
	
	double TextSize = 0.07;

	// -------------------------------------------------------------------------------------

//	TString PathToFiles = "/uboone/data/users/"+UserID+"/myEvents/OutputFiles/"+UBCodeVersion;
//	TString TrueSTVPath = PathToFiles+"/";
//	TString EfficiencyPath = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myEfficiencies/"+UBCodeVersion+"/"; 
//	TString PlotPath = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myPlots/"+UBCodeVersion+"/"; 

	// -------------------------------------------------------------------------------------

	int NEventsPassingSelectionCuts = 0;
	TString CutExtension = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
	VectorCuts.push_back("_NuScore");
	VectorCuts.push_back("_ThreePlaneLogChi2");
	VectorCuts.push_back("_Collinearity");

	int NCuts = (int)(VectorCuts.size());	
	for (int i = 0; i < NCuts; i++) { CutExtension = CutExtension + VectorCuts[i]; }

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
//	PlotNames.push_back("kMissPlot");
//	PlotNames.push_back("PMissPlot");
//	PlotNames.push_back("PMissMinusPlot");

	// -------------------------------------------------------------------------------------

	const int N1DPlots = PlotNames.size();
	//cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ------------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");
	Runs.push_back("Run2");
	Runs.push_back("Run3");
	Runs.push_back("Run4");
	Runs.push_back("Run5");				

	int NRuns = (int)(Runs.size());
	//cout << "Number of Runs = " << NRuns << endl;

	// -------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------------------------------------

		// To be removed when all the runs / systematics will be ready to process

		if (Runs[WhichRun] == "Run2") { continue; }
		if (Runs[WhichRun] == "Run4") { continue; }
		if (Runs[WhichRun] == "Run5") { continue; }

		if (Runs[WhichRun] == "Run1" && 
			(OverlaySample == "_LYAttenuation" || OverlaySample == "_WireModX" || OverlaySample == "_WireModYZ" || 
			OverlaySample == "_WireModThetaYZ" || OverlaySample == "_WireModThetaXZ" || OverlaySample == "_dEdx" || 
			OverlaySample == "_Recombination2" || OverlaySample == "_SCE" ) ) { continue; }

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------------------------------------------

		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
		vector<vector<TH1D*> > PlotsTrueReco; PlotsTrueReco.clear();

		gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(TextSize,"t");

		vector<TString> LabelsOfSamples;
		vector<TString> NameOfSamples;

		NameOfSamples.push_back("Overlay9");

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();
		vector<TFile*> TruthFileSample; TruthFileSample.clear();
		vector<TFile*> MigrationMatricesFile; MigrationMatricesFile.clear();

		TString Name = "";

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			TString STVPath = PathToFiles+"/"+CutExtension+"/";
			TString STVName = "STVStudies_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+CutExtension+".root";
			FileSample.push_back(TFile::Open(STVPath+STVName));
			
			TString TrueSTVName = "TruthSTVAnalysis_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
			TruthFileSample.push_back(TFile::Open(TrueSTVPath+TrueSTVName));

			TString MigrationMatrixName = "FileMigrationMatrices_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
			MigrationMatricesFile.push_back(TFile::Open(MigrationMatrixPath+MigrationMatrixName));

			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();
			vector<TH1D*> CurrentPlotsTrueReco; CurrentPlotsTrueReco.clear();
			vector<TH2D*> CurrentMigrationMatrix; CurrentMigrationMatrix.clear();

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

				TH1D* histTrue = (TH1D*)(TruthFileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
				CurrentPlotsTrue.push_back(histTrue);

				TH1D* histTrueReco = (TH1D*)(FileSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]));
				CurrentPlotsTrueReco.push_back(histTrueReco);

				TH2D* histMigrationMatrix = (TH2D*)(FileSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]+"2D"));
				CurrentMigrationMatrix.push_back(histMigrationMatrix);
		
			}

			PlotsTrue.push_back(CurrentPlotsTrue);
			PlotsTrueReco.push_back(CurrentPlotsTrueReco);

		}

		// Loop over the samples

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			TString EfficiencyName = "ForwardFoldEfficiences_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root"; 
			Name = FileEfficienciesPath + EfficiencyName;
			TFile* FileEfficiences = new TFile(Name,"recreate");

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				// ---------------------------------------------------------------------------------------------------------------------------	

				// Ratios to extract the forward folded efficiencies			

				// this is the line that has to be replaced by the relevant new function
				TH1D* pEffPlot = ForwardFold(PlotsTrue[WhichSample][WhichPlot],PlotsTrue[WhichSample][WhichPlot]);

				FileEfficiences->cd();
				pEffPlot->Write();

				// Plotting only for the nominal overlay sample
				if (WhichSample == 0 && OverlaySample == "") { 

					pEffPlot->SetLineWidth(3);
					pEffPlot->SetLineColor(kBlack);
					pEffPlot->SetMarkerStyle(20);

					pEffPlot->GetXaxis()->CenterTitle();
					pEffPlot->GetXaxis()->SetTitleFont(FontStyle);
					pEffPlot->GetXaxis()->SetLabelFont(FontStyle);
					pEffPlot->GetXaxis()->SetTitle(PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->GetTitle());
					pEffPlot->GetXaxis()->SetTitleSize(TextSize);
					pEffPlot->GetXaxis()->SetLabelSize(TextSize);
					pEffPlot->GetXaxis()->SetNdivisions(5);

					pEffPlot->GetYaxis()->CenterTitle();
					pEffPlot->GetYaxis()->SetTitleFont(FontStyle);
					pEffPlot->GetYaxis()->SetTitleSize(TextSize);
					pEffPlot->GetYaxis()->SetLabelFont(FontStyle);
					pEffPlot->GetYaxis()->SetNdivisions(6);
					pEffPlot->GetYaxis()->SetTitleSize(TextSize);
					pEffPlot->GetYaxis()->SetLabelSize(TextSize);
					pEffPlot->GetYaxis()->SetTitle("Efficiency");
					pEffPlot->GetYaxis()->SetRangeUser(0.,1.2*pEffPlot->GetMaximum());

					TString CanvasEffName = NameOfSamples[WhichSample]+"_"+"Eff"+PlotNames[WhichPlot];
					TCanvas* PlotEffCanvas = new TCanvas(CanvasEffName,CanvasEffName,205,34,1024,768);
					PlotEffCanvas->cd();
					PlotEffCanvas->SetBottomMargin(0.16);
					PlotEffCanvas->SetLeftMargin(0.18);

					PlotEffCanvas->cd();
					pEffPlot->Draw();

					TLatex *textEff = new TLatex();
					textEff->SetTextFont(FontStyle);
					textEff->SetTextSize(TextSize);
					textEff->DrawTextNDC(0.2, 0.8, Runs[WhichRun]);
				
					TString CanvasEffPath = PlotPath+NameOfSamples[WhichSample]+"/";
					TString CanvasEffRatioName = "Eff"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".pdf";
					PlotEffCanvas->SaveAs(CanvasEffPath + CanvasEffRatioName);

					delete PlotEffCanvas;
									
				}				

			} // End of the loop over the plots

			FileEfficiences->Close();

			std::cout << std::endl << "-----------------------------------------------------------------------" << std::endl << std::endl;
			std::cout << std::endl << "Efficiency file " << Name << " created" << std::endl << std::endl;
			std::cout << std::endl << "-----------------------------------------------------------------------" << std::endl << std::endl;

			FileSample[WhichSample]->Close();
			TruthFileSample[WhichSample]->Close();

		} // End of the loop over the samples

	} // End of the loop over the runs	

} // End of the program 
