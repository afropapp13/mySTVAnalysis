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

void StandardEfficiencies(TString OverlaySample) {

	// -------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;
	gStyle->SetOptStat(0);	
	TGaxis::SetMaxDigits(4);
	TGaxis::SetExponentOffset(-0.05, 0., "y");	
	
	double TextSize = 0.07;

	// -------------------------------------------------------------------------------------

	int NEventsPassingSelectionCuts = 0;
	TString CutExtension = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	// v52
	VectorCuts.push_back("");
	VectorCuts.push_back("_PID");
	VectorCuts.push_back("_NuScore");

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

	PlotNames.push_back("VertexXPlot");
	PlotNames.push_back("VertexYPlot");
	PlotNames.push_back("VertexZPlot");

	PlotNames.push_back("EvPlot");
	PlotNames.push_back("NuPlot");

	PlotNames.push_back("MuonTrueMomentumLongitudinalRatio");
	PlotNames.push_back("ProtonTrueMomentumLongitudinalRatio");
	PlotNames.push_back("MuonTrueMomentumTransverseRatio");
	PlotNames.push_back("ProtonTrueMomentumTransverseRatio");

	PlotNames.push_back("CCQEMuonMomentumPlot"); 
	PlotNames.push_back("CCQEMuonCosThetaPlot"); 
	PlotNames.push_back("CCQEProtonMomentumPlot"); 
	PlotNames.push_back("CCQEProtonCosThetaPlot");

	// -------------------------------------------------------------------------------------

	const int N1DPlots = PlotNames.size();
	//cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ------------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");				
	Runs.push_back("Combined");				

	int NRuns = (int)(Runs.size());
	//cout << "Number of Runs = " << NRuns << endl;

	// -------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

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

		TString Name = "";

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			TString STVPath = PathToFiles+"/"+CutExtension+"/";
			TString STVName = "STVStudies_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+CutExtension+".root";
			FileSample.push_back(TFile::Open(STVPath+STVName));
			
			TString TrueSTVName = "TruthSTVAnalysis_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
			TruthFileSample.push_back(TFile::Open(TrueSTVPath+TrueSTVName));

			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();
			vector<TH1D*> CurrentPlotsTrueReco; CurrentPlotsTrueReco.clear();

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

				TH1D* histTrue = (TH1D*)(TruthFileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
				CurrentPlotsTrue.push_back(histTrue);

				TH1D* histTrueReco = (TH1D*)(FileSample[WhichSample]->Get("CC1pTrue"+PlotNames[WhichPlot]));
				CurrentPlotsTrueReco.push_back(histTrueReco);
		
			}

			PlotsTrue.push_back(CurrentPlotsTrue);
			PlotsTrueReco.push_back(CurrentPlotsTrueReco);

		}

		// Loop over the samples

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			TString EfficiencyName = "FileStandardEfficiences_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root"; 
			Name = FileEfficienciesPath + EfficiencyName;
			TFile* FileEfficiences = new TFile(Name,"recreate");

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				// Number of event distributions for True CC1p and Reco CC1p

/*				if (WhichSample == 0 && OverlaySample == "") {

					PlotsTrue[WhichSample][WhichPlot]->SetLineColor(kRed);
					PlotsTrue[WhichSample][WhichPlot]->SetLineWidth(3);
					PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->CenterTitle();
					PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
					PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
					PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->SetTitleSize(TextSize);
					PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(TextSize);
					PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->SetNdivisions(5);

					PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->CenterTitle();
					PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
					PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(TextSize);
					PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
					PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetRangeUser(0.,1.2*PlotsTrue[WhichSample][WhichPlot]->GetMaximum());
					PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(6);
					PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(TextSize);
					PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(TextSize);
					PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetTitle("# events");

					PlotsTrueReco[WhichSample][WhichPlot]->SetLineColor(kBlue);
					PlotsTrueReco[WhichSample][WhichPlot]->SetLineWidth(3);

					TString CanvasName = NameOfSamples[WhichSample]+"_"+PlotNames[WhichPlot];
					TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
					PlotCanvas->cd();
					PlotCanvas->SetBottomMargin(0.16);
					PlotCanvas->SetLeftMargin(0.15);								

					TLegend* leg = new TLegend(0.25,0.92,0.9,1.);
					leg->SetBorderSize(0);
					leg->SetTextSize(TextSize);
					leg->SetTextFont(FontStyle);
					leg->SetNColumns(2);
					leg->SetMargin(0.15);

					// ----------------------------------------------------------------------------------------

					PlotsTrue[WhichSample][WhichPlot]->Draw();
					PlotsTrueReco[WhichSample][WhichPlot]->Draw("same");

					leg->AddEntry(PlotsTrue[WhichSample][WhichPlot],"True CC1p");
					leg->AddEntry(PlotsTrueReco[WhichSample][WhichPlot],"Reco CC1p");
					leg->Draw();

					TLatex *text = new TLatex();
					text->SetTextFont(FontStyle);
					text->SetTextSize(0.08);
					text->DrawTextNDC(0.18, 0.8, Runs[WhichRun]);

					// ----------------------------------------------------------------------------------------
				
					TString CanvasPath = PlotPath+NameOfSamples[WhichSample]+"/";
					TString CanvasPdfName = PlotNames[WhichPlot]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".pdf";
					PlotCanvas->SaveAs(CanvasPath + CanvasPdfName); 

					delete PlotCanvas;
					
				}
*/

				// ---------------------------------------------------------------------------------------------------------------------------	

				// Ratios to extract the effective efficiencies			

				TH1D* pEffPlot = (TH1D*)PlotsTrueReco[WhichSample][WhichPlot]->Clone();
				pEffPlot->Divide(PlotsTrue[WhichSample][WhichPlot]);

				FileEfficiences->cd();
				pEffPlot->Write();

				if (WhichSample == 0 && OverlaySample == "") { 

					pEffPlot->SetLineWidth(2);
					pEffPlot->SetLineColor(kBlack);
					pEffPlot->SetMarkerStyle(20);

					pEffPlot->GetXaxis()->CenterTitle();
					pEffPlot->GetXaxis()->SetTitleFont(FontStyle);
					pEffPlot->GetXaxis()->SetLabelFont(FontStyle);
					pEffPlot->GetXaxis()->SetTitle(PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->GetTitle());
					pEffPlot->GetXaxis()->SetTitleSize(TextSize);
					pEffPlot->GetXaxis()->SetLabelSize(TextSize);
					pEffPlot->GetXaxis()->SetNdivisions(9);

					pEffPlot->GetYaxis()->CenterTitle();
					pEffPlot->GetYaxis()->SetTitleFont(FontStyle);
					pEffPlot->GetYaxis()->SetTitleSize(TextSize);
					pEffPlot->GetYaxis()->SetLabelFont(FontStyle);
					pEffPlot->GetYaxis()->SetNdivisions(6);
					pEffPlot->GetYaxis()->SetTitleSize(TextSize);
					pEffPlot->GetYaxis()->SetLabelSize(TextSize);
					pEffPlot->GetYaxis()->SetTitle("Efficiency [%]");
					pEffPlot->Scale(100.);
					pEffPlot->GetYaxis()->SetRangeUser(0.,1.2*pEffPlot->GetMaximum());

					TString CanvasEffName = NameOfSamples[WhichSample]+"_"+"StandardEff"+PlotNames[WhichPlot]+"_"+Runs[WhichRun];
					TCanvas* PlotEffCanvas = new TCanvas(CanvasEffName,CanvasEffName,205,34,1024,768);
					PlotEffCanvas->cd();
					PlotEffCanvas->SetBottomMargin(0.16);
					PlotEffCanvas->SetLeftMargin(0.18);

					PlotEffCanvas->cd();

					pEffPlot->SetMarkerColor(kBlack);
					pEffPlot->SetMarkerStyle(20);
					pEffPlot->SetMarkerSize(2.);
					pEffPlot->SetLineColor(kBlack);

					pEffPlot->Draw("e1");

					TLatex *textEff = new TLatex();
					textEff->SetTextFont(FontStyle);
					textEff->SetTextSize(TextSize);
					textEff->DrawTextNDC(0.22, 0.8, Runs[WhichRun]);
				
					TString CanvasEffPath = PlotPath+NameOfSamples[WhichSample]+"/";
					TString CanvasEffRatioName = "StandardEff"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".pdf";
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
