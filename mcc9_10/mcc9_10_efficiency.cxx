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

#include "../../../generators/constants.h"
#include "../../../generators/helper_functions.cxx"

using namespace std;
using namespace constants;

void mcc9_10_efficiency(TString OverlaySample, bool DetVar = false) {

	// -------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);	
	TGaxis::SetMaxDigits(4);
	TGaxis::SetExponentOffset(-0.05, 0., "y");	
	TH1::AddDirectory(kFALSE); 
	TH2::AddDirectory(kFALSE); 
	
	double TextSize = 0.07;

	// -------------------------------------------------------------------------------------

	int NEventsPassingSelectionCuts = 0;
	TString CutExtension = "_nocuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	VectorCuts.push_back("");
	//VectorCuts.push_back("_PID");
	
	int NCuts = (int)(VectorCuts.size());	
	for (int i = 0; i < NCuts; i++) { CutExtension = CutExtension + VectorCuts[i]; }

	// -------------------------------------------------------------------------------------

	vector<TString> PlotNamesClone = PlotNames;

	// -------------------------------------------------------------------------------------

	const int N1DPlots = PlotNamesClone.size();
	//cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ------------------------------------------------------------------------------------------------------------------------------------------	

	int NRuns = (int)(xsec_Runs.size());
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

		NameOfSamples.push_back("mcc9_10_Overlay9");

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();
		vector<TFile*> TruthFileSample; TruthFileSample.clear();

		TString Name = "";

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			TString STVPath = event_selection_file_path + "/"+CutExtension+"/";
			TString STVName = "ncpi0_"+NameOfSamples[WhichSample]+"_"+xsec_Runs[WhichRun]+OverlaySample+CutExtension+".root";
			FileSample.push_back(TFile::Open(STVPath+STVName));
			
			TString TrueSTVName = "Truthncpi0_"+NameOfSamples[WhichSample]+"_"+xsec_Runs[WhichRun]+OverlaySample+".root";
			TruthFileSample.push_back(TFile::Open(event_selection_file_path+TrueSTVName));

			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();
			vector<TH1D*> CurrentPlotsTrueReco; CurrentPlotsTrueReco.clear();

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

				TH1D* histTrue = (TH1D*)(TruthFileSample[WhichSample]->Get("True"+PlotNamesClone[WhichPlot]));
				CurrentPlotsTrue.push_back(histTrue);


				TH1D* histTrueReco = (TH1D*)(FileSample[WhichSample]->Get("NCCOHTrue"+PlotNamesClone[WhichPlot]));
				CurrentPlotsTrueReco.push_back(histTrueReco);
		
			}

			PlotsTrue.push_back(CurrentPlotsTrue);
			PlotsTrueReco.push_back(CurrentPlotsTrueReco);

		}

		// Loop over the samples

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			TString EfficiencyName = "FileStandardEfficiences_"+NameOfSamples[WhichSample]+"_"+xsec_Runs[WhichRun]+OverlaySample+".root"; 
			Name = efficiency_path + EfficiencyName;
			TFile* FileEfficiences = new TFile(Name,"recreate");

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				// Number of event distributions for True NCCOH and Reco NCCOH

				if (WhichSample == 0 && OverlaySample == "" && DetVar == false) {

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

					TString CanvasName = NameOfSamples[WhichSample]+"_"+PlotNamesClone[WhichPlot];
					TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
					PlotCanvas->cd();
					PlotCanvas->SetBottomMargin(0.17);
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

					leg->AddEntry(PlotsTrue[WhichSample][WhichPlot],"True NCCOH");
					leg->AddEntry(PlotsTrueReco[WhichSample][WhichPlot],"Reco NCCOH");
					leg->Draw();

					TLatex *text = new TLatex();
					text->SetTextFont(FontStyle);
					text->SetTextSize(0.08);
					text->DrawTextNDC(0.18, 0.8, xsec_Runs[WhichRun]);

					// ----------------------------------------------------------------------------------------
				
					TString CanvasPath = plot_path + NameOfSamples[WhichSample]+"/";
					TString CanvasPdfName = PlotNamesClone[WhichPlot]+"_"+xsec_Runs[WhichRun]+OverlaySample+".pdf";
					PlotCanvas->SaveAs(CanvasPath + CanvasPdfName); 

					delete PlotCanvas;
					
				}

				// ---------------------------------------------------------------------------------------------------------------------------	

				// ratios to extract the efficiencies			

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
					pEffPlot->GetXaxis()->SetTitleOffset(1.1);
					pEffPlot->GetXaxis()->SetLabelSize(TextSize);
					pEffPlot->GetXaxis()->SetLabelOffset(0.02);
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

					TString CanvasEffName = NameOfSamples[WhichSample]+"_"+"StandardEff"+PlotNamesClone[WhichPlot]+"_"+xsec_Runs[WhichRun];
					TCanvas* PlotEffCanvas = new TCanvas(CanvasEffName,CanvasEffName,205,34,1024,768);
					PlotEffCanvas->cd();
					PlotEffCanvas->SetTopMargin(0.11);
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
					textEff->DrawLatexNDC(0.22, 0.93, LatexLabel[PlotNamesClone[WhichPlot]]);				

					//----------------------------------------//
	
					TString CanvasEffPath = plot_path+NameOfSamples[WhichSample]+"/";
					TString CanvasEffRatioName = "StandardEff"+PlotNamesClone[WhichPlot]+"_"+xsec_Runs[WhichRun]+OverlaySample+".pdf";
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
