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

#include "../myClasses/Constants.h"
#include "../myClasses/myFunctions.cpp"

using namespace std;
using namespace Constants;

void efficiency(TString OverlaySample, bool DetVar = false) {

	// -------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);	
	TGaxis::SetMaxDigits(4);
	TGaxis::SetExponentOffset(-0.05, 0., "y");	
	
	double TextSize = 0.07;

	// -------------------------------------------------------------------------------------

	int NEventsPassingSelectionCuts = 0;
	TString CutExtension = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	VectorCuts.push_back("");
	VectorCuts.push_back("_PID");
	VectorCuts.push_back("_NuScore");
	VectorCuts.push_back("_CRT");

	int NCuts = (int)(VectorCuts.size());	
	for (int i = 0; i < NCuts; i++) { CutExtension = CutExtension + VectorCuts[i]; }

	// -------------------------------------------------------------------------------------

	vector<TString> PlotNamesClone = PlotNames;

	// -------------------------------------------------------------------------------------

	const int N1DPlots = PlotNamesClone.size();
	//cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ------------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");				
	Runs.push_back("Combined");	

	// For runs 1-3, we used only run3 det vars

	/*if (DetVar) {

		Runs.clear();
		Runs.push_back("Run3");

	}*/			

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

				TH1D* histTrue = (TH1D*)(TruthFileSample[WhichSample]->Get("True"+PlotNamesClone[WhichPlot]));
				CurrentPlotsTrue.push_back(histTrue);

				TH1D* histTrueReco = (TH1D*)(FileSample[WhichSample]->Get("CC1pTrue"+PlotNamesClone[WhichPlot]));
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

					leg->AddEntry(PlotsTrue[WhichSample][WhichPlot],"True CC1p");
					leg->AddEntry(PlotsTrueReco[WhichSample][WhichPlot],"Reco CC1p");
					leg->Draw();

					TLatex *text = new TLatex();
					text->SetTextFont(FontStyle);
					text->SetTextSize(0.08);
					text->DrawTextNDC(0.18, 0.8, Runs[WhichRun]);

					// ----------------------------------------------------------------------------------------
				
					TString CanvasPath = PlotPath+NameOfSamples[WhichSample]+"/";
					TString CanvasPdfName = PlotNamesClone[WhichPlot]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".pdf";
					PlotCanvas->SaveAs(CanvasPath + CanvasPdfName); 

					delete PlotCanvas;
					
				}

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
					pEffPlot->GetXaxis()->SetTitleOffset(1.1);
					pEffPlot->GetXaxis()->SetLabelSize(TextSize);
					pEffPlot->GetXaxis()->SetLabelOffset(0.02);
					pEffPlot->GetXaxis()->SetNdivisions(9);
					bin_number_x_title(pEffPlot);

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

					TString CanvasEffName = NameOfSamples[WhichSample]+"_"+"StandardEff"+PlotNamesClone[WhichPlot]+"_"+Runs[WhichRun];
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
					if (Runs[WhichRun] != "Combined") { textEff->DrawTextNDC(0.22, 0.8, Runs[WhichRun]); }
					if ( !(string(PlotNamesClone[WhichPlot]).find("VertexZ") != std::string::npos) ) 
						{ textEff->DrawLatexNDC(0.22, 0.93, LatexLabel[PlotNamesClone[WhichPlot]]); }					
				
					//----------------------------------------//

					// Plot vertical lines
					// Add latex label with phase space limits

					if (string(PlotNames[WhichPlot]).find("Serial") != std::string::npos) {	

						TString clone_name = PlotNames[WhichPlot];
						clone_name.ReplaceAll("Reco","");
						vector<int> bin_break_points = get_2d_bin_break_points( map_to_2d_bin.at(clone_name) );

						int nbreaks = bin_break_points.size() - 1;
						vector<TLine*> line; line.resize(nbreaks);

						for (int ipoint = 0; ipoint < nbreaks; ipoint ++) {

							line.at(ipoint) = new TLine( bin_break_points.at(ipoint) + 0.5,0., bin_break_points.at(ipoint) + 0.5, pEffPlot->GetMaximum() );
							line.at(ipoint)->SetLineStyle(kDashed);
							line.at(ipoint)->Draw("same");

						}
	
						//----------------------------------------//

						vector<TLatex*> slice; slice.resize(nbreaks+1);

						for (int ipoint = 0; ipoint < nbreaks + 1; ipoint ++) {

			
							slice.at(ipoint) = new TLatex();
							slice.at(ipoint)->SetTextFont(FontStyle);
							slice.at(ipoint)->SetTextSize(0.02);
							TString phase_space = MapUncorCor[ clone_name + "_" + TString(std::to_string(ipoint) ) ];
							if (ipoint == 0) { slice.at(ipoint)->DrawLatex( bin_break_points.at(ipoint) / 5. , 0.9 * pEffPlot->GetMaximum(), LatexLabel[phase_space ]); }
							else { slice.at(ipoint)->DrawLatex( bin_break_points.at(ipoint - 1) + ( bin_break_points.at(ipoint) - bin_break_points.at(ipoint-1) ) / 5. , 0.9 * pEffPlot->GetMaximum(), LatexLabel[phase_space ]); }


						}

					}

					//----------------------------------------//
	
					TString CanvasEffPath = PlotPath+NameOfSamples[WhichSample]+"/";
					TString CanvasEffRatioName = "StandardEff"+PlotNamesClone[WhichPlot]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".pdf";
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
