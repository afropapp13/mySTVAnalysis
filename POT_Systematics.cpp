#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TMath.h>
#include <TLatex.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include  "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/ToString.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/MakeMyPlotPretty.cpp"

#include "./Constants.h"

using namespace std;
using namespace Constants;

void POT_Systematics() {

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;

	TString PathToFiles = "myXSec/";

	double POTUncertainty = 0.02; // 2% POT Uncertainty

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

	PlotNames.push_back("DeltaPTPlot"); 
	PlotNames.push_back("DeltaAlphaTPlot"); 
	PlotNames.push_back("DeltaPhiTPlot");
	PlotNames.push_back("MuonMomentumPlot"); 
	PlotNames.push_back("MuonCosThetaPlot"); 
	PlotNames.push_back("MuonPhiPlot");
	PlotNames.push_back("ProtonMomentumPlot"); 
	PlotNames.push_back("ProtonCosThetaPlot");
	PlotNames.push_back("ProtonPhiPlot");
	PlotNames.push_back("ECalPlot");
	PlotNames.push_back("EQEPlot"); 
	PlotNames.push_back("Q2Plot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		TFile* file = new TFile("mySystematics/"+UBCodeVersion+"/POT_Systematics_"+Runs[WhichRun]+".root","recreate");

		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();

		gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();

		vector<TString> NameOfSamples; NameOfSamples.clear();
		std::vector<int> Colors; Colors.clear();
		std::vector<int> Markers; Markers.clear();
	
		NameOfSamples.push_back(""); // Reference plot
		Colors.push_back(kBlack); Markers.push_back(20);

		// Detector Variations

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			FileSample.push_back(TFile::Open(PathToFiles+UBCodeVersion+"/ExtractedXSec_Overlay9_"+Runs[WhichRun]+NameOfSamples[WhichSample]+"_"+UBCodeVersion+".root"));

			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

				TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
				CurrentPlotsReco.push_back(histReco);
		
			}

			PlotsReco.push_back(CurrentPlotsReco);		

		}

		// -------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

			// ------------------------------------------------------------------------------------------------------------------------------------------------------

			// Clone the reference plot

			TH1D* SystPlot = (TH1D*)PlotsReco[0][WhichPlot]->Clone();

			int NBins = SystPlot->GetXaxis()->GetNbins();

			for (int WhichBin = 1; WhichBin <= NBins; WhichBin++){

				SystPlot->SetBinContent(WhichBin,0.);

			}

			// ------------------------------------------------------------------------------------------------------------------------------------------------------
	
			TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+Runs[WhichRun],PlotNames[WhichPlot]+Runs[WhichRun],205,34,1024,768);
			PlotCanvas->cd();

			TPad *midPad = new TPad("midPad", "", 0.005, 0., 0.995, 0.995);
			midPad->SetTopMargin(0.12);
			midPad->SetBottomMargin(0.13);
			midPad->SetLeftMargin(0.17);
			midPad->Draw();

			TLegend* leg = new TLegend(0.03,0.89,0.57,0.99);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.06);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(3);

			// ----------------------------------------------------------------------------------------------------------------------------------------------

			double max = -99.; 

			// Drawing data plots using the efficiencies from CV & files for systematics

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				MakeMyPlotPretty(PlotsReco[WhichSample][WhichPlot]);
				PlotsReco[WhichSample][WhichPlot]->SetLineColor(Colors[WhichSample]);
				PlotsReco[WhichSample][WhichPlot]->SetMarkerStyle(Markers[WhichSample]);
				PlotsReco[WhichSample][WhichPlot]->SetMarkerColor(Colors[WhichSample]);
				PlotsReco[WhichSample][WhichPlot]->SetMarkerSize(2.);

				PlotsReco[WhichSample][WhichPlot]->GetXaxis()->SetTitleSize(0.06);
				PlotsReco[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(0.06);
				PlotsReco[WhichSample][WhichPlot]->GetXaxis()->SetTitleOffset(1.);

				PlotsReco[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(1.27);
				PlotsReco[WhichSample][WhichPlot]->GetYaxis()->SetTitle(PlotXAxis[WhichPlot]);
				PlotsReco[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
				PlotsReco[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(0.06);
				PlotsReco[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(0.06);
				PlotsReco[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(3);

				double LocalMax = PlotsReco[WhichSample][WhichPlot]->GetMaximum();
				max = TMath::Max(LocalMax,max);
				PlotsReco[0][WhichPlot]->GetYaxis()->SetRangeUser(0,1.2*max);

				midPad->cd();

				// Add POT Uncertainty to the statistical one in quadrature

				TH1D* ClonePlot = (TH1D*)PlotsReco[0][WhichPlot]->Clone();

				for (int WhichBin = 1; WhichBin <= NBins; WhichBin++) {

					double CurrentBinContent = PlotsReco[0][WhichPlot]->GetBinContent(WhichBin);
					double CurrentBinError = PlotsReco[0][WhichPlot]->GetBinError(WhichBin);
					double POTError = POTUncertainty*CurrentBinContent;
					double TotalError = TMath::Sqrt( POTError*POTError + CurrentBinError*CurrentBinError );

					ClonePlot->SetBinError(WhichBin,TotalError);

					SystPlot->SetBinContent(WhichBin,POTError);
					SystPlot->SetBinError(WhichBin,0.);

				}

				ClonePlot->SetLineColor(kRed);
				ClonePlot->Draw("e1x0 same");
				leg->AddEntry(ClonePlot,"Total","ep");

				PlotsReco[WhichSample][WhichPlot]->Draw("e1x0 same");
				leg->AddEntry(PlotsReco[WhichSample][WhichPlot],"Statistical","ep");

			} // End of the loop over the overlay sample

			leg->Draw();	

			TLatex latex;
			latex.SetTextFont(FontStyle);
			latex.SetTextSize(0.06);

			double tor860_wcut = -99.;

			if (Runs[WhichRun] == "Run1") { tor860_wcut = tor860_wcut_Run1; }
//			if (Runs[WhichRun] == "Run2") { tor860_wcut = tor860_wcut_Run2; }
//			if (Runs[WhichRun] == "Run3") { tor860_wcut = tor860_wcut_Run3; }
//			if (Runs[WhichRun] == "Run4") { tor860_wcut = tor860_wcut_Run4; }
//			if (Runs[WhichRun] == "Run5") { tor860_wcut = tor860_wcut_Run5; }

			TString Label = Runs[WhichRun] + " " +ToString(tor860_wcut)+" POT";

			latex.DrawLatexNDC(0.63,0.94, Label);	

			// -----------------------------------------------------------------------------------------------------------------------------------

			// Saving the canvas where the CV & SystVar predictions have been overlaid

			PlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/BeamOn9/POT_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
			delete PlotCanvas;

			// -------------------------------------------------------------------------------------------------------------------------------------------------

			// Store the extracted systematic uncertainty

			file->cd();
			SystPlot->Write(PlotNames[WhichPlot]);

		} // End of the loop over the plots

	} // End of the loop over the runs	

} // End of the program 
