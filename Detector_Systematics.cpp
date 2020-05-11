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

void Detector_Systematics() {

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;

	TString PathToFiles = "myXSec/";

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

		TFile* file = new TFile("mySystematics/"+UBCodeVersion+"/Detector_Systematics_"+Runs[WhichRun]+".root","recreate");

		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();

		gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();

		vector<TString> NameOfSamples; NameOfSamples.clear();
		std::vector<int> Colors; Colors.clear();
		std::vector<int> Markers; Markers.clear();
	
		NameOfSamples.push_back("CV"); // Reference plot
		Colors.push_back(kBlack); Markers.push_back(20);

		// Detector Variations

		NameOfSamples.push_back("X"); Colors.push_back(kOrange+7); Markers.push_back(21);
		NameOfSamples.push_back("YZ"); Colors.push_back(kGreen+2); Markers.push_back(22);
		NameOfSamples.push_back("LY"); Colors.push_back(kBlue); Markers.push_back(23);
		NameOfSamples.push_back("LYRayleigh"); Colors.push_back(kMagenta); Markers.push_back(29);

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			FileSample.push_back(TFile::Open(PathToFiles+UBCodeVersion+"/ExtractedXSec_Overlay9_"+Runs[WhichRun]+"_"+NameOfSamples[WhichSample]+"_"+UBCodeVersion+".root"));

			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();
			vector<TH1D*> CurrentPlotsCC1pReco; CurrentPlotsCC1pReco.clear();
			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

				TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
				CurrentPlotsReco.push_back(histReco);

				TH1D* histCC1pReco = (TH1D*)(FileSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]));
				CurrentPlotsCC1pReco.push_back(histCC1pReco);

				TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
				CurrentPlotsTrue.push_back(histTrue);
		
			}

			PlotsReco.push_back(CurrentPlotsReco);		
			PlotsCC1pReco.push_back(CurrentPlotsCC1pReco);
			PlotsTrue.push_back(CurrentPlotsTrue);		

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
				PlotsReco[WhichSample][WhichPlot]->Draw("hist p0 same");
				leg->AddEntry(PlotsReco[WhichSample][WhichPlot],NameOfSamples[WhichSample],"p");

				// Take differences to CV, store them and add them in quadrature

				TH1D* ClonePlot = (TH1D*)PlotsReco[0][WhichPlot]->Clone();
				ClonePlot->Add(PlotsReco[WhichSample][WhichPlot],-1);

				for (int WhichBin = 1; WhichBin <= NBins; WhichBin++){

					double CurrentBinContent = SystPlot->GetBinContent(WhichBin);
					double ExtraBinContent = ClonePlot->GetBinContent(WhichBin);
					double CombinedBinContent = TMath::Sqrt( CurrentBinContent*CurrentBinContent + ExtraBinContent*ExtraBinContent );

					SystPlot->SetBinContent(WhichBin,CombinedBinContent);
					SystPlot->SetBinError(WhichBin,0.);

				}

			} // End of the loop over the detector variation samples 

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

			PlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/BeamOn9/Detector_Variations_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");

			delete PlotCanvas;

			// -------------------------------------------------------------------------------------------------------------------------------------------------

			// Store the extracted systematic uncertainty

			file->cd();
			SystPlot->Write(PlotNames[WhichPlot]);

		} // End of the loop over the plots

	} // End of the loop over the runs	

} // End of the program 
