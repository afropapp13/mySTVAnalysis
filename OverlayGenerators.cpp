#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TEfficiency.h>
#include <TMath.h>
#include <TLatex.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>
//#include <math>

#include  "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/myFunctions.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/MakeMyPlotPretty.cpp"

#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

void OverlayGenerators() {

	int DecimalAccuracy = 2;

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;

	TString PathToFiles = "myXSec/";

	// ---------------------------------------------------------------------------------------------------------------------------

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

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
//	Runs.push_back("Run1");
	Runs.push_back("Run3");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		vector<vector<TH1D*> > PlotsTotalReco; PlotsTotalReco.clear();
		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();

		gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();

		vector<TString> NameOfSamples; NameOfSamples.clear();
	
		// CV

		NameOfSamples.push_back("Overlay9");

		NameOfSamples.push_back("Genie_v3_0_6_Out_Of_The_Box");
		NameOfSamples.push_back("Genie_v3_0_6_uB_Tune_1");
		NameOfSamples.push_back("SuSav2");
		NameOfSamples.push_back("NuWro");
		NameOfSamples.push_back("GiBUU");
		NameOfSamples.push_back("GENIEv2");
		NameOfSamples.push_back("NEUT");
		NameOfSamples.push_back("GENIEv3_0_4");

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();

		// -------------------------------------------------------------------------------------------------------------------

		// Open the files and grap the relevant plots

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			vector<TH1D*> CurrentPlotsTotalReco; CurrentPlotsTotalReco.clear();
			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();
			vector<TH1D*> CurrentPlotsCC1pReco; CurrentPlotsCC1pReco.clear();
			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();

			// CV With Statistical Uncertainties

			if (NameOfSamples[WhichSample] == "Overlay9") { // CV with statistical uncertainties only for now

				FileSample.push_back(TFile::Open(PathToFiles+UBCodeVersion+"/ExtractedXSec_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root","readonly")); 

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TH1D* histTotalReco = (TH1D*)(FileSample[WhichSample]->Get("TotalReco"+PlotNames[WhichPlot]));
					CurrentPlotsTotalReco.push_back(histTotalReco);

					TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
					CurrentPlotsReco.push_back(histReco);

					TH1D* histCC1pReco = (TH1D*)(FileSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]));
					CurrentPlotsCC1pReco.push_back(histCC1pReco);

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);
		
				}

			} 

			else {

				if (
					NameOfSamples[WhichSample] == "Genie_v3_0_6_Out_Of_The_Box" || 
					NameOfSamples[WhichSample] == "Genie_v3_0_6_uB_Tune_1" || 
					NameOfSamples[WhichSample] == "SuSav2" ||
					NameOfSamples[WhichSample] == "GENIEv2" ||
					NameOfSamples[WhichSample] == "GENIEv3_0_4"
				) {
					FileSample.push_back(TFile::Open("../myGenieAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root")); 
				}

				if (NameOfSamples[WhichSample] == "NuWro") 
					{ FileSample.push_back(TFile::Open("../myNuWroAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root")); }

				if (NameOfSamples[WhichSample] == "GiBUU") 
					{ FileSample.push_back(TFile::Open("../myGiBUUAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root")); }

				if (NameOfSamples[WhichSample] == "NEUT") 
					{ FileSample.push_back(TFile::Open("../myNEUTAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root")); }

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

					TH1D* histTotalReco = nullptr;
					CurrentPlotsTotalReco.push_back(histTotalReco);

					TH1D* histReco = nullptr;
					CurrentPlotsReco.push_back(histReco);

					TH1D* histCC1pReco = nullptr;
					CurrentPlotsCC1pReco.push_back(histCC1pReco);

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);
		
				}

			}

			PlotsTotalReco.push_back(CurrentPlotsTotalReco);		
			PlotsReco.push_back(CurrentPlotsReco);		
			PlotsCC1pReco.push_back(CurrentPlotsCC1pReco);
			PlotsTrue.push_back(CurrentPlotsTrue);		

		}

		// -----------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

			TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+"_"+Runs[WhichRun],PlotNames[WhichPlot]+"_"+Runs[WhichRun],205,34,1024,768);
			PlotCanvas->cd();

			TPad *midPad = new TPad("midPad", "",0.005, 0., 0.995, 0.995);
			midPad->SetBottomMargin(0.14);
			midPad->SetTopMargin(0.12);
			midPad->SetLeftMargin(0.16);
			midPad->Draw();

			TLegend* leg = new TLegend(0.02,0.89,0.98,0.99);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.04);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(4);

			// ------------------------------------------------------------------------------------------------------------------

			// Make the canvas/plots pretty

			midPad->cd();

			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetTitleSize(0.06);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetLabelSize(0.06);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetTitleOffset(1.05);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetNdivisions(5);


			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetNdivisions(4);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitleOffset(1.2);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitleSize(0.06);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetLabelSize(0.06);

			// -----------------------------------------------------------------------------------------------------------------

			// Overlay

			PlotsCC1pReco[0][WhichPlot]->SetLineColor(OverlayColor);
			PlotsCC1pReco[0][WhichPlot]->SetFillColor(OverlayColor);
			PlotsCC1pReco[0][WhichPlot]->Draw("e2 same");		

			// -----------------------------------------------------------------------------------------------------------------

			// BeamOn Statistical Uncertainty

			PlotsReco[0][WhichPlot]->SetLineWidth(5);
			PlotsReco[0][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsReco[0][WhichPlot]->SetMarkerColor(BeamOnColor);
			PlotsReco[0][WhichPlot]->SetMarkerSize(1.5);
			PlotsReco[0][WhichPlot]->SetMarkerStyle(20);
//			PlotsReco[0][WhichPlot]->Draw("e1x0 same");

			// -----------------------------------------------------------------------------------------------------------------

			// BeamOn Total Uncertainty

			PlotsTotalReco[0][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsTotalReco[0][WhichPlot]->SetMarkerColor(BeamOnColor);
			PlotsTotalReco[0][WhichPlot]->SetMarkerSize(1.5);
			PlotsTotalReco[0][WhichPlot]->SetMarkerStyle(20);
//			PlotsTotalReco[0][WhichPlot]->Draw("e1x0 same");

			// -----------------------------------------------------------------------------------------------------------------

			// Genie Overlay

			PlotsTrue[0][WhichPlot]->SetLineColor(GenieOverlayColor);
			PlotsTrue[0][WhichPlot]->SetLineWidth(3);
			PlotsTrue[0][WhichPlot]->Draw("e same");

			// -----------------------------------------------------------------------------------------------------------------

			// Genie v3.0.6 Out Of The Box

			PlotsTrue[1][WhichPlot]->SetLineColor(Geniev3OutOfTheBoxColor);
			PlotsTrue[1][WhichPlot]->SetLineWidth(3);
			PlotsTrue[1][WhichPlot]->Draw("e same");

			// -----------------------------------------------------------------------------------------------------------------

			// Genie v3.0.6 MicroBooNE Tune v1

			PlotsTrue[2][WhichPlot]->SetLineColor(GenieColor);
			PlotsTrue[2][WhichPlot]->SetLineWidth(3);
			PlotsTrue[2][WhichPlot]->Draw("e same");

			// -----------------------------------------------------------------------------------------------------------------

			// Genie SuSav2

			PlotsTrue[3][WhichPlot]->SetLineColor(SuSav2Color);
			PlotsTrue[3][WhichPlot]->SetLineWidth(3);
			PlotsTrue[3][WhichPlot]->Draw("e same");

			// -----------------------------------------------------------------------------------------------------------------

			// NuWro

			PlotsTrue[4][WhichPlot]->SetLineColor(NuWroColor);
			PlotsTrue[4][WhichPlot]->SetLineWidth(3);
			PlotsTrue[4][WhichPlot]->Draw("e same");

			// -----------------------------------------------------------------------------------------------------------------

			// GiBUU

			PlotsTrue[5][WhichPlot]->SetLineColor(GiBUUColor);
			PlotsTrue[5][WhichPlot]->SetLineWidth(3);
			PlotsTrue[5][WhichPlot]->Draw("e same");

			// -----------------------------------------------------------------------------------------------------------------

			// GENIE v2

			PlotsTrue[6][WhichPlot]->SetLineColor(GENIEv2Color);
			PlotsTrue[6][WhichPlot]->SetLineWidth(3);
			PlotsTrue[6][WhichPlot]->Draw("e same");

			// -----------------------------------------------------------------------------------------------------------------

			// NEUT

			PlotsTrue[7][WhichPlot]->SetLineColor(NEUTColor);
			PlotsTrue[7][WhichPlot]->SetLineWidth(3);
			PlotsTrue[7][WhichPlot]->Draw("e same");

			// -----------------------------------------------------------------------------------------------------------------

			// GENIE v3.0.4

			PlotsTrue[8][WhichPlot]->SetLineColor(GENIEv3_0_4_Color);
			PlotsTrue[8][WhichPlot]->SetLineWidth(3);
			PlotsTrue[8][WhichPlot]->Draw("e same");

			// ---------------------------------------------------------------------------------------------------------

			PlotsReco[0][WhichPlot]->Draw("e1x0 same"); // BeamOn Stat Unc Only

			// ---------------------------------------------------------------------------------------------------------

			// Legend & Run / POT

			double tor860_wcut = -99.;
			if (Runs[WhichRun] == "Run1") { tor860_wcut = tor860_wcut_Run1; }
			if (Runs[WhichRun] == "Run2") { tor860_wcut = tor860_wcut_Run2; }
			if (Runs[WhichRun] == "Run3") { tor860_wcut = tor860_wcut_Run3; }
			if (Runs[WhichRun] == "Run4") { tor860_wcut = tor860_wcut_Run4; }
			if (Runs[WhichRun] == "Run5") { tor860_wcut = tor860_wcut_Run5; }
			TString Label = ToString(tor860_wcut)+" POT";

			TLegendEntry* lGenie_GENIEv2 = leg->AddEntry(PlotsTrue[6][WhichPlot],"v2.12.10","l");
			lGenie_GENIEv2->SetTextColor(GENIEv2Color);

			TLegendEntry* lGenie_GENIEv3_0_4 = leg->AddEntry(PlotsTrue[8][WhichPlot],"v3.0.4","l");
			lGenie_GENIEv3_0_4->SetTextColor(GENIEv3_0_4_Color);

			TLegendEntry* lGenie = leg->AddEntry(PlotsTrue[1][WhichPlot],"v3.0.6","l");
			lGenie->SetTextColor(Geniev3OutOfTheBoxColor);

			TLegendEntry* lGenie_uBTunev1 = leg->AddEntry(PlotsTrue[2][WhichPlot],"v3.0.6 (uB Tune v1)","l");
			lGenie_uBTunev1->SetTextColor(GenieColor);

			TLegendEntry* lGenie_SuSav2 = leg->AddEntry(PlotsTrue[3][WhichPlot],"SuSav2","l");
			lGenie_SuSav2->SetTextColor(SuSav2Color);

			TLegendEntry* lGenie_NuWro = leg->AddEntry(PlotsTrue[4][WhichPlot],"NuWro","l");
			lGenie_NuWro->SetTextColor(NuWroColor);

			TLegendEntry* lGenie_GiBUU = leg->AddEntry(PlotsTrue[5][WhichPlot],"GiBUU","l");
			lGenie_GiBUU->SetTextColor(GiBUUColor);

			TLegendEntry* lGenie_NEUT = leg->AddEntry(PlotsTrue[7][WhichPlot],"NEUT","l");
			lGenie_NEUT->SetTextColor(NEUTColor);

			TLegendEntry* lGenie_GenieOverlay = leg->AddEntry(PlotsTrue[0][WhichPlot],"Genie Overlay","l");
			lGenie_GenieOverlay->SetTextColor(GenieOverlayColor);

			TLegendEntry* lMC = leg->AddEntry(PlotsCC1pReco[0][WhichPlot],"MC","f");
			lMC->SetTextColor(OverlayColor);

//			leg->AddEntry(PlotsReco[0][WhichPlot],"MicroBooNE Data " + Runs[WhichRun] + " " + Label,"ep");
			leg->AddEntry(PlotsReco[0][WhichPlot],"MicroBooNE Data","ep");
			leg->AddEntry(PlotsReco[0][WhichPlot],Runs[WhichRun] + " " + Label,"");

			leg->Draw();

			// ----------------------------------------------------------------------------------------------

			// Saving the canvas with the data (stat uncertainties) vs overlay & generator predictions

			PlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/BeamOn9/Generator_Data_XSections_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");

			// ----------------------------------------------------------------------------------------------

			// Saving the canvas with the data (total uncertainties) vs overlay & generator predictions

			PlotsTotalReco[0][WhichPlot]->Draw("e1x0 same"); // BeamOn Total Unc
			PlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/BeamOn9/Generator_TotalUnc_Data_XSections_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");

			// ----------------------------------------------------------------------------------------------

		} // End of the loop over the plots

	} // End of the loop over the runs	

} // End of the program 
