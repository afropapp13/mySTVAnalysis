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
#include <TGaxis.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../Secondary_Code/CenterAxisTitle.cpp"
#include "../Secondary_Code/SetOffsetAndSize.cpp"
#include "../Secondary_Code/MakeMyPlotPretty.cpp"
#include "../Secondary_Code/myFunctions.cpp"

#include "../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

void Detector_Systematics_LY() {

	TH1D::SetDefaultSumw2();
	TGaxis::SetMaxDigits(3);
	TGaxis::SetExponentOffset(-0.1, 1., "y");	
	
	vector<TString> PlotNames;

	// ---------------------------------------------------------------------------------------------------------------------------------------	

	TString PathToSystematics = "/uboone/data/users/"+UserID+"/mySTVAnalysis/mySystematics/"+UBCodeVersion+"/";
	TString PathToFiles = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myXSec/"+UBCodeVersion+"/";
	TString PlotPath = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myPlots/"+UBCodeVersion+"/"; 

	// ---------------------------------------------------------------------------------------------------------------------------------------	

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

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ---------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");
	Runs.push_back("Run2");
	Runs.push_back("Run3");
	Runs.push_back("Run4");
	Runs.push_back("Run5");				

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// ------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------------------------------------------------------------------------

		// To be removed when the rest of the runs are ready

		if (Runs[WhichRun] == "Run2") { continue; }
		if (Runs[WhichRun] == "Run4") { continue; }
		if (Runs[WhichRun] == "Run5") { continue; }

		// ------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------------------------------------------------------------------------

		TFile* file = new TFile(PathToSystematics+"/Detector_Systematics_LY_"+Runs[WhichRun]+".root","recreate");

		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();

		gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();

		vector<TString> NameOfSamples; NameOfSamples.clear();
		std::vector<int> Colors; Colors.clear();
		std::vector<int> Markers; Markers.clear();

		// ------------------------------------------------------------------------------------------------------------------
	
		NameOfSamples.push_back("CV"); // Reference plot
		Colors.push_back(kBlack); Markers.push_back(20);

		// Detector Variations

		// Runs 1 & 3

		NameOfSamples.push_back("LYDown"); Colors.push_back(kBlue); Markers.push_back(23);
		NameOfSamples.push_back("LYRayleigh"); Colors.push_back(kMagenta); Markers.push_back(29);
		
		// Run 3

//		if (Runs[WhichRun] == "Run3") {	

//			NameOfSamples.push_back("LYAttenuation"); Colors.push_back(kGreen+2); Markers.push_back(22);
//			NameOfSamples.push_back("WireModX"); Colors.push_back(kGreen+2); Markers.push_back(22);
//			NameOfSamples.push_back("WireModYZ"); Colors.push_back(kBlue); Markers.push_back(23);
//			NameOfSamples.push_back("WireModThetaYZ"); Colors.push_back(kMagenta); Markers.push_back(29);
//			NameOfSamples.push_back("WireModThetaXZ"); Colors.push_back(kOrange+7); Markers.push_back(47);
//	//		NameOfSamples.push_back("dEdx"); Colors.push_back(410); Markers.push_back(48);
//			NameOfSamples.push_back("Recombination2"); Colors.push_back(610); Markers.push_back(49);
//			NameOfSamples.push_back("SCE"); Colors.push_back(kCyan-7); Markers.push_back(33);

//		}	

		// ------------------------------------------------------------------------------------------------------------------				

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();

		TFile* NominalFile = TFile::Open(PathToFiles+"ExtractedXSec_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root");

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			FileSample.push_back(TFile::Open(PathToFiles+"ExtractedXSec_Overlay9_"+Runs[WhichRun]+"_"+NameOfSamples[WhichSample]+"_"+UBCodeVersion+".root"));

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

		// ----------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

			// ----------------------------------------------------------------------------------------------------------------

			// Nominal Plot

			TH1D* NominalPlot = (TH1D*)(NominalFile->Get("Reco"+PlotNames[WhichPlot]));

			// Clone the reference plot

			TH1D* SystPlot = (TH1D*)PlotsReco[0][WhichPlot]->Clone();

			int NBins = SystPlot->GetXaxis()->GetNbins();

			for (int WhichBin = 1; WhichBin <= NBins; WhichBin++){

				SystPlot->SetBinContent(WhichBin,0.);

			}

			// ----------------------------------------------------------------------------------------------------------------------
	
			TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+Runs[WhichRun],PlotNames[WhichPlot]+Runs[WhichRun],205,34,1024,768);
			PlotCanvas->cd();

			TPad *midPad = new TPad("midPad", "", 0.005, 0., 0.995, 0.995);
			midPad->SetTopMargin(0.13);
			midPad->SetBottomMargin(0.13);
			midPad->SetLeftMargin(0.17);
			midPad->Draw();

			TLegend* leg = new TLegend(0.0,0.87,0.98,0.98);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.06);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(3);

			// ------------------------------------------------------------------------------------------------------------------

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
//				PlotsReco[WhichSample][WhichPlot]->Draw("e same");
				leg->AddEntry(PlotsReco[WhichSample][WhichPlot],NameOfSamples[WhichSample],"p");

				// Take differences/2 to CV, store them and add them in quadrature

				TH1D* ClonePlot = (TH1D*)PlotsReco[0][WhichPlot]->Clone();
				ClonePlot->Add(PlotsReco[WhichSample][WhichPlot],-1);
				ClonePlot->Scale(0.5);

				for (int WhichBin = 1; WhichBin <= NBins; WhichBin++){

					double CurrentBinContent = SystPlot->GetBinContent(WhichBin);
					double ExtraBinContent = ClonePlot->GetBinContent(WhichBin);
					double CombinedBinContent = TMath::Sqrt( CurrentBinContent*CurrentBinContent + ExtraBinContent*ExtraBinContent );

					SystPlot->SetBinContent(WhichBin,CombinedBinContent);
					SystPlot->SetBinError(WhichBin,0.);

				}

			} // End of the loop over the detector variation samples 

//			midPad->cd();
//			MakeMyPlotPretty(NominalPlot);
//			NominalPlot->SetLineColor(kRed);
//			NominalPlot->SetMarkerStyle(22);
//			NominalPlot->SetMarkerColor(kRed);
//			NominalPlot->SetMarkerSize(2.);
//			NominalPlot->Draw("hist p0 same");
////			NominalPlot->Draw("e same");
//			leg->AddEntry(NominalPlot,"Nominal","p");

			leg->Draw();	

			TLatex latex;
			latex.SetTextFont(FontStyle);
			latex.SetTextSize(0.06);

			double tor860_wcut = -99.;

			if (Runs[WhichRun] == "Run1") { tor860_wcut = tor860_wcut_Run1; }
			if (Runs[WhichRun] == "Run2") { tor860_wcut = tor860_wcut_Run2; }
			if (Runs[WhichRun] == "Run3") { tor860_wcut = tor860_wcut_Run3; }
			if (Runs[WhichRun] == "Run4") { tor860_wcut = tor860_wcut_Run4; }
			if (Runs[WhichRun] == "Run5") { tor860_wcut = tor860_wcut_Run5; }

			TString Label = Runs[WhichRun] + " " +ToString(tor860_wcut)+" POT";

//			latex.DrawLatexNDC(0.63,0.94, Label);	
			latex.DrawLatexNDC(0.53,0.78, Label);				

			// -------------------------------------------------------------------------------------------------

			// Saving the canvas where the CV & SystVar predictions have been overlaid

			PlotCanvas->SaveAs(PlotPath+"BeamOn9/Detector_Variations_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");

			//delete PlotCanvas;

			// -----------------------------------------------------------------------------------------------------------

			// Store the extracted systematic uncertainty

			file->cd();
			SystPlot->Write(PlotNames[WhichPlot]);

		} // End of the loop over the plots

	} // End of the loop over the runs	

} // End of the program 
