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

#include "../Secondary_Code/myFunctions.cpp"
#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "ubana/AnalysisCode/Secondary_Code/GlobalSettings.cpp"

void Systematics() {

	// ---------------------------------------------------------------------------------------------------------------------------

	GlobalSettings();

	// ---------------------------------------------------------------------------------------------------------------------------

	int DecimalAccuracy = 2;

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;

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

	PlotNames.push_back("CCQEMuonMomentumPlot"); 
	PlotNames.push_back("CCQEMuonCosThetaPlot"); 
	PlotNames.push_back("CCQEProtonMomentumPlot"); 
	PlotNames.push_back("CCQEProtonCosThetaPlot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ------------------------------------------------------------------------------------------------------------------------------

	//vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();

//		gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();

		vector<TString> NameOfSamples; NameOfSamples.clear();
	
		// CV

		NameOfSamples.push_back("Overlay9");

		// Systematics

		NameOfSamples.push_back("Detector_Systematics_LY");
		NameOfSamples.push_back("Detector_Systematics_TPC");
		NameOfSamples.push_back("POT_Systematics");
		NameOfSamples.push_back("NTarget_Systematics");
		NameOfSamples.push_back("G4_Systematics");
		NameOfSamples.push_back("Flux_Systematics");
		NameOfSamples.push_back("Genie_Systematics");		

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();

		// -------------------------------------------------------------------------------------------------------------------

		// Open the files and grap the relevant plots

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();
			vector<TH1D*> CurrentPlotsCC1pReco; CurrentPlotsCC1pReco.clear();
			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();

			// CV With Statistical Uncertainties

			if (WhichSample == 0) { // CV with statistical uncertainties

				FileSample.push_back(TFile::Open(PathToExtractedXSec+"ExtractedXSec_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root","update")); 

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
					CurrentPlotsReco.push_back(histReco);

					TH1D* histCC1pReco = (TH1D*)(FileSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]));
					CurrentPlotsCC1pReco.push_back(histCC1pReco);

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);
		
				}

			} 

			// Systematic Uncertainties

			else {

				FileSample.push_back(TFile::Open(PathToSystematics+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+".root")); 

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

					TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get(PlotNames[WhichPlot]));
					CurrentPlotsReco.push_back(histReco);

					TH1D* histCC1pReco = new TH1D();
					CurrentPlotsCC1pReco.push_back(histCC1pReco);

					TH1D* histTrue = new TH1D();
					CurrentPlotsTrue.push_back(histTrue);
		
				}

			}

			PlotsReco.push_back(CurrentPlotsReco);		
			PlotsCC1pReco.push_back(CurrentPlotsCC1pReco);
			PlotsTrue.push_back(CurrentPlotsTrue);		

		}

		// -----------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

			// Systematic Uncertainties

			TCanvas* PlotCanvasSyst = new TCanvas("Syst"+PlotNames[WhichPlot]+Runs[WhichRun],"Syst"+PlotNames[WhichPlot]+Runs[WhichRun],205,34,1024,768);
			PlotCanvasSyst->cd();

			TPad *midPadSyst = new TPad("midPadSyst", "",0.005, 0., 0.995, 0.995);
			midPadSyst->SetBottomMargin(0.14);
			midPadSyst->SetTopMargin(0.12);
			midPadSyst->SetLeftMargin(0.18);
			midPadSyst->Draw();

			TLegend* leg = new TLegend(0.17,0.89,0.55,0.99);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.06);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(2);

			// --------------------------------------------------------------------------------------------------------------------

			int NBinsX = PlotsReco[0][WhichPlot]->GetNbinsX();

			TH1D* SystDataplot = (TH1D*)(PlotsReco[0][WhichPlot]->Clone("TotalReco"+PlotNames[WhichPlot]));

			// Loop over bins to get the total uncertainty

			for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

				double CurrentDataEntry = PlotsReco[0][WhichPlot]->GetBinContent(WhichXBin+1); // CV
				double CurrentDataError = PlotsReco[0][WhichPlot]->GetBinError(WhichXBin+1); // Statistical Uncertainty
				double CurrentDataErrorSquared = CurrentDataError * CurrentDataError; // Statistical Uncertainty Squared

				double SystUncertaintySquared = 0;

				for (int WhichSystFile = 1; WhichSystFile < NSamples; WhichSystFile++) {

					double CurrentSystDataEntry = PlotsReco[WhichSystFile][WhichPlot]->GetBinContent(WhichXBin+1); // CV of SystFile

					SystUncertaintySquared += TMath::Power(CurrentSystDataEntry,2.);

				}
				

// Safety check
//if (SystUncertaintySquared > 10 * CurrentDataErrorSquared) { SystUncertaintySquared = 0.; }			

				double TotalUncertainty = sqrt(SystUncertaintySquared + CurrentDataErrorSquared); // Total Uncertainty

				SystDataplot->SetBinContent(WhichXBin+1,CurrentDataEntry);
				SystDataplot->SetBinError(WhichXBin+1,TotalUncertainty);

			} // End of the loop over the bins

			// ------------------------------------------------------------------------------------------------------------------

			// Make the canvas/plots pretty

			midPadSyst->cd();

			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetTitleSize(0.06);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetLabelSize(0.06);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetTitleOffset(1.05);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetNdivisions(5);


			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetNdivisions(4);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitleOffset(1.3);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitleSize(0.06);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetLabelSize(0.06);

			// -----------------------------------------------------------------------------------------------------------------

			// Overlay

			PlotsCC1pReco[0][WhichPlot]->SetLineColor(OverlayColor);
			PlotsCC1pReco[0][WhichPlot]->SetFillColor(OverlayColor);
			PlotsCC1pReco[0][WhichPlot]->Draw("e2 same");

			// BeamOn Statistical Uncertainty

			PlotsReco[0][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsReco[0][WhichPlot]->SetMarkerColor(BeamOnColor);
			PlotsReco[0][WhichPlot]->SetMarkerSize(1.5);
			PlotsReco[0][WhichPlot]->SetMarkerStyle(20);
			PlotsReco[0][WhichPlot]->Draw("e1x0 same");

			// BeamOn Total Uncertainty

			SystDataplot->SetLineColor(BeamOnColor);
			SystDataplot->SetMarkerColor(BeamOnColor);
			SystDataplot->SetMarkerSize(2.0);
			SystDataplot->SetMarkerStyle(20);
			SystDataplot->Draw("e1x0 same");

			// Genie Overlay for closure test

			PlotsTrue[0][WhichPlot]->SetLineColor(GenieColor);
			PlotsTrue[0][WhichPlot]->SetLineWidth(3);
			PlotsTrue[0][WhichPlot]->Draw("e same");

			// ---------------------------------------------------------------------------------------------------------

			// Legend & Run / POT

			double tor860_wcut = ReturnBeamOnRunPOT(Runs[WhichRun]);
			TString Label = ToString(tor860_wcut)+" POT";

			TLegendEntry* lMC = leg->AddEntry(PlotsCC1pReco[0][WhichPlot],"MC","f");
			lMC->SetTextColor(OverlayColor);

			TLegendEntry* lGenie = leg->AddEntry(PlotsTrue[0][WhichPlot],"GENIE Overlay (uB Tune v2)","l");
			lGenie->SetTextColor(GenieColor);

			leg->AddEntry(SystDataplot,"MicroBooNE Data " + Runs[WhichRun] + " " + Label,"ep");
			leg->Draw();

			// Add the total xsec (sta + syst) in root file

			FileSample[0]->cd();
			SystDataplot->Write();
			//FileSample[0]->Write("",TObject::kOverwrite);

			// ----------------------------------------------------------------------------------------------

			// Saving the canvas with the data (total uncertainties) vs overlay predictions

			PlotCanvasSyst->SaveAs(PlotPath+"BeamOn9/TotalUnc_Data_XSections_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");

			delete PlotCanvasSyst;

			// ------------------------------------------------------------------------------------------------------------------------

		} // End of the loop over the plots

		//FileSample[0]->Close();

	} // End of the loop over the runs	

} // End of the program 
