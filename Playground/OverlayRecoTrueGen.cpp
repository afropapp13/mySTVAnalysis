#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TLatex.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "../../myClasses/Util.h"

//----------------------------------------//

TString ToString(double num) {

	std::ostringstream start;
	start << std::fixed << std::setprecision(2) << num;
	TString start1 = TString(start.str());
	return start1;

}

//----------------------------------------//

int LocateBinWithValue(TH1D* h, double Value) {

	int NBins = h->GetXaxis()->GetNbins();

	for (int i = 1; i <= NBins; i++) {

		double CurrentEntry = h->GetBinContent(i);
		if (CurrentEntry == Value) { return i; } 

	}

	return -99;

}

//----------------------------------------//

void PrettyPlot(TH1D* h) {

	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetTitleSize(0.06);
	h->GetXaxis()->SetLabelSize(0.06);
	h->GetXaxis()->SetTitleOffset(1.05);
	h->GetXaxis()->SetNdivisions(9);


	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetNdivisions(8);
	h->GetYaxis()->SetTitleOffset(1.2);
	h->GetYaxis()->SetTitleSize(0.06);
	h->GetYaxis()->SetLabelSize(0.06);

}

//----------------------------------------//

void OverlayRecoTrueGen() {

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(4);		

	TString PathToFiles = "../myXSec/";

	//----------------------------------------//

	vector<TString> PlotNames;
	PlotNames.push_back("DeltaPTPlot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		vector<vector<TH2D*> > MigrationMatrix; MigrationMatrix.clear();
		vector<vector<TH1D*> > PlotsTotalReco; PlotsTotalReco.clear();
		vector<vector<TH1D*> > PlotsNormOnly; PlotsNormOnly.clear();		
		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
		vector<vector<TH1D*> > PlotsAbsTrue; PlotsAbsTrue.clear();		

		gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();

		vector<TString> NameOfSamples; NameOfSamples.clear();
		vector<int> Colors; Colors.clear();		
		vector<TString> Labels; Labels.clear();

		// CV

		NameOfSamples.push_back("Overlay9"); Colors.push_back(OverlayColor); Labels.push_back("MC uB Tune");                     

		// -------------------------------------------------------------------------------------------------------------------

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();
		vector<TFile*> MigrationSample; MigrationSample.clear();		

		// -------------------------------------------------------------------------------------------------------------------

		// Open the files and grap the relevant plots

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			vector<TH2D*> CurrentMigrationMatrix; CurrentMigrationMatrix.clear();
			vector<TH1D*> CurrentPlotsTotalReco; CurrentPlotsTotalReco.clear();
			vector<TH1D*> CurrentPlotsNormOnly; CurrentPlotsNormOnly.clear();			
			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();
			vector<TH1D*> CurrentPlotsCC1pReco; CurrentPlotsCC1pReco.clear();
			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();
			vector<TH1D*> CurrentPlotsAbsTrue; CurrentPlotsAbsTrue.clear();


			// CV With Statistical Uncertainties

			if (NameOfSamples[WhichSample] == "Overlay9") { // CV with statistical uncertainties only for now

				TString FileSampleName = PathToFiles+UBCodeVersion+"/WienerSVD_ExtractedXSec_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; 
				FileSample.push_back(TFile::Open(FileSampleName,"readonly"));

				TString MigrationName = "../myMigrationMatrices/"+UBCodeVersion+"/FileMigrationMatrices_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; 
				MigrationSample.push_back(TFile::Open(MigrationName,"readonly"));				 

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TH1D* histTotalReco = (TH1D*)(FileSample[WhichSample]->Get("StatReco"+PlotNames[WhichPlot]));
					CurrentPlotsTotalReco.push_back(histTotalReco);

					TH2D* histMigrationMatrix = (TH2D*)(MigrationSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]+"2D"));
					CurrentMigrationMatrix.push_back(histMigrationMatrix);					

					TH1D* histNormOnly = (TH1D*)(FileSample[WhichSample]->Get("NormOnlyReco"+PlotNames[WhichPlot]));
					CurrentPlotsNormOnly.push_back(histNormOnly);					

					TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("RecoFullUnc"+PlotNames[WhichPlot]));
					CurrentPlotsReco.push_back(histReco);

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);

					TH1D* histAbsTrue = (TH1D*)(FileSample[WhichSample]->Get("NoSmearTrue"+PlotNames[WhichPlot]));
					CurrentPlotsAbsTrue.push_back(histAbsTrue);					
		
				}

			} else if (NameOfSamples[WhichSample] == "Overlay9NuWro") {

				TString FileSampleName = PathToFiles+UBCodeVersion+"/"+NameOfSamples[WhichSample]+"WienerSVD_ExtractedXSec_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; 
				FileSample.push_back(TFile::Open(FileSampleName,"readonly")); 

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					//TH1D* histTotalReco = (TH1D*)(FileSample[WhichSample]->Get("StatReco"+PlotNames[WhichPlot]));
					CurrentPlotsTotalReco.push_back(nullptr);
					CurrentPlotsNormOnly.push_back(nullptr);					

					TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
					CurrentPlotsReco.push_back(histReco);

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);
		
				}

			}

			else {

				if (
					NameOfSamples[WhichSample] == "Genie_v3_0_6_Out_Of_The_Box" || 
					NameOfSamples[WhichSample] == "Genie_v3_0_6_uB_Tune_1" || 
					NameOfSamples[WhichSample] == "Genie_v3_0_6_Nominal" || 
					NameOfSamples[WhichSample] == "Genie_v3_0_6_NoFSI" || 
					NameOfSamples[WhichSample] == "Genie_v3_0_6_NoRPA" || 
					NameOfSamples[WhichSample] == "Genie_v3_0_6_NoCoulomb" || 
					NameOfSamples[WhichSample] == "Genie_v3_0_6_hN2018" ||
					NameOfSamples[WhichSample] == "Genie_v3_0_6_RFG" ||  
					NameOfSamples[WhichSample] == "Genie_v3_0_6_EffSF" ||  
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

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TH1D* histTotalReco = nullptr;
					CurrentPlotsTotalReco.push_back(histTotalReco);

					TH1D* histNormOnly = nullptr;
					CurrentPlotsNormOnly.push_back(histNormOnly);					

					TH1D* histReco = nullptr;
					CurrentPlotsReco.push_back(histReco);

					TH1D* histCC1pReco = nullptr;
					CurrentPlotsCC1pReco.push_back(histCC1pReco);

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);
		
				}

			}

			MigrationMatrix.push_back(CurrentMigrationMatrix);
			PlotsTotalReco.push_back(CurrentPlotsTotalReco);
			PlotsNormOnly.push_back(CurrentPlotsNormOnly);					
			PlotsReco.push_back(CurrentPlotsReco);		
			PlotsCC1pReco.push_back(CurrentPlotsCC1pReco);
			PlotsTrue.push_back(CurrentPlotsTrue);
			PlotsAbsTrue.push_back(CurrentPlotsAbsTrue);					

		}

		// -----------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

			TH2D* Ac = (TH2D*)FileSample[0]->Get("Ac"+PlotNames[WhichPlot]);
			TH2D* Cov = (TH2D*)FileSample[0]->Get("UnfCov"+PlotNames[WhichPlot]);			

			// -----------------------------------------------------------------------------------------------------------------------------

			TH2D* CovClone = (TH2D*)(Cov->Clone()); 

			int n = Cov->GetXaxis()->GetNbins();

			for (int ix = 1; ix <= n; ix++) {

				for (int iy = 1; iy <= n; iy++) {

					double WidthX = Cov->GetXaxis()->GetBinWidth(ix);
					double WidthY = Cov->GetYaxis()->GetBinWidth(iy);

					double TwoDWidth = WidthX * WidthY;
					double BinContent = Cov->GetBinContent(ix,iy);

					CovClone->SetBinContent(ix,iy,BinContent/TwoDWidth);

				}					

			}	

			// -----------------------------------------------------------------------------------------------------------------------------			

			TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+"_"+Runs[WhichRun],PlotNames[WhichPlot]+"_"+Runs[WhichRun],205,34,1024,768);
			PlotCanvas->cd();

			TPad *midPad = new TPad("midPad", "",0.005, 0.2, 0.995, 0.995);
//			TPad *midPad = new TPad("midPad", "",0.005, 0., 0.995, 0.995);
			midPad->SetBottomMargin(0.14);
			midPad->SetTopMargin(0.12);
			midPad->SetLeftMargin(0.17);
			midPad->Draw();

			TPad *botPad = new TPad("botPad", "",0.005, 0., 0.995, 0.2);
			botPad->SetBottomMargin(0.14);
			botPad->SetTopMargin(0.);
			botPad->SetLeftMargin(0.17);
			botPad->Draw();

			TLegend* leg = new TLegend(0.5,0.58,0.63,0.85);
			if (PlotNames[WhichPlot] == "MuonPhiPlot" || PlotNames[WhichPlot] == "ProtonPhiPlot" || PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot") { leg = new TLegend(0.25,0.2,0.8,0.4); }
			if (PlotNames[WhichPlot] == "DeltaAlphaTPlot" || PlotNames[WhichPlot] == "MuonCosThetaPlot" || PlotNames[WhichPlot] == "ProtonCosThetaPlot" || PlotNames[WhichPlot] == "DeltaPtyPlot") 
				{ leg = new TLegend(0.2,0.58,0.5,0.85); }

			leg->SetBorderSize(0);
			leg->SetTextSize(0.04);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(1);
			if (PlotNames[WhichPlot] == "MuonPhiPlot" || PlotNames[WhichPlot] == "ProtonPhiPlot") { leg->SetNColumns(2); }
			leg->SetMargin(0.15);

			// ------------------------------------------------------------------------------------------------------------------

			// BeamOn Total Uncertainty

			PrettyPlot(PlotsReco[0][WhichPlot]);

			double MaxValue = PlotsReco[0][WhichPlot]->GetMaximum();
			int MaxValueBin = LocateBinWithValue(PlotsReco[0][WhichPlot],MaxValue);
			double MaxValueError = PlotsReco[0][WhichPlot]->GetBinError(MaxValueBin);

			double MinValue = PlotsReco[0][WhichPlot]->GetMinimum();
													
			PlotsReco[0][WhichPlot]->GetYaxis()->SetRangeUser(XSecRange[PlotNames[WhichPlot]].first,XSecRange[PlotNames[WhichPlot]].second);
			PlotsReco[0][WhichPlot]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);			

			PlotsReco[0][WhichPlot]->SetLineColor(kWhite);
			PlotsReco[0][WhichPlot]->SetMarkerColor(kWhite);
			PlotsReco[0][WhichPlot]->SetMarkerSize(1.);
			PlotsReco[0][WhichPlot]->SetMarkerStyle(20);

			midPad->cd();

			PlotsReco[0][WhichPlot]->Draw("e1x0 same"); // Total Unc

			PlotsTotalReco[0][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsTotalReco[0][WhichPlot]->SetMarkerColor(BeamOnColor);
			//PlotsTotalReco[0][WhichPlot]->Draw("e1x0 same"); // Stat Only

			// -----------------------------------------------------------------------------------------------------------------

			// Overlay Ac * Truth

			PrettyPlot(PlotsTrue[0][WhichPlot]);
			PlotsTrue[0][WhichPlot]->SetLineColor(Colors[0]);
			PlotsTrue[0][WhichPlot]->SetMarkerColor(Colors[0]);
			PlotsTrue[0][WhichPlot]->Draw("hist same");

			// -----------------------------------------------------------------------------------------------------------------

			// Overlay Truth

			PrettyPlot(PlotsAbsTrue[0][WhichPlot]);
			PlotsAbsTrue[0][WhichPlot]->SetLineColor(kOrange+7);
			PlotsAbsTrue[0][WhichPlot]->SetMarkerColor(kOrange+7);
			PlotsAbsTrue[0][WhichPlot]->Draw("hist same");			

			//----------------------------------------//

			// Grab the true signal and multiply it by the migration matrix

			int m = PlotsAbsTrue[0][WhichPlot]->GetNbinsX();
			double Nuedges[m+1];
			    
			for (int i = 0; i < m+1; i++) { Nuedges[i] = PlotsAbsTrue[0][WhichPlot]->GetBinLowEdge(i+1); }			

			TVectorD signal(m);
			TMatrixD migration(m, m);

			H2V(PlotsAbsTrue[0][WhichPlot], signal);
			H2M(MigrationMatrix[0][WhichPlot], migration, kTRUE); // X axis: True, Y axis: Reco	

			TH1D* RecoUnf = new TH1D("RecoUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";",n,Nuedges);
			TVectorD RecoUnfold = migration * signal;
			V2H(RecoUnfold, RecoUnf);

			RecoUnf->SetLineColor(kGreen-3);
			RecoUnf->SetMarkerColor(kGreen-3);	
			RecoUnf->Draw("hist same");									

			//----------------------------------------//			

			// Legend & Run / POT

			double tor860_wcut = -99.;
			if (Runs[WhichRun] == "Run1") { tor860_wcut = Fulltor860_wcut_Run1; }
			if (Runs[WhichRun] == "Run2") { tor860_wcut = Fulltor860_wcut_Run2; }
			if (Runs[WhichRun] == "Run3") { tor860_wcut = Fulltor860_wcut_Run3; }
			if (Runs[WhichRun] == "Run4") { tor860_wcut = Fulltor860_wcut_Run4; }
			if (Runs[WhichRun] == "Run5") { tor860_wcut = Fulltor860_wcut_Run5; }
			if (Runs[WhichRun] == "Combined") { tor860_wcut = Fulltor860_wcut_Combined; }
			//TString Label = ToString(tor860_wcut)+" POT";

			// ---------------------------------------------------------------------------------------------------------
			// ---------------------------------------------------------------------------------------------------------

			//PlotsReco[0][WhichPlot]->Draw("e1x0 same"); // BeamOn Stat Total

			//leg->AddEntry(PlotsReco[0][WhichPlot],"MicroBooNE Data","ep");
			//leg->AddEntry(PlotsReco[0][WhichPlot],Label,"");
			leg->AddEntry(PlotsTrue[0][WhichPlot],"Ac*Truth ("+ToString(PlotsTrue[0][WhichPlot]->Integral("Width"))+" x10^{-38}) cm^{2}","l");
			leg->AddEntry(PlotsAbsTrue[0][WhichPlot],"Truth ("+ToString(PlotsAbsTrue[0][WhichPlot]->Integral("Width"))+" x10^{-38}) cm^{2}","l");
			leg->AddEntry(RecoUnf,"Reco = Smear*Truth ("+ToString(RecoUnf->Integral("Width"))+" x10^{-38} cm^{2})","l");

			leg->Draw();

			// ----------------------------------------------------------------------------------------------

			// Ratio plot

			botPad->cd();
			botPad->cd();
			botPad->SetGridx();
			botPad->SetGridy();						

			TH1D* CloneTrue = (TH1D*)(PlotsAbsTrue[0][WhichPlot]->Clone());
			CloneTrue->Add(PlotsTrue[0][WhichPlot],-1);
			CloneTrue->Divide(PlotsTrue[0][WhichPlot]);	

			CloneTrue->GetXaxis()->SetLabelSize(0.);
			CloneTrue->GetXaxis()->SetTitleSize(0.);

			CloneTrue->GetYaxis()->CenterTitle();
			CloneTrue->GetYaxis()->SetTitle("Residual");
			CloneTrue->GetYaxis()->SetLabelSize(0.15);
			CloneTrue->GetYaxis()->SetTitleSize(0.2);												
			CloneTrue->GetYaxis()->SetTitleOffset(0.2);
			CloneTrue->GetYaxis()->SetRangeUser(-0.19,0.69);			

			CloneTrue->Draw("hist same");	

			TH1D* CloneReco = (TH1D*)(RecoUnf->Clone());
			CloneReco->Add(PlotsTrue[0][WhichPlot],-1);
			CloneReco->Divide(PlotsTrue[0][WhichPlot]);	

			CloneReco->Draw("hist same");									

			// ----------------------------------------------------------------------------------------------			

			// Saving the canvas with the data (total uncertainties) vs overlay & generator predictions

			//PlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/BeamOn9/RecoTrue_WienerSVD_Generator_TotalUnc_Data_XSections_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");

			//delete PlotCanvas;

			// ----------------------------------------------------------------------------------------------

		} // End of the loop over the plots

	} // End of the loop over the runs	

} // End of the program 
