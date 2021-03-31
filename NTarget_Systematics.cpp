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

#include "../Secondary_Code/CenterAxisTitle.cpp"
#include "../Secondary_Code/SetOffsetAndSize.cpp"
#include "../Secondary_Code/MakeMyPlotPretty.cpp"
#include "../Secondary_Code/myFunctions.cpp"

#include "../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

void NTarget_Systematics() {

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;

	// -----------------------------------------------------------------------------------------------------------------------------------------

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

	// ------------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");
//	Runs.push_back("Run2");
	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// ------------------------------------------------------------------------------------------------------------------------------

	// Covariance matrices for each run / sample / plot

	std::vector< std::vector< std::vector<TH2D*> > > NTargetCovarianceMatrix; NTargetCovarianceMatrix.clear(); 

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	cout << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// ------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------------------------------------------------------------------------

		TString FileName = PathToSystematics+"NTarget_Systematics_"+Runs[WhichRun]+".root";
		TFile* SystFile = new TFile(FileName,"recreate");

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

			FileSample.push_back(TFile::Open(PathToExtractedXSec+"/ExtractedXSec_Overlay9_"+Runs[WhichRun]+NameOfSamples[WhichSample]+"_"+UBCodeVersion+".root"));

			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

				TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
				CurrentPlotsReco.push_back(histReco);
		
			}

			PlotsReco.push_back(CurrentPlotsReco);		

		}

		// ------------------------------------------------------------------------------------------------------------------

		// Covariance matrices

		std::vector< std::vector<TH2D*> > RunNTargetCovarianceMatrix; RunNTargetCovarianceMatrix.clear();

		// -----------------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

			// -----------------------------------------------------------------------------------------------------------------------

			// Clone the reference plot

			TH1D* SystPlot = (TH1D*)PlotsReco[0][WhichPlot]->Clone();

			int NBins = SystPlot->GetXaxis()->GetNbins();

			for (int WhichBin = 1; WhichBin <= NBins; WhichBin++){

				SystPlot->SetBinContent(WhichBin,0.);

			}

			// ---------------------------------------------------------------------------------------------------------------------------
	
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

			// ------------------------------------------------------------------------------------------------------------

			// Covariance matrices

			std::vector<TH2D*> PlotRunNTargetCovarianceMatrix; PlotRunNTargetCovarianceMatrix.clear();

			// ----------------------------------------------------------------------------------------------------------------------------	

			double max = -99.; 

			// Drawing data plots using the efficiencies from CV & files for systematics

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				// -------------------------------------------------------------------------------------------------

				// Covariance matrix array for specific run / EventWeightLabel / plot

				double ArrayXSecDiff[NBins][NBins];
				// initialize 2D array to 0
				// https://stackoverflow.com/questions/3082914/c-compile-error-variable-sized-object-may-not-be-initialized
				memset( ArrayXSecDiff, 0, NBins*NBins*sizeof(double) );

				// ---------------------------------------------------------------------------------------------------

				PlotCanvas->cd();
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

				// Add NTarget Uncertainty to the statistical one in quadrature

				TH1D* ClonePlot = (TH1D*)PlotsReco[0][WhichPlot]->Clone();

				for (int WhichBin = 1; WhichBin <= NBins; WhichBin++) {

					double CurrentBinContent = PlotsReco[0][WhichPlot]->GetBinContent(WhichBin);
					double CurrentBinError = PlotsReco[0][WhichPlot]->GetBinError(WhichBin);
					double NTargetError = NTargetUncertainty*CurrentBinContent;
					double TotalError = TMath::Sqrt( NTargetError*NTargetError + CurrentBinError*CurrentBinError );

					ClonePlot->SetBinError(WhichBin,TotalError);

					SystPlot->SetBinContent(WhichBin,NTargetError);
					SystPlot->SetBinError(WhichBin,0.);

					// Covariance matrix
					// Uncorrelated NTarget systematic errors

					ArrayXSecDiff[WhichBin-1][WhichBin-1] = TMath::Power(NTargetError,2.);

				}

				ClonePlot->SetLineColor(kRed);
				midPad->cd();
				ClonePlot->Draw("e1x0 same");
				leg->AddEntry(ClonePlot,"Total","ep");

				PlotsReco[WhichSample][WhichPlot]->Draw("e1x0 same");
				leg->AddEntry(PlotsReco[WhichSample][WhichPlot],"Statistical","ep");

				// ----------------------------------------------------------------------------------------------------

				// Covariance matrices

				TString TMatrixName = "NTargetCoveriantMatrix_"+Runs[WhichRun]+"_"+PlotNames[WhichPlot];	
				TString CovTitleAndLabels = TString(PlotsReco[WhichSample][WhichPlot]->GetXaxis()->GetTitle())+" "+Runs[WhichRun]+";Bin # ;Bin #";
				TH2D* LocalMatrix = new TH2D("Local"+TMatrixName,CovTitleAndLabels,NBins,0.5,NBins-0.5,NBins,0.5,NBins-0.5);

				for (int WhichXBin = 0; WhichXBin < NBins; WhichXBin++) {

					for (int WhichYBin = 0; WhichYBin < NBins; WhichYBin++) {

						LocalMatrix->SetBinContent(WhichXBin+1,WhichYBin+1,ArrayXSecDiff[WhichXBin][WhichYBin]);
					
					}

				}	

				PlotRunNTargetCovarianceMatrix.push_back( LocalMatrix );	
			
				SystFile->cd();
				PlotRunNTargetCovarianceMatrix[WhichSample]->Write(TMatrixName);

				// ----------------------------------------------------------------------------------------------------

				// Covariance matrices
				// Store them in pdf format

				TCanvas* PlotCov = new TCanvas("Cov"+PlotNames[WhichPlot]+Runs[WhichRun],"Cov"+PlotNames[WhichPlot]+Runs[WhichRun],205,34,1024,768);
				PlotCov->cd();
				gStyle->SetTitleFont(FontStyle,"t");
				PlotCov->SetRightMargin(0.15);
				PlotRunNTargetCovarianceMatrix[WhichSample]->GetXaxis()->CenterTitle();
				PlotRunNTargetCovarianceMatrix[WhichSample]->GetYaxis()->CenterTitle();
				PlotRunNTargetCovarianceMatrix[WhichSample]->SetMarkerColor(kWhite);				
				PlotRunNTargetCovarianceMatrix[WhichSample]->SetMarkerSize(1.2);
				PlotRunNTargetCovarianceMatrix[WhichSample]->Draw("text coltz");
				PlotCov->SaveAs(PlotPath+"BeamOn9/CovMatrix_NTarget_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
				delete PlotCov;

				// ----------------------------------------------------------------------------------------------------

			} // End of the loop over the xsec samples

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

			TString Label = Runs[WhichRun] + " " +ToString(tor860_wcut)+" NTarget";

			latex.DrawLatexNDC(0.63,0.94, Label);	

			// ---------------------------------------------------------------------------------------------------

			// Saving the canvas where the CV & SystVar predictions have been overlaid

			PlotCanvas->SaveAs(PlotPath+"BeamOn9/NTarget_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
			delete PlotCanvas;

			// ---------------------------------------------------------------------------------------------------------------------

			// Store the extracted systematic uncertainty

			SystFile->cd();
			SystPlot->Write(PlotNames[WhichPlot]);

			// ---------------------------------------------------------------------------------------------------------------------

			// Covariance matrices

			RunNTargetCovarianceMatrix.push_back(PlotRunNTargetCovarianceMatrix);

		} // End of the loop over the plots

		// ----------------------------------------------------------------------------------

		cout << endl << "Systematics file " << FileName << " has been created" << endl << endl;

		// ----------------------------------------------------------------------------------------------------------------

		// Covariance matrices		

		NTargetCovarianceMatrix.push_back(RunNTargetCovarianceMatrix);	

		// ----------------------------------------------------------------------------------

	} // End of the loop over the runs	

	// --------------------------------------------------------------------------------------------------------------

} // End of the program 
