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

void Stat_Systematics() {

	TH1D::SetDefaultSumw2();
//	vector<TString> PlotNames;

	// -----------------------------------------------------------------------------------------------------------------------------------------

//	PlotNames.push_back("DeltaPTPlot"); 
//	PlotNames.push_back("DeltaAlphaTPlot"); 
//	PlotNames.push_back("DeltaPhiTPlot");
//	PlotNames.push_back("MuonMomentumPlot"); 
//	PlotNames.push_back("MuonCosThetaPlot"); 
//	PlotNames.push_back("MuonPhiPlot");
//	PlotNames.push_back("ProtonMomentumPlot"); 
//	PlotNames.push_back("ProtonCosThetaPlot");
//	PlotNames.push_back("ProtonPhiPlot");
//	PlotNames.push_back("ECalPlot");
//	PlotNames.push_back("EQEPlot"); 
//	PlotNames.push_back("Q2Plot");

//	PlotNames.push_back("CCQEMuonMomentumPlot"); 
//	PlotNames.push_back("CCQEMuonCosThetaPlot"); 
//	PlotNames.push_back("CCQEProtonMomentumPlot"); 
//	PlotNames.push_back("CCQEProtonCosThetaPlot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ------------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// ------------------------------------------------------------------------------------------------------------------------------

	// Covariance matrices for each run / sample / plot

	std::vector< std::vector< std::vector<TH2D*> > > StatCovarianceMatrix; StatCovarianceMatrix.clear(); 

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	cout << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// ------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------------------------------------------------------------------------

		TString FileName = PathToSystematics+"Stat_Systematics_"+Runs[WhichRun]+".root";
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

		std::vector< std::vector<TH2D*> > RunStatCovarianceMatrix; RunStatCovarianceMatrix.clear();

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
	
//			TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+Runs[WhichRun],PlotNames[WhichPlot]+Runs[WhichRun],205,34,1024,768);
//			PlotCanvas->cd();

//			TPad *midPad = new TPad("midPad", "", 0.005, 0., 0.995, 0.995);
//			midPad->SetTopMargin(0.12);
//			midPad->SetBottomMargin(0.13);
//			midPad->SetLeftMargin(0.17);
//			midPad->Draw();

//			TLegend* leg = new TLegend(0.03,0.89,0.57,0.99);
//			leg->SetBorderSize(0);
//			leg->SetTextSize(0.06);
//			leg->SetTextFont(FontStyle);
//			leg->SetNColumns(3);

			// ------------------------------------------------------------------------------------------------------------

			// Covariance matrices

			std::vector<TH2D*> PlotRunStatCovarianceMatrix; PlotRunStatCovarianceMatrix.clear();

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

				// Add Stat Uncertainty to the statistical one in quadrature

				for (int WhichBin = 1; WhichBin <= NBins; WhichBin++) {

					double CurrentBinError = PlotsReco[0][WhichPlot]->GetBinError(WhichBin);

					SystPlot->SetBinContent(WhichBin,CurrentBinError);
					SystPlot->SetBinError(WhichBin,0.);

					// Covariance matrix
					// Uncorrelated Stat systematic errors

					ArrayXSecDiff[WhichBin-1][WhichBin-1] = TMath::Power(CurrentBinError,2.);

				}

				// ----------------------------------------------------------------------------------------------------

				// Covariance matrices

				TString TMatrixName = "StatCoveriantMatrix_"+Runs[WhichRun]+"_"+PlotNames[WhichPlot];	
				TString CovTitleAndLabels = TString(PlotsReco[WhichSample][WhichPlot]->GetXaxis()->GetTitle())+" "+Runs[WhichRun]+";Bin # ;Bin #";

				TH2D* LocalMatrix = nullptr;
				if (NBins == 1) { LocalMatrix = new TH2D("Local"+TMatrixName,CovTitleAndLabels,1,0.5,1,1,0.5,1); }
				else { LocalMatrix = new TH2D("Local"+TMatrixName,CovTitleAndLabels,NBins,0.5,NBins-0.5,NBins,0.5,NBins-0.5); }

				for (int WhichXBin = 0; WhichXBin < NBins; WhichXBin++) {

					for (int WhichYBin = 0; WhichYBin < NBins; WhichYBin++) {

						LocalMatrix->SetBinContent(WhichXBin+1,WhichYBin+1,ArrayXSecDiff[WhichXBin][WhichYBin]);
					
					}

				}	

				PlotRunStatCovarianceMatrix.push_back( LocalMatrix );	
			
				SystFile->cd();
				PlotRunStatCovarianceMatrix[WhichSample]->Write(TMatrixName);

				// ----------------------------------------------------------------------------------------------------

				// Covariance matrices
				// Store them in pdf format

				TCanvas* PlotCov = new TCanvas("Cov"+PlotNames[WhichPlot]+Runs[WhichRun],"Cov"+PlotNames[WhichPlot]+Runs[WhichRun],205,34,1024,768);
				PlotCov->cd();
				PlotCov->SetRightMargin(0.15);
				gStyle->SetTitleFont(FontStyle,"t");
				PlotRunStatCovarianceMatrix[WhichSample]->GetXaxis()->CenterTitle();
				PlotRunStatCovarianceMatrix[WhichSample]->GetYaxis()->CenterTitle();
				PlotRunStatCovarianceMatrix[WhichSample]->SetMarkerColor(kWhite);				
				PlotRunStatCovarianceMatrix[WhichSample]->SetMarkerSize(1.2);
				PlotRunStatCovarianceMatrix[WhichSample]->Draw("text coltz");
				PlotCov->SaveAs(PlotPath+"BeamOn9/CovMatrix_Stat_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
				delete PlotCov;

				// ----------------------------------------------------------------------------------------------------

			} // End of the loop over the xsec samples

			// ---------------------------------------------------------------------------------------------------------------------

			// Store the extracted systematic uncertainty

			SystFile->cd();
			SystPlot->Write(PlotNames[WhichPlot]);

			// ---------------------------------------------------------------------------------------------------------------------

			// Covariance matrices

			RunStatCovarianceMatrix.push_back(PlotRunStatCovarianceMatrix);

		} // End of the loop over the plots

		// ----------------------------------------------------------------------------------

		cout << endl << "Systematics file " << FileName << " has been created" << endl << endl;

		// ----------------------------------------------------------------------------------------------------------------

		// Covariance matrices		

		StatCovarianceMatrix.push_back(RunStatCovarianceMatrix);	

		// ----------------------------------------------------------------------------------

	} // End of the loop over the runs	

	// --------------------------------------------------------------------------------------------------------------

} // End of the program 
