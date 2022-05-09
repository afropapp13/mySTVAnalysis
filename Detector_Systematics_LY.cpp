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
	
//	vector<TString> PlotNames;
//	PlotNames.push_back("DeltaPTPlot"); 
//	PlotNames.push_back("DeltaAlphaTPlot"); 

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ---------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");				

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// ----------------------------------------------------------------------------------------------------------------------

	// Covariance matrices for each run / detector variation sample / plot

	std::vector< std::vector< std::vector<TH2D*> > > DetLYCovarianceMatrix; DetLYCovarianceMatrix.clear(); 

	// -----------------------------------------------------------------------------------------------------------------------

	cout << endl;

	// -----------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// ------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------------------------------------------------------------------------

		TString FileName = PathToSystematics+"Detector_Systematics_LY_"+Runs[WhichRun]+".root";
		TFile* SystFile = new TFile(FileName,"recreate");

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

		// -----------------------------------------------------------------------------------------				

		// LY Detector Variations

		NameOfSamples.push_back("LYDown"); Colors.push_back(kBlue); Markers.push_back(23);
		NameOfSamples.push_back("LYRayleigh"); Colors.push_back(kMagenta); Markers.push_back(29);
		NameOfSamples.push_back("LYAttenuation"); Colors.push_back(kRed+1); Markers.push_back(33);

		// ----------------------------------------------------------------------------------------

		// Covariance matrices

		std::vector< std::vector<TH2D*> > RunDetLYCovarianceMatrix; RunDetLYCovarianceMatrix.clear();

		// -----------------------------------------------------------------------------------------				

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();

		TFile* NominalFile = TFile::Open(PathToExtractedXSec+"ExtractedXSec_Overlay9_"+Runs[WhichRun]+"_CV_"+UBCodeVersion+".root");

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			FileSample.push_back(TFile::Open(PathToExtractedXSec+"ExtractedXSec_Overlay9_"+Runs[WhichRun]+"_"+NameOfSamples[WhichSample]+"_"+UBCodeVersion+".root"));

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

		} // End of the loop over the detector variation samples

		// ------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

			// -------------------------------------------------------------------------------------------------------

			// Covariance matrices

			std::vector<TH2D*> PlotRunDetLYCovarianceMatrix; PlotRunDetLYCovarianceMatrix.clear();

			// -------------------------------------------------------------------------------------------------------

			// Nominal Plot

			TH1D* NominalPlot = (TH1D*)(NominalFile->Get("Reco"+PlotNames[WhichPlot]));

			// Clone the reference plot

			TH1D* SystPlot = (TH1D*)PlotsReco[0][WhichPlot]->Clone();

			int NBins = SystPlot->GetXaxis()->GetNbins();

			// -----------------------------------------------------------------------------------------------------

			for (int WhichBin = 1; WhichBin <= NBins; WhichBin++){

				SystPlot->SetBinContent(WhichBin,0.);

			}

			// ----------------------------------------------------------------------------------------------------
	
			TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+Runs[WhichRun],PlotNames[WhichPlot]+Runs[WhichRun],205,34,1024,768);
			PlotCanvas->cd();

			TPad *midPad = new TPad("midPad", "", 0.005, 0., 0.995, 0.995);
			midPad->SetTopMargin(0.15);
			midPad->SetBottomMargin(0.13);
			midPad->SetLeftMargin(0.17);
			midPad->Draw();

			TLegend* leg = new TLegend(0.0,0.89,0.98,0.98);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.06);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(3);

			// --------------------------------------------------------------------------------------------------

			double max = -99.; 

			// Drawing data plots using the efficiencies from CV & files for systematics

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				// ------------------------------------------------------------------------------------------------

				// Covariance matrix array for specific run / detector variation sample / plot

				double ArrayXSecDiff[NBins][NBins];
				// initialize 2D array to 0
				// https://stackoverflow.com/questions/3082914/c-compile-error-variable-sized-object-may-not-be-initialized
				memset( ArrayXSecDiff, 0, NBins*NBins*sizeof(double) );

				// -----------------------------------------------------------------------------------------------

				MakeMyPlotPretty(PlotsReco[WhichSample][WhichPlot]);
				PlotsReco[WhichSample][WhichPlot]->SetLineColor(Colors[WhichSample]);
				PlotsReco[WhichSample][WhichPlot]->SetMarkerStyle(Markers[WhichSample]);
				PlotsReco[WhichSample][WhichPlot]->SetMarkerColor(Colors[WhichSample]);
				PlotsReco[WhichSample][WhichPlot]->SetMarkerSize(2.);

				PlotsReco[WhichSample][WhichPlot]->GetXaxis()->SetTitleSize(0.06);
				PlotsReco[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(0.06);
				PlotsReco[WhichSample][WhichPlot]->GetXaxis()->SetTitleOffset(1.);

				PlotsReco[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(1.27);
				PlotsReco[WhichSample][WhichPlot]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);
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

				for (int WhichBin = 1; WhichBin <= NBins; WhichBin++) {

					double CurrentBinContent = SystPlot->GetBinContent(WhichBin);
					double ExtraBinContent = ClonePlot->GetBinContent(WhichBin);
					double CombinedBinContent = TMath::Sqrt( CurrentBinContent*CurrentBinContent + ExtraBinContent*ExtraBinContent );

					SystPlot->SetBinContent(WhichBin,CombinedBinContent);
					SystPlot->SetBinError(WhichBin,0.);

					// Covariance matrix
					// Loop over all the other bin entries & take the relevant differences 

					double XSecDiffBin = ( PlotsReco[0][WhichPlot]->GetBinContent(WhichBin) - PlotsReco[WhichSample][WhichPlot]->GetBinContent(WhichBin) ) / 2.;

					for (int WhichOtherBin = 0; WhichOtherBin < NBins; WhichOtherBin++) {

						// Covariance Matrix
						// Take the xsec difference in loop over other bins

						double CVSampleOtherBin = PlotsReco[0][WhichPlot]->GetBinContent(WhichOtherBin+1);
						double VariationSampleOtherBin = PlotsReco[WhichSample][WhichPlot]->GetBinContent(WhichOtherBin+1);
						double XSecDiffOtherBin = (CVSampleOtherBin - VariationSampleOtherBin) / 2.;

						// Multisim approach, don't forget to divide by the number of universes
						double ArrayXSecDiffEntry = XSecDiffBin * XSecDiffOtherBin;

						ArrayXSecDiff[WhichBin-1][WhichOtherBin] += ArrayXSecDiffEntry;

					}

				}

				// ----------------------------------------------------------------------------------------------------

				// Covariance matrices

				TString TMatrixName = "DetLYCoveriantMatrix_"+Runs[WhichRun]+"_"+NameOfSamples[WhichSample]+"_"+PlotNames[WhichPlot];	
				TString CovTitleAndLabels = TString(PlotsReco[0][WhichPlot]->GetXaxis()->GetTitle())+" "+Runs[WhichRun]+";Bin # ;Bin #";

				TH2D* LocalMatrix = nullptr;
				if (NBins == 1) { LocalMatrix = new TH2D("Local"+TMatrixName,CovTitleAndLabels,1,0.5,1,1,0.5,1); }
				else { LocalMatrix = new TH2D("Local"+TMatrixName,CovTitleAndLabels,NBins,0.5,NBins-0.5,NBins,0.5,NBins-0.5); }

				for (int WhichXBin = 0; WhichXBin < NBins; WhichXBin++) {

					for (int WhichYBin = 0; WhichYBin < NBins; WhichYBin++) {

						LocalMatrix->SetBinContent(WhichXBin+1,WhichYBin+1,ArrayXSecDiff[WhichXBin][WhichYBin]);
					
					}

				}	

				PlotRunDetLYCovarianceMatrix.push_back( LocalMatrix );

				SystFile->cd();
				PlotRunDetLYCovarianceMatrix[WhichSample]->Write(TMatrixName);

				// ----------------------------------------------------------------------------------------------------

				// Covariance matrices
				// Store them in pdf format

				TCanvas* PlotCov = new TCanvas("Cov"+PlotNames[WhichPlot]+Runs[WhichRun],"Cov"+PlotNames[WhichPlot]+Runs[WhichRun],205,34,1024,768);
				PlotCov->cd();
				PlotCov->SetRightMargin(0.15);
				PlotRunDetLYCovarianceMatrix[WhichSample]->GetXaxis()->CenterTitle();
				PlotRunDetLYCovarianceMatrix[WhichSample]->GetYaxis()->CenterTitle();
				PlotRunDetLYCovarianceMatrix[WhichSample]->SetMarkerColor(kWhite);				
				PlotRunDetLYCovarianceMatrix[WhichSample]->SetMarkerSize(1.2);
				PlotRunDetLYCovarianceMatrix[WhichSample]->Draw("text coltz");
				PlotCov->SaveAs(PlotPath+"BeamOn9/CovMatrix_DetLY_"+PlotNames[WhichPlot]+"_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
				delete PlotCov;

				// -------------------------------------------------------------------------------------------

			} // End of the loop over the detector variation samples 		

			// -------------------------------------------------------------------------------------------

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

			if (Runs[WhichRun] == "Run1") { tor860_wcut = Fulltor860_wcut_Run1; }
			if (Runs[WhichRun] == "Run2") { tor860_wcut = Fulltor860_wcut_Run2; }
			if (Runs[WhichRun] == "Run3") { tor860_wcut = Fulltor860_wcut_Run3; }
			if (Runs[WhichRun] == "Run4") { tor860_wcut = Fulltor860_wcut_Run4; }
			if (Runs[WhichRun] == "Run5") { tor860_wcut = Fulltor860_wcut_Run5; }
			if (Runs[WhichRun] == "Combined") { tor860_wcut = Fulltor860_wcut_Combined; }			

			TString Label = Runs[WhichRun] + " " +ToString(tor860_wcut)+" POT";

//			latex.DrawLatexNDC(0.63,0.94, Label);	
			latex.DrawLatexNDC(0.53,0.78, Label);				

			// -------------------------------------------------------------------------------------------------

			// Saving the canvas where the CV & SystVar predictions have been overlaid

			PlotCanvas->SaveAs(PlotPath+"BeamOn9/LY_Detector_Variations_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");

			delete PlotCanvas;

			// -----------------------------------------------------------------------------------------------------------

			// Store the extracted systematic uncertainty

			SystFile->cd();
			SystPlot->Write(PlotNames[WhichPlot]);

			// ----------------------------------------------------------------------------------

			// Covariance matrices

			RunDetLYCovarianceMatrix.push_back( PlotRunDetLYCovarianceMatrix );

			// Covariance matrices	
			// Sum of detector variation contributions
			// To obtain total covariance matrix

			TString TMatrixName = "DetLYCoveriantMatrix_"+Runs[WhichRun]+"_"+NameOfSamples[0]+"_"+PlotNames[WhichPlot];
			TH2D* OverallLYDetCovMatrix = (TH2D*)(SystFile->Get(TMatrixName));

			for (int WhichSample = 1; WhichSample < NSamples; WhichSample ++) {	

				TString LocalTMatrixName = "DetLYCoveriantMatrix_"+Runs[WhichRun]+"_"+NameOfSamples[WhichSample]+"_"+PlotNames[WhichPlot];	
				TH2D* LocalDetLYCovMatrix = (TH2D*)(SystFile->Get(LocalTMatrixName));
				OverallLYDetCovMatrix->Add(LocalDetLYCovMatrix);

			}	

			OverallLYDetCovMatrix->Write("OverallDetLYCovMatrix_"+PlotNames[WhichPlot]);

			// ----------------------------------------------------------------------------------------------------

			// Overall Covariance matrix
			// Store them in pdf format

			TCanvas* OverallPlotCov = new TCanvas("OverallCov"+PlotNames[WhichPlot]+Runs[WhichRun],"Cov"+PlotNames[WhichPlot]+Runs[WhichRun],205,34,1024,768);
			OverallPlotCov->cd();
			OverallPlotCov->SetRightMargin(0.15);
			OverallLYDetCovMatrix->GetXaxis()->CenterTitle();
			OverallLYDetCovMatrix->GetYaxis()->CenterTitle();
			OverallLYDetCovMatrix->SetMarkerColor(kWhite);
			OverallLYDetCovMatrix->SetMarkerSize(1.2);
			OverallLYDetCovMatrix->Draw("text coltz");
			OverallPlotCov->SaveAs(PlotPath+"BeamOn9/OverallCovMatrix_DetLY_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
			delete OverallPlotCov;

			// ----------------------------------------------------------------------------------

		} // End of the loop over the plots

		// --------------------------------------------------------------------------------------------------------

		// Covariance matrices
			
		DetLYCovarianceMatrix.push_back(RunDetLYCovarianceMatrix);

		// ----------------------------------------------------------------------------------------------------

		cout << endl << "Systematics file " << FileName << " has been created" << endl << endl;

		// ----------------------------------------------------------------------------------------------------

	} // End of the loop over the runs	

} // End of the program 
