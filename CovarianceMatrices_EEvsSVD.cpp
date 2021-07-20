#include <TFile.h>
#include <TF1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TLegend.h>

#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "ubana/AnalysisCode/Secondary_Code/GlobalSettings.cpp"
#include "ubana/AnalysisCode/Secondary_Code/myFunctions.cpp"

void CovarianceMatrices_EEvsSVD(TString BaseMC = "Overlay9") {

	// -------------------------------------------------------------------------------------

	GlobalSettings();
	TGaxis::SetMaxDigits(3);

	// -------------------------------------------------------------------------------------

	vector<TString> PlotNames;
	PlotNames.push_back("DeltaPTPlot"); 
	PlotNames.push_back("DeltaAlphaTPlot"); 
	PlotNames.push_back("DeltaPhiTPlot"); 
	PlotNames.push_back("MuonMomentumPlot"); 
	PlotNames.push_back("MuonPhiPlot"); 
	PlotNames.push_back("MuonCosThetaPlot");
	PlotNames.push_back("ProtonMomentumPlot"); 
	PlotNames.push_back("ProtonPhiPlot"); 
	PlotNames.push_back("ProtonCosThetaPlot");

	const int N1DPlots = PlotNames.size();
		
	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
//	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");			
	Runs.push_back("Combined");

	const int NRuns = (int)(Runs.size());

	// -------------------------------------------------------------------------------------

	// XSec Samples

	vector<TFile*> EE_BeamOnFileSample; EE_BeamOnFileSample.resize(NRuns);
	vector<TFile*> SVD_BeamOnFileSample; SVD_BeamOnFileSample.resize(NRuns);

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	// XSec Plots

	vector <TH1D*> EE_BeamOnPlots; EE_BeamOnPlots.resize(N1DPlots);
	vector <TH1D*> SVD_BeamOnPlots; SVD_BeamOnPlots.resize(N1DPlots);

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	// Covariances

	vector<vector <TH2D*> > Covariances; Covariances.resize(NRuns,vector<TH2D*>(N1DPlots));

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------

		TString FileName = MigrationMatrixPath+"WienerSVD_UnfoldingTechnique_CovarianceMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TFile* FileCovarianceMatrices = new TFile(FileName,"recreate");
		
		TString SVD_FileName = PathToExtractedXSec+"WienerSVD_ExtractedXSec_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		SVD_BeamOnFileSample[WhichRun] = TFile::Open(SVD_FileName,"readonly");

		TString EE_FileName = PathToExtractedXSec+"ExtractedXSec_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		EE_BeamOnFileSample[WhichRun] = TFile::Open(EE_FileName,"readonly");

		// -------------------------------------------------------------------------------------

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

			// -------------------------------------------------------------------------------------------------------

			// Grab XSec plots

			SVD_BeamOnPlots[WhichPlot] = (TH1D*)(SVD_BeamOnFileSample[WhichRun]->Get("Reco"+PlotNames[WhichPlot]));
			EE_BeamOnPlots[WhichPlot] = (TH1D*)(EE_BeamOnFileSample[WhichRun]->Get("Reco"+PlotNames[WhichPlot]));

			// -------------------------------------------------------------------------------------------------------
			
			int NBins = SVD_BeamOnPlots[WhichPlot]->GetXaxis()->GetNbins();
			const double* ArrayBins = SVD_BeamOnPlots[WhichPlot]->GetXaxis()->GetXbins()->GetArray();
			TString XTitle = SVD_BeamOnPlots[WhichPlot]->GetXaxis()->GetTitle();

			// -------------------------------------------------------------------------------------------------------

			// Declare the matrix & initialize the entries to 0

			Covariances[WhichRun][WhichPlot] = new TH2D("UnfoldingTechnique_Covariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";i bin "+XTitle+";j bin "+XTitle,NBins,ArrayBins,NBins,ArrayBins);

			for (int WhichXBin = 0; WhichXBin < NBins; WhichXBin++) { 

				for (int WhichYBin = 0; WhichYBin < NBins; WhichYBin++) {

					Covariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,0.);
					Covariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,0.);

				}  // end of the loop over bin Y

			} // end of the loop over bin X

			// -------------------------------------------------------------------------------------------------------

			for (int WhichXBin = 0; WhichXBin < NBins; WhichXBin++) { 

				for (int WhichYBin = 0; WhichYBin < NBins; WhichYBin++) {

					double CovEntry = 0.;
					double CovError = 0.;

					if (WhichXBin == WhichYBin) {

						// X Bin entry / error

						double DataEntryX = SVD_BeamOnPlots[WhichPlot]->GetBinContent(WhichXBin+1);
						double DataErrorX = SVD_BeamOnPlots[WhichPlot]->GetBinError(WhichXBin+1);

						double AltDataEntryX = EE_BeamOnPlots[WhichPlot]->GetBinContent(WhichXBin+1);
						double AltDataErrorX = EE_BeamOnPlots[WhichPlot]->GetBinError(WhichXBin+1);

						// Y Bin entry / error

						double DataEntryY = SVD_BeamOnPlots[WhichPlot]->GetBinContent(WhichYBin+1);
						double DataErrorY = SVD_BeamOnPlots[WhichPlot]->GetBinError(WhichYBin+1);

						double AltDataEntryY = EE_BeamOnPlots[WhichPlot]->GetBinContent(WhichYBin+1);
						double AltDataErrorY = EE_BeamOnPlots[WhichPlot]->GetBinError(WhichYBin+1);

						// -------------------------------------------------------------------------------------------------------

						// Setting the elements of the Fractional Cov Matrix

						CovEntry = TMath::Max( ((AltDataEntryX - DataEntryX) / DataEntryX) * ( (AltDataEntryY - DataEntryY) / DataEntryY) / TMath::Sqrt(12),1E-8);

	//					LocalCovError = TMath::Max( TMath::Sqrt( 
	//						TMath::Power(DataEntryY - AltDataEntryY,2.) * ( TMath::Power(DataErrorX,2.) + TMath::Power(AltDataErrorX,2.) ) +
	//						TMath::Power(DataEntryX - AltDataEntryX,2.) * ( TMath::Power(DataErrorY,2.) + TMath::Power(AltDataErrorY,2.) ) ), 1E-10) ;

						CovError = 1E-8;

					}

					Covariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovEntry);
					Covariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovError);

					// -------------------------------------------------------------------------------------------------------

				}  // end of the loop over bin Y

			} // end of the loop over bin X
			
			FileCovarianceMatrices->cd();
			Covariances[WhichRun][WhichPlot]->Write();
			
			// ---------------------------------------------------------------------------------------	

			// Plot the total covariance matrix		
		
			TString CanvasName = "UnfoldingTechnique_"+PlotNames[WhichPlot]+BaseMC+"_"+Runs[WhichRun];
			TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.16);
			PlotCanvas->SetLeftMargin(0.15);
			PlotCanvas->SetRightMargin(0.25);			
			
			gStyle->SetMarkerSize(1.5);
			gStyle->SetPaintTextFormat("4.3f");			
			
			Covariances[WhichRun][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
			Covariances[WhichRun][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
			Covariances[WhichRun][WhichPlot]->GetXaxis()->SetTitleSize(TextSize);
			Covariances[WhichRun][WhichPlot]->GetXaxis()->SetLabelSize(TextSize);			
			Covariances[WhichRun][WhichPlot]->GetXaxis()->CenterTitle();
			Covariances[WhichRun][WhichPlot]->GetXaxis()->SetNdivisions(8);
			
			Covariances[WhichRun][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
			Covariances[WhichRun][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
			Covariances[WhichRun][WhichPlot]->GetYaxis()->SetTitleSize(TextSize);
			Covariances[WhichRun][WhichPlot]->GetYaxis()->SetLabelSize(TextSize);			
			Covariances[WhichRun][WhichPlot]->GetYaxis()->CenterTitle();
			Covariances[WhichRun][WhichPlot]->GetYaxis()->SetNdivisions(5);
			Covariances[WhichRun][WhichPlot]->GetYaxis()->SetTitleOffset(1.);						

			Covariances[WhichRun][WhichPlot]->SetTitle(Runs[WhichRun] + " Unfolding Technique");	

			double CovMax = TMath::Min(1.,1.05*Covariances[WhichRun][WhichPlot]->GetMaximum());
			double CovMin = TMath::Min(0.,1.05*Covariances[WhichRun][WhichPlot]->GetMinimum());
			Covariances[WhichRun][WhichPlot]->GetZaxis()->SetRangeUser(CovMin,CovMax);
//			Covariances[WhichRun][WhichPlot]->GetZaxis()->SetTitle("[x10^{-76} cm^{4}]");
			Covariances[WhichRun][WhichPlot]->GetZaxis()->CenterTitle();
			Covariances[WhichRun][WhichPlot]->GetZaxis()->SetTitleFont(FontStyle);
			Covariances[WhichRun][WhichPlot]->GetZaxis()->SetTitleSize(TextSize);
			Covariances[WhichRun][WhichPlot]->GetZaxis()->SetLabelFont(FontStyle);
			Covariances[WhichRun][WhichPlot]->GetZaxis()->SetLabelSize(TextSize-0.01);
			Covariances[WhichRun][WhichPlot]->GetZaxis()->SetNdivisions(5);

			Covariances[WhichRun][WhichPlot]->SetMarkerColor(kWhite);			
			Covariances[WhichRun][WhichPlot]->SetMarkerSize(1.5);
//			Covariances[WhichRun][WhichPlot]->Draw("text colz e"); 
			Covariances[WhichRun][WhichPlot]->Draw("colz");
			
			PlotCanvas->SaveAs(PlotPath+BaseMC+"/WienerSVD_UnfoldingTechnique_CovarianceMatrices_"+PlotNames[WhichPlot]+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
			
			delete PlotCanvas;	

			// ---------------------------------------------------------------------------------------			

		} // End of the loop over the plots

		// ---------------------------------------------------------------------------------------	

		// Close covariance matrix file & base files

		FileCovarianceMatrices->Close();

		SVD_BeamOnFileSample[WhichRun]->Close();
		EE_BeamOnFileSample[WhichRun]->Close();

		// ---------------------------------------------------------------------------------------	

		cout << endl << "File with UnfoldingTechnique Covariance matrices " << FileName << " has been created" << endl << endl;

		// ---------------------------------------------------------------------------------------	

	} // End of the loop over the runs	

} // End of the program
