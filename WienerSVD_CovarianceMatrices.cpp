#include <TFile.h>
#include <TF1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TGaxis.h>

#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "ubana/AnalysisCode/Secondary_Code/GlobalSettings.cpp"
#include "ubana/AnalysisCode/Secondary_Code/myFunctions.cpp"

// TString Syst = "NTarget"

void WienerSVD_CovarianceMatrices(TString Syst = "None",TString BaseMC = "Overlay9",TString BeamOnSample = "BeamOn9",TString BeamOffSample = "ExtBNB9",TString DirtSample = "OverlayDirt9") {

	// -------------------------------------------------------------------------------------

	if (Syst == "None") { cout << "What the heck are you doing ? Specify the systematic uncertainty that you want to obtain!" << endl; return ;}

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
	Runs.push_back("Run1");
//	Runs.push_back("Run2");
	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");			

	const int NRuns = (int)(Runs.size());

	// -------------------------------------------------------------------------------------

	// Base Samples

	vector<TFile*> MCFileSample; MCFileSample.resize(NRuns);
	vector<TFile*> BeamOnFileSample; BeamOnFileSample.resize(NRuns);
	vector<TFile*> BeamOffFileSample; BeamOffFileSample.resize(NRuns);
	vector<TFile*> DirtFileSample; DirtFileSample.resize(NRuns);

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	// Base Plots

	vector <TH1D*> CC1pPlots; CC1pPlots.resize(N1DPlots); 
	vector <TH1D*> NonCC1pPlots; NonCC1pPlots.resize(N1DPlots);
	vector <TH1D*> BeamOnPlots; BeamOnPlots.resize(N1DPlots);
	vector <TH1D*> BeamOffPlots; BeamOffPlots.resize(N1DPlots);
	vector <TH1D*> DirtPlots; DirtPlots.resize(N1DPlots);
	vector <TH2D*> Covariances; Covariances.resize(N1DPlots);

	// -------------------------------------------------------------------------------------------------------------------------------------
	
	// CV Flux File

	TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------
	
		double DataPOT = ReturnBeamOnRunPOT(Runs[WhichRun]);						
		
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface) * (SoftFidSurface / Nominal_UB_XY_Surface);

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------

		TString FileName = MigrationMatrixPath+"WienerSVD_"+Syst+"_CovarianceMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TFile* FileCovarianceMatrices = new TFile(FileName,"recreate");
		
		// Open base files

		TString ExactFileLocation = PathToFiles+CutExtension;
		MCFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+BaseMC+"_"+Runs[WhichRun]+CutExtension+".root","readonly");
		BeamOnFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOnSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");
		BeamOffFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOffSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");
		DirtFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+DirtSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");

		// -------------------------------------------------------------------------------------

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

			// Grab base plots

			CC1pPlots[WhichPlot] = (TH1D*)(MCFileSample[WhichRun]->Get("CC1pReco"+PlotNames[WhichPlot]));
			NonCC1pPlots[WhichPlot] = (TH1D*)(MCFileSample[WhichRun]->Get("NonCC1pReco"+PlotNames[WhichPlot]));
			BeamOnPlots[WhichPlot] = (TH1D*)(BeamOnFileSample[WhichRun]->Get("Reco"+PlotNames[WhichPlot]));
			BeamOffPlots[WhichPlot] = (TH1D*)(BeamOffFileSample[WhichRun]->Get("Reco"+PlotNames[WhichPlot]));
			DirtPlots[WhichPlot] = (TH1D*)(DirtFileSample[WhichRun]->Get("Reco"+PlotNames[WhichPlot]));

			// -------------------------------------------------------------------------------------------------------

			TH1D* DataPlot = (TH1D*)BeamOnPlots[WhichPlot]->Clone();

			if (BeamOnSample == "BeamOn9") { // Beam On Sample, need to subtract ALL backgrounds

				DataPlot->Add(NonCC1pPlots[WhichPlot],-1.);
				DataPlot->Add(BeamOffPlots[WhichPlot],-1.);
				DataPlot->Add(DirtPlots[WhichPlot],-1.);

			} else { // Fake Data Studies, working with MC CC1p signal events as BeamOn

				DataPlot = CC1pPlots[WhichPlot];

			}

			// -------------------------------------------------------------------------------------------------------
			
			int NBins = BeamOnPlots[WhichPlot]->GetXaxis()->GetNbins();
			const double* ArrayBins = BeamOnPlots[WhichPlot]->GetXaxis()->GetXbins()->GetArray();
			TString XTitle = BeamOnPlots[WhichPlot]->GetXaxis()->GetTitle();

			// Declare the matrix & initialize the entries to 0

			Covariances[WhichPlot] = new TH2D("Covariance_"+PlotNames[WhichPlot],";i bin "+XTitle+";j bin "+XTitle,NBins,ArrayBins,NBins,ArrayBins);

			for (int WhichXBin = 0; WhichXBin < NBins; WhichXBin++) { 

				for (int WhichYBin = 0; WhichYBin < NBins; WhichYBin++) {

					Covariances[WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,0.);
					Covariances[WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,0.);

				}  // end of the loop over bin Y

			} // end of the loop over bin X

			// -------------------------------------------------------------------------------------------------------

			for (int WhichXBin = 0; WhichXBin < NBins; WhichXBin++) { 

				for (int WhichYBin = 0; WhichYBin < NBins; WhichYBin++) {

					// X Bin entry / error

					double DataEntryX = DataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
					double DataErrorX = DataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

					double AltDataEntryX = DataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
					double AltDataErrorX = DataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

					// Y Bin entry / error

					double DataEntryY = DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
					double DataErrorY = DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

					double AltDataEntryY = DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
					double AltDataErrorY = DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

					// -------------------------------------------------------------------------------------------------------
					// -------------------------------------------------------------------------------------------------------

					// Based on the type of systematic unceratainty, choose how to handle it

					if (Syst == "NTarget") {

						AltDataEntryX = (1+NTargetUncertainty) * DataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorX = (1+NTargetUncertainty) * DataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

						AltDataEntryY = (1+NTargetUncertainty) * DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorY = (1+NTargetUncertainty) * DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

					} else if (Syst == "POT") {

						AltDataEntryX = (1+POTUncertainty) * DataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorX = (1+POTUncertainty) * DataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

						AltDataEntryY = (1+POTUncertainty) * DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorY = (1+POTUncertainty) * DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

					} else if (Syst == "Stat") {

						if (WhichXBin == WhichYBin) {

							DataEntryX = DataErrorX;
							DataEntryY = DataErrorY;

							AltDataEntryX = 0.;
							AltDataErrorX = 0.;

							AltDataEntryY = 0.;
							AltDataErrorY = 0.;

						} else {

							DataEntryX = 0.;
							DataEntryY = 0.;

							AltDataEntryX = 0.;
							AltDataErrorX = 0.;

							AltDataEntryY = 0.;
							AltDataErrorY = 0.;

						}

					}

					else { cout << "Don't know how to handle this systematic uncentainty! Abort !"<< endl << endl; return; }

					// -------------------------------------------------------------------------------------------------------
					// -------------------------------------------------------------------------------------------------------

					// CV

					double CovEntry = TMath::Max((AltDataEntryX - DataEntryX) * (AltDataEntryY - DataEntryY),1E-8);
					Covariances[WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovEntry);

					// Error

					double CovError = TMath::Max( TMath::Sqrt( 
							TMath::Power(DataEntryY - AltDataEntryY,2.) * ( TMath::Power(DataErrorX,2.) + TMath::Power(AltDataErrorX,2.) ) +
							TMath::Power(DataEntryX - AltDataEntryX,2.) * ( TMath::Power(DataErrorY,2.) + TMath::Power(AltDataErrorY,2.) ) ), 1E-10) ;

					Covariances[WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovError);

				}  // end of the loop over bin Y

			} // end of the loop over bin X
			
			FileCovarianceMatrices->cd();
			Covariances[WhichPlot]->Write();
			
			// ---------------------------------------------------------------------------------------	

			// Plot the total covariance matrix		
		
			TString CanvasName = Syst+"_"+PlotNames[WhichPlot]+BaseMC+"_"+Runs[WhichRun];
			TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.16);
			PlotCanvas->SetLeftMargin(0.15);
			PlotCanvas->SetRightMargin(0.25);			
			
			gStyle->SetMarkerSize(1.5);
			gStyle->SetPaintTextFormat("4.3f");			
			
			Covariances[WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
			Covariances[WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
			Covariances[WhichPlot]->GetXaxis()->SetTitleSize(TextSize);
			Covariances[WhichPlot]->GetXaxis()->SetLabelSize(TextSize);			
			Covariances[WhichPlot]->GetXaxis()->CenterTitle();
			Covariances[WhichPlot]->GetXaxis()->SetNdivisions(5);
			
			Covariances[WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
			Covariances[WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
			Covariances[WhichPlot]->GetYaxis()->SetTitleSize(TextSize);
			Covariances[WhichPlot]->GetYaxis()->SetLabelSize(TextSize);			
			Covariances[WhichPlot]->GetYaxis()->CenterTitle();
			Covariances[WhichPlot]->GetYaxis()->SetNdivisions(5);
			Covariances[WhichPlot]->GetYaxis()->SetTitleOffset(1.);						

			Covariances[WhichPlot]->SetTitle(Runs[WhichRun] + " " + Syst);	

			double CovMax = 1.05*Covariances[WhichPlot]->GetMaximum();
			double CovMin = TMath::Min(0.,1.05*Covariances[WhichPlot]->GetMinimum());
			Covariances[WhichPlot]->GetZaxis()->SetRangeUser(CovMin,CovMax);
			Covariances[WhichPlot]->GetZaxis()->SetTitle("[x10^{-76} cm^{4}]");
			Covariances[WhichPlot]->GetZaxis()->CenterTitle();
			Covariances[WhichPlot]->GetZaxis()->SetTitleFont(FontStyle);
			Covariances[WhichPlot]->GetZaxis()->SetTitleSize(TextSize);
			Covariances[WhichPlot]->GetZaxis()->SetLabelFont(FontStyle);
			Covariances[WhichPlot]->GetZaxis()->SetLabelSize(TextSize-0.01);
			Covariances[WhichPlot]->GetZaxis()->SetNdivisions(5);

			Covariances[WhichPlot]->SetMarkerColor(kWhite);			
			Covariances[WhichPlot]->SetMarkerSize(1.5);
//			Covariances[WhichPlot]->Draw("text colz e"); 
			Covariances[WhichPlot]->Draw("colz");
			
			PlotCanvas->SaveAs(PlotPath+BaseMC+"/WienerSVD_"+Syst+"_CovarianceMatrices_"+PlotNames[WhichPlot]+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
			
			delete PlotCanvas;			

		} // End of the loop over the plots

		// ---------------------------------------------------------------------------------------	

		// Close covariance matrix file & base files

		FileCovarianceMatrices->Close();

		MCFileSample[WhichRun]->Close();
		BeamOnFileSample[WhichRun]->Close();
		BeamOffFileSample[WhichRun]->Close();
		DirtFileSample[WhichRun]->Close();

		cout << endl << "File with "+Syst+" Covariance matrices " << FileName << " has been created" << endl << endl;

		// ---------------------------------------------------------------------------------------	

	} // End of the loop over the runs	

} // End of the program
