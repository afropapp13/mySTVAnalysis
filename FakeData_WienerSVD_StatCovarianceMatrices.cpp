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
#include <TMatrixD.h>
#include <TVectorD.h>

#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "ubana/AnalysisCode/Secondary_Code/GlobalSettings.cpp"
#include "ubana/AnalysisCode/Secondary_Code/myFunctions.cpp"

#include "ubana/myClasses/Util.h"

//----------------------------------------//

void FakeData_WienerSVD_StatCovarianceMatrices(TString Syst = "None",TString BaseMC = "Overlay9",TString BeamOnSample = "BeamOn9",TString BeamOffSample = "ExtBNB9",TString DirtSample = "OverlayDirt9", TString Tune = "") {

	// -------------------------------------------------------------------------------------

	if (Syst == "None") { cout << "What the heck are you doing ? Specify the smearing / efficiency systematic uncertainty that you want to obtain!" << endl; return ;}

	// -------------------------------------------------------------------------------------

	GlobalSettings();
	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	TGaxis::SetMaxDigits(3);

	// -------------------------------------------------------------------------------------

//	vector<TString> PlotNames;
//	PlotNames.push_back("DeltaPTPlot"); 
//	PlotNames.push_back("DeltaAlphaTPlot"); 
//	PlotNames.push_back("DeltaPhiTPlot"); 
//	PlotNames.push_back("MuonMomentumPlot"); 
//	PlotNames.push_back("MuonPhiPlot"); 
//	PlotNames.push_back("MuonCosThetaPlot");
//	PlotNames.push_back("MuonCosThetaSingleBinPlot");
//	PlotNames.push_back("ProtonMomentumPlot"); 
//	PlotNames.push_back("ProtonPhiPlot"); 
//	PlotNames.push_back("ProtonCosThetaPlot");

//	PlotNames.push_back("CCQEMuonMomentumPlot"); 
//	PlotNames.push_back("CCQEMuonCosThetaPlot"); 
//	PlotNames.push_back("CCQEProtonMomentumPlot"); 
//	PlotNames.push_back("CCQEProtonCosThetaPlot");

	const int N1DPlots = PlotNames.size();
		
	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");			
	Runs.push_back("Combined");

	const int NRuns = (int)(Runs.size());

	// -------------------------------------------------------------------------------------

	// Base Samples

	vector<TFile*> TrueMCFileSample; TrueMCFileSample.resize(NRuns);
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
	vector<vector <TH2D*> > FracCovariances; FracCovariances.resize(NRuns,vector<TH2D*>(N1DPlots));
	vector<vector <TH2D*> > Covariances; Covariances.resize(NRuns,vector<TH2D*>(N1DPlots));

	// -------------------------------------------------------------------------------------------------------------------------------------
	
	// CV Flux File

	TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------
	
		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);						
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface);

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------

		TString FileName = MigrationMatrixPath+Tune+"WienerSVD_"+Syst+"_CovarianceMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		if (BeamOnSample != "BeamOn9") { FileName = MigrationMatrixPath+BeamOnSample+"WienerSVD_"+Syst+"_CovarianceMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }

		TFile* FileCovarianceMatrices = new TFile(FileName,"recreate");
		
		// Open base files

		TString ExactFileLocation = PathToFiles+CutExtension;
		TString TStringBaseMC = ExactFileLocation+"/"+Tune+"STVStudies_"+BaseMC+"_"+Runs[WhichRun]+CutExtension+".root";
		TString TrueTStringBaseMC = PathToFiles+"/"+Tune+"TruthSTVAnalysis_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TString ResponseFileName = MigrationMatrixPath+Tune+"FileResponseMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";

		TFile* FileResponseMatrices = new TFile(ResponseFileName,"readonly");
		TrueMCFileSample[WhichRun] = TFile::Open(TrueTStringBaseMC,"readonly");
		MCFileSample[WhichRun] = TFile::Open(TStringBaseMC,"readonly");
		BeamOffFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOffSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");
		DirtFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/"+Tune+"STVStudies_"+DirtSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");

		if (BeamOnSample == "NoTuneOverlay9") 
			{ BeamOnFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/NoTuneSTVStudies_Overlay9_"+Runs[WhichRun]+CutExtension+".root","readonly"); }
		else { BeamOnFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOnSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly"); }

		// -------------------------------------------------------------------------------------

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

			// -------------------------------------------------------------------------------------------------------

			// Grab base plots

			CC1pPlots[WhichPlot] = (TH1D*)(MCFileSample[WhichRun]->Get("CC1pReco"+PlotNames[WhichPlot]));
			NonCC1pPlots[WhichPlot] = (TH1D*)(MCFileSample[WhichRun]->Get("NonCC1pReco"+PlotNames[WhichPlot]));
			BeamOnPlots[WhichPlot] = (TH1D*)(BeamOnFileSample[WhichRun]->Get("Reco"+PlotNames[WhichPlot]));
			BeamOffPlots[WhichPlot] = (TH1D*)(BeamOffFileSample[WhichRun]->Get("Reco"+PlotNames[WhichPlot]));
			DirtPlots[WhichPlot] = (TH1D*)(DirtFileSample[WhichRun]->Get("Reco"+PlotNames[WhichPlot]));

			// -------------------------------------------------------------------------------------------------------

			TH1D* DataPlot = (TH1D*)BeamOnPlots[WhichPlot]->Clone();
			DataPlot->Add(BeamOffPlots[WhichPlot],-1.);
			DataPlot->Add(DirtPlots[WhichPlot],-1.);
			DataPlot->Add(NonCC1pPlots[WhichPlot],-1.);	

			// -------------------------------------------------------------------------------------------------------
			// -------------------------------------------------------------------------------------------------------
			
			int NBins = BeamOnPlots[WhichPlot]->GetXaxis()->GetNbins();
			const double* ArrayBins = BeamOnPlots[WhichPlot]->GetXaxis()->GetXbins()->GetArray();
			TString XTitle = BeamOnPlots[WhichPlot]->GetXaxis()->GetTitle();

			// -------------------------------------------------------------------------------------------------------

			// Declare the matrix & initialize the entries to 0

			if (PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot") {

				FracCovariances[WhichRun][WhichPlot] = new TH2D(Syst+"_FracCovariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";i bin "+XTitle+";j bin "+XTitle,1,-1.,1.,1,-1.,1.);
				Covariances[WhichRun][WhichPlot] = new TH2D(Syst+"_Covariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";i bin "+XTitle+";j bin "+XTitle,1,-1.,1.,1,-1.,1.);

			} else {

				FracCovariances[WhichRun][WhichPlot] = new TH2D(Syst+"_FracCovariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";i bin "+XTitle+";j bin "+XTitle,NBins,ArrayBins,NBins,ArrayBins);
				Covariances[WhichRun][WhichPlot] = new TH2D(Syst+"_Covariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";i bin "+XTitle+";j bin "+XTitle,NBins,ArrayBins,NBins,ArrayBins);

			}

			for (int WhichXBin = 0; WhichXBin < NBins; WhichXBin++) { 

				for (int WhichYBin = 0; WhichYBin < NBins; WhichYBin++) {

					FracCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,0.);
					FracCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,0.);

					Covariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,0.);
					Covariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,0.);

				}  // end of the loop over bin Y

			} // end of the loop over bin X

			// -------------------------------------------------------------------------------------------------------

			for (int WhichXBin = 0; WhichXBin < NBins; WhichXBin++) { 

				for (int WhichYBin = 0; WhichYBin < NBins; WhichYBin++) {

					// X Bin entry / error

					double DataEntryX = DataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
					double DataErrorX = DataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

					// Y Bin entry / error

					double DataEntryY = DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
					double DataErrorY = DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

					// -------------------------------------------------------------------------------------------------------
					// -------------------------------------------------------------------------------------------------------

					double CovFracEntry = -99.;
					double CovFracError = -99.;

					double CovEntry = -99.;
					double CovError = -99.;

					double DataEntryXCV = DataEntryX;
					double DataEntryYCV = DataEntryY;

					if (WhichXBin == WhichYBin) {

						DataEntryX = DataErrorX;
						DataEntryY = DataErrorY;

					} else {

						DataEntryX = 0.;
						DataEntryY = 0.;

					} 

					CovFracEntry = TMath::Max( (DataEntryX / DataEntryXCV) * (DataEntryY / DataEntryYCV) ,1E-8);
					CovEntry = TMath::Max( DataEntryX * DataEntryY,1E-8);

					CovFracError = 1E-8;
					CovError = 1E-8;					

					// -------------------------------------------------------------------------------------------------------

					// Setting the elements of the Cov Matrix

					FracCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovFracEntry);
					FracCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovFracError);

					Covariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovEntry);
					Covariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovError);

					// -------------------------------------------------------------------------------------------------------

				}  // end of the loop over bin Y

			} // end of the loop over bin X
			
			FileCovarianceMatrices->cd();
			FracCovariances[WhichRun][WhichPlot]->Write();
			Covariances[WhichRun][WhichPlot]->Write();	

		} // End of the loop over the plots

		// ---------------------------------------------------------------------------------------	

		// Close covariance matrix file & base files

		FileCovarianceMatrices->Close();

		MCFileSample[WhichRun]->Close();
		BeamOnFileSample[WhichRun]->Close();
		BeamOffFileSample[WhichRun]->Close();
		DirtFileSample[WhichRun]->Close();

		// ---------------------------------------------------------------------------------------	

		cout << endl << "File with "+Syst+" Covariance matrices " << FileName << " has been created" << endl << endl;

		// ---------------------------------------------------------------------------------------	

	} // End of the loop over the runs	

} // End of the program