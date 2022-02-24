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

// --------------------------------------------------------------------------------------------------------------------------------------------

void StoreCanvas(TH2D* h, TString Label, TString Syst, TString PlotNames, TString BaseMC, TString Runs, TString Tune) {

	TString CanvasName = Syst+"_"+PlotNames+BaseMC+"_"+Runs;
	TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
	PlotCanvas->cd();
	PlotCanvas->SetBottomMargin(0.16);
	PlotCanvas->SetLeftMargin(0.15);
	PlotCanvas->SetRightMargin(0.25);			
	
	gStyle->SetMarkerSize(1.5);
	gStyle->SetPaintTextFormat("4.3f");			
	
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleSize(TextSize);
	h->GetXaxis()->SetLabelSize(TextSize);			
	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetNdivisions(8);
	
	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetTitleSize(TextSize);
	h->GetYaxis()->SetLabelSize(TextSize);			
	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetNdivisions(5);
	h->GetYaxis()->SetTitleOffset(1.);			

	h->SetTitle(Runs + " " + Syst);	

	double CovMax = TMath::Min(1.,1.05 * h->GetMaximum());
	double CovMin = TMath::Min(0.,1.05 * h->GetMinimum());
	h->GetZaxis()->SetRangeUser(CovMin,CovMax);
	h->GetZaxis()->CenterTitle();
	h->GetZaxis()->SetTitleFont(FontStyle);
	h->GetZaxis()->SetTitleSize(TextSize);
	h->GetZaxis()->SetLabelFont(FontStyle);
	h->GetZaxis()->SetLabelSize(TextSize-0.01);
	h->GetZaxis()->SetNdivisions(5);

	h->SetMarkerColor(kWhite);			
	h->SetMarkerSize(1.5);
	//h->Draw("text colz e"); 
	h->Draw("colz");
	
	PlotCanvas->SaveAs(PlotPath+BaseMC+"/"+Tune+"WienerSVD_"+Syst+"_"+Label+"CovarianceMatrices_"+PlotNames+BaseMC+"_"+Runs+"_"+UBCodeVersion+".pdf");
	
	delete PlotCanvas;	

}

// -----------------------------------------------------------------------------------------------------------------------------------

TH1D* Multiply(TH1D* True, TH2D* SmearMatrix) {

	TH1D* TrueClone = (TH1D*)(True->Clone());

	int XBins = SmearMatrix->GetXaxis()->GetNbins();
	int YBins = SmearMatrix->GetYaxis()->GetNbins();

	if (XBins != YBins) { std::cout << "Not symmetric matrix" << std::endl; }

	TVectorD signal(XBins);
	TMatrixD response(XBins,XBins);

	H2V(True, signal);
	H2M(SmearMatrix, response, kFALSE); // X axis: Reco, Y axis: True

	TVectorD RecoSpace = response * signal;
	V2H(RecoSpace, TrueClone);	

	return TrueClone;

}

// ------------------------------------------------------------------------------------------------------------------------------

// TString Syst = "Stat" "POT" "NTarget" "LY" "TPC" "SCERecomb2" "XSec" "G4" "Flux" "Dirt" "MC_Stat"

void WienerSVD_CovarianceMatrices(TString Syst = "None",TString BaseMC = "Overlay9",TString BeamOnSample = "BeamOn9",TString BeamOffSample = "ExtBNB9",TString DirtSample = "OverlayDirt9", TString Tune = "") {

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

	const int N1DPlots = PlotNames.size();
		
	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");			
	Runs.push_back("Combined");

	if (Syst == "LY" || Syst == "TPC" || Syst == "SCERecomb2" || Syst == "MC_LY" || Syst == "MC_TPC" || Syst == "MC_SCERecomb2" 
	|| Syst == "SmEff_LY" || Syst == "SmEff_TPC" || Syst == "SmEff_SCERecomb2" ) {

		Runs.clear();
		Runs.push_back("Run3");

	}

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

	vector <TH1D*> TrueCC1pPlots; TrueCC1pPlots.resize(N1DPlots); 
	vector <TH2D*> CC1pResponseMatrix; CC1pResponseMatrix.resize(N1DPlots); 
	vector <TH1D*> ClosureTestCC1pPlots; ClosureTestCC1pPlots.resize(N1DPlots); 
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

	// -------------------------------------------------------------------------------------------------------------------------------------

	// Alternative MC Models

	std::vector<TString> AltModels;
	std::vector<int> AltUniverses;
	std::vector<int> Colors;
	std::vector<int> Universes;
	std::vector<TString> UniAltModels;
	std::vector<TH1D*> UniHistoFlux;
	std::vector<double> UniFlux;

	//----------------------------------------//

	if (Syst == "LY" || Syst == "MC_LY" || Syst == "SmEff_LY") {

		AltModels.push_back("_LYDown"); Colors.push_back(kRed+1); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_LYRayleigh"); Colors.push_back(kGreen+2); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_LYAttenuation"); Colors.push_back(kOrange+1); Universes.push_back(1); AltUniverses.push_back(1);

	}

	//----------------------------------------//	

	if (Syst == "TPC" || Syst == "MC_TPC" || Syst == "SmEff_TPC") {

		AltModels.push_back("_X"); Colors.push_back(kRed+1); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_YZ"); Colors.push_back(kGreen+2); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_ThetaXZ"); Colors.push_back(kOrange+1); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_ThetaYZ"); Colors.push_back(kBlue-3); Universes.push_back(1); AltUniverses.push_back(1);

	}

	//----------------------------------------//	

	if (Syst == "SCERecomb2" || Syst == "MC_SCERecomb2" || Syst == "SmEff_SCERecomb2") {

		AltModels.push_back("_SCE"); Colors.push_back(kBlue); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_Recombination2"); Colors.push_back(kMagenta); Universes.push_back(1); AltUniverses.push_back(1);

	}

	//----------------------------------------//

	if (Syst == "XSec" || Syst == "MC_XSec" || Syst == "SmEff_XSec") {

		UniAltModels.push_back("_AxFFCCQEshape_UBGenie"); Universes.push_back(2); 
		UniAltModels.push_back("_DecayAngMEC_UBGenie"); Universes.push_back(2);
		UniAltModels.push_back("_NormCCCOH_UBGenie"); Universes.push_back(2);
		UniAltModels.push_back("_NormNCCOH_UBGenie"); Universes.push_back(2);
		UniAltModels.push_back("_RPA_CCQE_UBGenie"); Universes.push_back(2);
		UniAltModels.push_back("_ThetaDelta2NRad_UBGenie"); Universes.push_back(2);
		UniAltModels.push_back("_Theta_Delta2Npi_UBGenie"); Universes.push_back(2);
		UniAltModels.push_back("_VecFFCCQEshape_UBGenie"); Universes.push_back(2);
		UniAltModels.push_back("_XSecShape_CCMEC_UBGenie"); Universes.push_back(2);
		UniAltModels.push_back("_All_UBGenie"); Universes.push_back(100);

		for (int UniAlt = 0; UniAlt < (int)(UniAltModels.size()); UniAlt++ ) {

			for (int Uni = 0; Uni < Universes[UniAlt]; Uni++ ) {

				AltModels.push_back(UniAltModels[UniAlt]+"_"+TString(std::to_string(Uni))); Colors.push_back(kGreen+2);
	 			AltUniverses.push_back(Universes[UniAlt]);

			}

		}

	}

	//----------------------------------------//	

	if (Syst == "G4" || Syst == "MC_G4" || Syst == "SmEff_G4") {

		UniAltModels.push_back("_reinteractions"); Universes.push_back(100);

		for (int UniAlt = 0; UniAlt < (int)(UniAltModels.size()); UniAlt++ ) {

			for (int Uni = 0; Uni < Universes[UniAlt]; Uni++ ) {

				AltModels.push_back(UniAltModels[UniAlt]+"_"+TString(std::to_string(Uni))); Colors.push_back(kGreen+2);
				AltUniverses.push_back(Universes[UniAlt]);
			
			}

		}

	}

	//----------------------------------------//	

	if (Syst == "Flux" || Syst == "MC_Flux" || Syst == "SmEff_Flux") {

		UniAltModels.push_back("_fluxes"); Universes.push_back(100);

		int NUniAltModels = (int)(UniAltModels.size());

		for (int UniAlt = 0; UniAlt < NUniAltModels; UniAlt++ ) {

			for (int Uni = 0; Uni < Universes[UniAlt]; Uni++ ) {

				AltModels.push_back(UniAltModels[UniAlt]+"_"+TString(std::to_string(Uni))); Colors.push_back(kGreen+2);
				AltUniverses.push_back(Universes[UniAlt]);

				TString FluxHistoName = "numu_ms_total/hEnumu_ms_"+TString(std::to_string(Uni));
				UniHistoFlux.push_back( (TH1D*)(FluxFile->Get(FluxHistoName) ) );
	
			}

		}

	}

	//----------------------------------------//	

	if (Syst == "MC_Stat") {

		UniAltModels.push_back("_MC_Stat"); Universes.push_back(100);

		for (int UniAlt = 0; UniAlt < (int)(UniAltModels.size()); UniAlt++ ) {

			for (int Uni = 0; Uni < Universes[UniAlt]; Uni++ ) {

				AltModels.push_back(UniAltModels[UniAlt]+"_"+TString(std::to_string(Uni))); Colors.push_back(kGreen+2);
				AltUniverses.push_back(Universes[UniAlt]);
			
			}

		}

	}	

	//----------------------------------------//

	int NAltModels = AltModels.size();

	vector< vector<TFile*> > AltMCFileSample; AltMCFileSample.resize(NRuns,vector<TFile*>(NAltModels));
	vector< vector<TFile*> > AltMCFileSampleResponseMatrix; AltMCFileSampleResponseMatrix.resize(NRuns,vector<TFile*>(NAltModels));

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	// Alternative Plots

	vector <vector <TH1D*> > AltForwardFoldedCC1pPlots; AltForwardFoldedCC1pPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));
	vector <vector <TH1D*> > AltCC1pPlots; AltCC1pPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));
	vector <vector <TH1D*> > AltNonCC1pPlots; AltNonCC1pPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));
	vector <vector <TH1D*> > AltBeamOnPlots; AltBeamOnPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));
	vector <vector <TH1D*> > AltBeamOffPlots; AltBeamOffPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------

//		// Until Run3 NuWro is produced
//		if ( Runs[WhichRun] == "Run3" && (BaseMC == "Overlay9NuWro" || BeamOnSample == "Overlay9NuWro") ) 
//			{ continue;}

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------
	
		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);						
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface);

		if (Syst == "Flux" || Syst == "MC_Flux" || Syst == "SmEff_Flux") {

			UniFlux.clear();

			for (int UniAlt = 0; UniAlt < NAltModels; UniAlt++ ) {

				double UniFluxPOT = (UniHistoFlux[UniAlt]->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface);
				UniFlux.push_back(UniFluxPOT);

			}

		}

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------

		TString FileName = MigrationMatrixPath+Tune+"WienerSVD_"+Syst+"_CovarianceMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		if (BeamOnSample != "BeamOn9") { FileName = MigrationMatrixPath+BeamOnSample+"WienerSVD_"+Syst+"_CovarianceMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }

		TFile* FileCovarianceMatrices = new TFile(FileName,"recreate");
		
		// Open base files

		TString ExactFileLocation = PathToFiles+CutExtension;
		TString TStringBaseMC = ExactFileLocation+"/"+Tune+"STVStudies_"+BaseMC+"_"+Runs[WhichRun]+CutExtension+".root";
		TString TrueTStringBaseMC = PathToFiles+"/"+Tune+"TruthSTVAnalysis_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TString ResponseFileName = MigrationMatrixPath+Tune+"FileResponseMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";

		if (Syst == "LY" || Syst == "TPC" || Syst == "MC_LY" || Syst == "MC_TPC"  || Syst == "SmEff_LY" || Syst == "SmEff_TPC") { 

			TStringBaseMC = ExactFileLocation+"/"+Tune+"STVStudies_"+BaseMC+"_"+Runs[WhichRun]+"_CV"+CutExtension+".root"; 
			TrueTStringBaseMC = PathToFiles+"/"+Tune+"TruthSTVAnalysis_"+BaseMC+"_"+Runs[WhichRun]+"_CV_"+UBCodeVersion+".root"; 
			ResponseFileName = MigrationMatrixPath+Tune+"FileResponseMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_CV_"+UBCodeVersion+".root";

		}

		if (Syst == "SCERecomb2" || Syst == "MC_SCERecomb2"  || Syst == "SmEff_SCERecomb2") { 

			TStringBaseMC = ExactFileLocation+"/"+Tune+"STVStudies_"+BaseMC+"_"+Runs[WhichRun]+"_CVextra"+CutExtension+".root"; 
			TrueTStringBaseMC = PathToFiles+"/"+Tune+"TruthSTVAnalysis_"+BaseMC+"_"+Runs[WhichRun]+"_CVextra_"+UBCodeVersion+".root"; 
			ResponseFileName = MigrationMatrixPath+Tune+"FileResponseMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_CVextra_"+UBCodeVersion+".root";

		}

		TFile* FileResponseMatrices = new TFile(ResponseFileName,"readonly");
		TrueMCFileSample[WhichRun] = TFile::Open(TrueTStringBaseMC,"readonly");
		MCFileSample[WhichRun] = TFile::Open(TStringBaseMC,"readonly");
		BeamOnFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOnSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");
		BeamOffFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOffSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");
		DirtFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/"+Tune+"STVStudies_"+DirtSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");

		// -------------------------------------------------------------------------------------

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

			// -------------------------------------------------------------------------------------------------------

			// Grab base plots

			CC1pResponseMatrix[WhichPlot] = (TH2D*)(FileResponseMatrices->Get("POTScaledCC1pReco"+PlotNames[WhichPlot]+"2D"));
			TrueCC1pPlots[WhichPlot] = (TH1D*)(TrueMCFileSample[WhichRun]->Get("True"+PlotNames[WhichPlot]));
			CC1pPlots[WhichPlot] = (TH1D*)(MCFileSample[WhichRun]->Get("CC1pReco"+PlotNames[WhichPlot]));
			NonCC1pPlots[WhichPlot] = (TH1D*)(MCFileSample[WhichRun]->Get("NonCC1pReco"+PlotNames[WhichPlot]));
			BeamOnPlots[WhichPlot] = (TH1D*)(BeamOnFileSample[WhichRun]->Get("Reco"+PlotNames[WhichPlot]));
			BeamOffPlots[WhichPlot] = (TH1D*)(BeamOffFileSample[WhichRun]->Get("Reco"+PlotNames[WhichPlot]));
			DirtPlots[WhichPlot] = (TH1D*)(DirtFileSample[WhichRun]->Get("Reco"+PlotNames[WhichPlot]));
			ClosureTestCC1pPlots[WhichPlot] = Multiply(TrueCC1pPlots[WhichPlot],CC1pResponseMatrix[WhichPlot]);

			// -------------------------------------------------------------------------------------------------------

			// Grab alternative plots

			for (int alt = 0; alt < NAltModels; alt++ ) {

				if ( (Syst == "LY" || Syst == "MC_LY" || Syst == "SmEff_LY") && Runs[WhichRun] == "Run1" && AltModels[alt] == "_LYAttenuation") 
					{ continue;}

				// Open Alternative MC files & ReAlternative Response Matrices

				TString TStringAltBaseMC = ExactFileLocation+"/"+Tune+"STVStudies_"+BaseMC+"_"+Runs[WhichRun]+AltModels[alt]+CutExtension+".root";
				AltMCFileSample[WhichRun][alt] = TFile::Open(TStringAltBaseMC,"readonly");

				AltCC1pPlots[WhichPlot][alt] = (TH1D*)(AltMCFileSample[WhichRun][alt]->Get("CC1pReco"+PlotNames[WhichPlot]));
				AltNonCC1pPlots[WhichPlot][alt] = (TH1D*)(AltMCFileSample[WhichRun][alt]->Get("NonCC1pReco"+PlotNames[WhichPlot]));

				TString TStringAltBaseMCResponseMatrix = MigrationMatrixPath+Tune+"FileResponseMatrices_"+BaseMC+"_"+Runs[WhichRun]+AltModels[alt]+"_"+UBCodeVersion+".root";
				AltMCFileSampleResponseMatrix[WhichRun][alt] = TFile::Open(TStringAltBaseMCResponseMatrix,"readonly");

				TH2D* AltResponseMatrix = (TH2D*)(AltMCFileSampleResponseMatrix[WhichRun][alt]->Get("POTScaledCC1pReco"+PlotNames[WhichPlot]+"2D"));
				AltForwardFoldedCC1pPlots[WhichPlot][alt] = Multiply(TrueCC1pPlots[WhichPlot],AltResponseMatrix);

				AltCC1pPlots[WhichPlot][alt]->SetDirectory(0); // to decouple it from the open file directory
				AltNonCC1pPlots[WhichPlot][alt]->SetDirectory(0); // to decouple it from the open file directory

				AltMCFileSample[WhichRun][alt]->Close();

				AltForwardFoldedCC1pPlots[WhichPlot][alt]->SetDirectory(0); // to decouple it from the open file directory

				AltMCFileSampleResponseMatrix[WhichRun][alt]->Close();

			}

			// -------------------------------------------------------------------------------------------------------

			// Given that we want to construct fractional covariance matrices
			// We have to be very careful as to what we are normalizing to
			
			// For the statistical uncertainties, normalize to BeamOn - ExtBNB - Dirt
			// For the systematic uncertainties, normalize to CV (MC signal + MC bkg)

			TH1D* DataPlot = (TH1D*)CC1pPlots[WhichPlot]->Clone();
			DataPlot->Add(NonCC1pPlots[WhichPlot]);

			TH1D* AltDataPlot = (TH1D*)CC1pPlots[WhichPlot]->Clone();
			AltDataPlot->Add(NonCC1pPlots[WhichPlot]);

			if ( string(Syst).find("Stat") != std::string::npos ) {

				if (Syst == "Stat") {

					// Don't forget to subtract the cosmic part, the dirt and the NonCC1p bkgs
					// The MC uncertainties have already been taken care of with the MC event rate covariances

					DataPlot = (TH1D*)BeamOnPlots[WhichPlot]->Clone();
					AltDataPlot = (TH1D*)BeamOnPlots[WhichPlot]->Clone();

					DataPlot->Add(BeamOffPlots[WhichPlot],-1.);
					DataPlot->Add(DirtPlots[WhichPlot],-1.);
					DataPlot->Add(NonCC1pPlots[WhichPlot],-1.);

				}

				if (Syst == "MC_Stat") {

					DataPlot = (TH1D*)CC1pPlots[WhichPlot]->Clone();
					AltDataPlot = (TH1D*)CC1pPlots[WhichPlot]->Clone();

				}

			}


			if (BeamOnSample == "BeamOn9" || BeamOnSample == "Overlay9") { // Beam On Sample, need to subtract ALL backgrounds

				for (int alt = 0; alt < NAltModels; alt++ ) {

					if ( (Syst == "LY" || Syst == "MC_LY"  || Syst == "SmEff_LY") && Runs[WhichRun] == "Run1" && AltModels[alt] == "_LYAttenuation") 
						{ continue;}

					// Make a copy of the alternative "CC1p" plot
					// Which was derived via the multiplication of the 
					// CV truth plot times the corresponding response matrix for a given universe 
					// and then add the beam related backgrounds in the given universe

					AltBeamOnPlots[WhichPlot][alt] = (TH1D*)AltForwardFoldedCC1pPlots[WhichPlot][alt]->Clone();
					AltBeamOnPlots[WhichPlot][alt]->Add(AltNonCC1pPlots[WhichPlot][alt]);
					//AltBeamOnPlots[WhichPlot][alt]->Add(BeamOffPlots[WhichPlot]);
					//AltBeamOnPlots[WhichPlot][alt]->Add(DirtPlots[WhichPlot]);

					// // Now subtract the CV bkgs
					//AltBeamOnPlots[WhichPlot][alt]->Add(NonCC1pPlots[WhichPlot],-1.);
					//AltBeamOnPlots[WhichPlot][alt]->Add(BeamOffPlots[WhichPlot],-1.);
					//AltBeamOnPlots[WhichPlot][alt]->Add(DirtPlots[WhichPlot],-1.);

					if (Syst == "MC_Stat") {

						AltBeamOnPlots[WhichPlot][alt] = (TH1D*)AltForwardFoldedCC1pPlots[WhichPlot][alt]->Clone();

					}					

				}

				if (Syst == "Dirt" || Syst == "MC_Dirt" || Syst == "SmEff_Dirt") {

					// Start from the CV reco plot (MC signal + MC bkg)
					DataPlot = (TH1D*)BeamOnPlots[WhichPlot]->Clone();
					// Add the CV dirt sample
					DataPlot->Add(DirtPlots[WhichPlot]);

					// Start from the modified reco plot (MC signal + MC bkg)
					AltDataPlot = (TH1D*)BeamOnPlots[WhichPlot]->Clone();
					// Scale the dirt sample down by 25%
					TH1D* DirtClone = (TH1D*)(DirtPlots[WhichPlot]->Clone());
					DirtClone->Scale(0.75);
					// Add the alternative dirt sample
					AltDataPlot->Add(DirtClone);

					//AltDataPlot->Add(NonCC1pPlots[WhichPlot]);
					//AltDataPlot->Add(BeamOffPlots[WhichPlot]);

					// // Now subtract the CV bkgs
					//AltDataPlot->Add(NonCC1pPlots[WhichPlot],-1.);
					//AltDataPlot->Add(BeamOffPlots[WhichPlot],-1.);
					//AltDataPlot->Add(DirtPlots[WhichPlot],-1.);

				}				

			} 

			// -------------------------------------------------------------------------------------------------------
			// -------------------------------------------------------------------------------------------------------

			double BeamOnMax = 1.25*DataPlot->GetMaximum();

			// GENIE CV LY / TPC / SCERecomb2 overlays

			if (	BaseMC == "Overlay9" 
				&& (Syst == "LY" || Syst == "TPC" || Syst == "SCERecomb2" 
				|| Syst == "SmEff_LY" || Syst == "SmEff_TPC" || Syst == "SmEff_SCERecomb2" 
				|| Syst == "MC_LY" || Syst == "MC_TPC" || Syst == "MC_SCERecomb2") 
			) {

				TString EventRatePlotCanvasName = Runs[WhichRun]+"_"+PlotNames[WhichPlot]+"_"+Syst;
				TCanvas* EventRatePlotCanvas = new TCanvas(EventRatePlotCanvasName,EventRatePlotCanvasName,205,34,1024,768);
				EventRatePlotCanvas->SetBottomMargin(0.16);
				EventRatePlotCanvas->SetLeftMargin(0.15);

				TLegend* leg = new TLegend(0.15,0.91,0.9,0.99);
				leg->SetBorderSize(0);
				leg->SetTextFont(FontStyle);
				leg->SetTextSize(TextSize-0.02);
				//if (Syst == "TPC" || Syst == "MC_TPC" || Syst == "SmEff_TPC") { leg->SetTextSize(TextSize-0.03); }
				leg->SetNColumns(3);

				DataPlot->GetXaxis()->CenterTitle();
				DataPlot->GetXaxis()->SetTitleFont(FontStyle);
				DataPlot->GetXaxis()->SetTitleSize(TextSize);
				DataPlot->GetXaxis()->SetLabelFont(FontStyle);
				DataPlot->GetXaxis()->SetLabelSize(TextSize);
				DataPlot->GetXaxis()->SetNdivisions(6);

				DataPlot->GetYaxis()->SetTitle("Flux Ave Events / "+ToString(DataPOT));
				DataPlot->GetYaxis()->CenterTitle();
				DataPlot->GetYaxis()->SetTitleFont(FontStyle);
				DataPlot->GetYaxis()->SetTitleSize(TextSize);
				DataPlot->GetYaxis()->SetLabelFont(FontStyle);
				DataPlot->GetYaxis()->SetLabelSize(TextSize);
				DataPlot->GetYaxis()->SetNdivisions(6);
				DataPlot->GetYaxis()->SetRangeUser(0.,BeamOnMax);
				DataPlot->GetYaxis()->SetTitleOffset(0.95);

				DataPlot->SetMarkerColor(kBlack);
				DataPlot->SetMarkerStyle(20);
				DataPlot->SetMarkerSize(2.);
				DataPlot->SetLineColor(kBlack);
				DataPlot->SetLineWidth(3);
				DataPlot->Draw("p0 hist same");

				leg->AddEntry(DataPlot,"CV","p");

				EventRatePlotCanvas->cd();

				for (int alt = 0; alt < NAltModels; alt++ ) {

					if ( (Syst == "LY" || Syst == "MC_LY"  || Syst == "SmEff_LY") && Runs[WhichRun] == "Run1" && AltModels[alt] == "_LYAttenuation") 
						{ continue;}

					AltBeamOnPlots[WhichPlot][alt]->SetMarkerColor(Colors[alt]);
					AltBeamOnPlots[WhichPlot][alt]->SetMarkerStyle(20);
					AltBeamOnPlots[WhichPlot][alt]->SetMarkerSize(2.);
					AltBeamOnPlots[WhichPlot][alt]->SetLineColor(Colors[alt]);
					AltBeamOnPlots[WhichPlot][alt]->SetLineWidth(3);
					AltBeamOnPlots[WhichPlot][alt]->Draw("p0 hist same");

					TString LabelInit = AltModels[alt];
					TString Label = LabelInit.ReplaceAll("_","");
					leg->AddEntry(AltBeamOnPlots[WhichPlot][alt],Label,"p");

				}

				DataPlot->Draw("p0 hist same");

				// // Closure Test 
				//ClosureTestCC1pPlots[WhichPlot]->SetMarkerColor(kMagenta);
				//ClosureTestCC1pPlots[WhichPlot]->Draw("p0 hist same");

				leg->Draw();

				TString EventRateCanvasName = "/"+Tune+"EventRate_WienerSVD_"+Syst+"_CovarianceMatrices_"+PlotNames[WhichPlot]+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf";
				EventRatePlotCanvas->SaveAs(PlotPath+BaseMC+EventRateCanvasName);
				delete EventRatePlotCanvas;

			} // End of the cases where we overlay the event rates for GENIE CV  LY / TPC / SCERecomb systematics

			// -------------------------------------------------------------------------------------------------------

			// GENIE CV multisim overlays

			if (BaseMC == "Overlay9" && (Syst == "XSec" || Syst == "G4" || Syst == "Flux" || Syst == "SmEff_XSec" || Syst == "SmEff_G4" || Syst == "SmEff_Flux" || Syst == "MC_Stat") ) {

				for (int unialt = 0; unialt < (int)(UniAltModels.size()); unialt++ ) {

					TString EventRatePlotCanvasName = Runs[WhichRun]+"_"+PlotNames[WhichPlot]+"_"+Syst;
					TCanvas* EventRatePlotCanvas = new TCanvas(EventRatePlotCanvasName,EventRatePlotCanvasName,205,34,1024,768);
					EventRatePlotCanvas->SetBottomMargin(0.16);
					EventRatePlotCanvas->SetLeftMargin(0.15);

					TLegend* leg = new TLegend(0.32,0.91,0.9,0.99);
					leg->SetBorderSize(0);
					leg->SetTextFont(FontStyle);
					leg->SetTextSize(TextSize-0.02);
					leg->SetNColumns(3);

					DataPlot->GetXaxis()->CenterTitle();
					DataPlot->GetXaxis()->SetTitleFont(FontStyle);
					DataPlot->GetXaxis()->SetTitleSize(TextSize);
					DataPlot->GetXaxis()->SetLabelFont(FontStyle);
					DataPlot->GetXaxis()->SetLabelSize(TextSize);
					DataPlot->GetXaxis()->SetNdivisions(6);

					DataPlot->GetYaxis()->SetTitle("Flux Ave Events / "+ToString(DataPOT));
					DataPlot->GetYaxis()->CenterTitle();
					DataPlot->GetYaxis()->SetTitleFont(FontStyle);
					DataPlot->GetYaxis()->SetTitleSize(TextSize);
					DataPlot->GetYaxis()->SetLabelFont(FontStyle);
					DataPlot->GetYaxis()->SetLabelSize(TextSize);
					DataPlot->GetYaxis()->SetNdivisions(6);
					DataPlot->GetYaxis()->SetRangeUser(0.,BeamOnMax);
					DataPlot->GetYaxis()->SetTitleOffset(0.95);

					DataPlot->SetMarkerColor(kBlack);
					DataPlot->SetMarkerStyle(20);
					DataPlot->SetMarkerSize(2.);
					DataPlot->SetLineColor(kBlack);
					DataPlot->SetLineWidth(3);
					EventRatePlotCanvas->cd();
					DataPlot->Draw("p0 hist same");

					leg->AddEntry(DataPlot,"CV","p");

					int FirstEntry = 0;

					for (int alt = 0; alt < NAltModels; alt++ ) {

						if (string(AltModels[alt]).find(UniAltModels[unialt]) != std::string::npos) {

							AltBeamOnPlots[WhichPlot][alt]->SetMarkerColor(Colors[alt]);
							AltBeamOnPlots[WhichPlot][alt]->SetMarkerStyle(20);
							AltBeamOnPlots[WhichPlot][alt]->SetMarkerSize(2.);
							AltBeamOnPlots[WhichPlot][alt]->SetLineColor(Colors[alt]);
							AltBeamOnPlots[WhichPlot][alt]->SetLineWidth(3);
							AltBeamOnPlots[WhichPlot][alt]->Draw("p0 hist same");

							TString LabelInit = AltModels[alt];
							TString Label = LabelInit.ReplaceAll("_","").ReplaceAll("0","");
							if (FirstEntry == 0) { leg->AddEntry(AltBeamOnPlots[WhichPlot][alt],Label,"p"); }

							FirstEntry++;

						} else { continue; }

					} // End of the loop over all the universes

					DataPlot->Draw("p0 hist same");

					// Closure Test 
					ClosureTestCC1pPlots[WhichPlot]->SetMarkerColor(kMagenta);
					//ClosureTestCC1pPlots[WhichPlot]->Draw("p0 hist same");

					leg->Draw();

					TString EventRateCanvasName = "/"+Tune+"EventRate_WienerSVD_"+Syst+UniAltModels[unialt]+"_CovarianceMatrices_"+PlotNames[WhichPlot]+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf";
					EventRatePlotCanvas->SaveAs(PlotPath+BaseMC+EventRateCanvasName);
					delete EventRatePlotCanvas;

				} // End of the loop over a specific universe name

			} // End of the cases where we overlay the event rates for GENIE CV  XSec/G4/Flux systematics	

			// -------------------------------------------------------------------------------------------------------
			// -------------------------------------------------------------------------------------------------------
			
			int NBins = BeamOnPlots[WhichPlot]->GetXaxis()->GetNbins();
			const double* ArrayBins = BeamOnPlots[WhichPlot]->GetXaxis()->GetXbins()->GetArray();
			TString XTitle = BeamOnPlots[WhichPlot]->GetXaxis()->GetTitle();

			// -------------------------------------------------------------------------------------------------------

			// Declare the matrix & initialize the entries to 0

			if (PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot" || NBins == 1) {

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

					double AltDataEntryX = DataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
					double AltDataErrorX = DataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

					// Y Bin entry / error

					double DataEntryY = DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
					double DataErrorY = DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

					double AltDataEntryY = DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
					double AltDataErrorY = DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

					// -------------------------------------------------------------------------------------------------------
					// -------------------------------------------------------------------------------------------------------

					double CovFracEntry = -99.;
					double CovFracError = -99.;

					double CovEntry = -99.;
					double CovError = -99.;

					// Based on the type of systematic unceratainty, choose how to handle it

					if (Syst == "NTarget" || Syst == "MC_NTarget"  || Syst == "SmEff_NTarget") {

						AltDataEntryX = (1+NTargetUncertainty) * DataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorX = (1+NTargetUncertainty) * DataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

						AltDataEntryY = (1+NTargetUncertainty) * DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorY = (1+NTargetUncertainty) * DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

						CovFracEntry = TMath::Max( ( (AltDataEntryX - DataEntryX) / DataEntryX) * ( (AltDataEntryY - DataEntryY) / DataEntryY),1E-8);
						CovEntry = TMath::Max( (AltDataEntryX - DataEntryX) * (AltDataEntryY - DataEntryY),1E-8);

						// CovError = TMath::Max( TMath::Sqrt( 
						// 	TMath::Power(DataEntryY - AltDataEntryY,2.) * ( TMath::Power(DataErrorX,2.) + TMath::Power(AltDataErrorX,2.) ) +
						// 	TMath::Power(DataEntryX - AltDataEntryX,2.) * ( TMath::Power(DataErrorY,2.) + TMath::Power(AltDataErrorY,2.) ) ), 1E-10) ;

						CovFracError = 1E-8;
						CovError = 1E-8;

					} else if (Syst == "POT" || Syst == "MC_POT" || Syst == "SmEff_POT") {

						AltDataEntryX = (1+POTUncertainty) * DataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorX = (1+POTUncertainty) * DataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

						AltDataEntryY = (1+POTUncertainty) * DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorY = (1+POTUncertainty) * DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

						CovFracEntry = TMath::Max( ((AltDataEntryX - DataEntryX) / DataEntryX) * ( (AltDataEntryY - DataEntryY) / DataEntryY ),1E-8);
						CovEntry = TMath::Max( (AltDataEntryX - DataEntryX) * (AltDataEntryY - DataEntryY),1E-8);

//						CovError = TMath::Max( TMath::Sqrt( 
//							TMath::Power(DataEntryY - AltDataEntryY,2.) * ( TMath::Power(DataErrorX,2.) + TMath::Power(AltDataErrorX,2.) ) +
//							TMath::Power(DataEntryX - AltDataEntryX,2.) * ( TMath::Power(DataErrorY,2.) + TMath::Power(AltDataErrorY,2.) ) ), 1E-10) ;

						CovFracError = 1E-8;
						CovError = 1E-8;


					} else if (Syst == "Dirt" || Syst == "MC_Dirt" || Syst == "SmEff_Dirt") {

						AltDataEntryX = AltDataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorX = AltDataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

						AltDataEntryY = AltDataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorY = AltDataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

						CovFracEntry = TMath::Max( ((AltDataEntryX - DataEntryX) / DataEntryX ) * ( (AltDataEntryY - DataEntryY) / DataEntryY ),1E-8);
						CovEntry = TMath::Max( (AltDataEntryX - DataEntryX) * (AltDataEntryY - DataEntryY),1E-8);

//						CovError = TMath::Max( TMath::Sqrt( 
//							TMath::Power(DataEntryY - AltDataEntryY,2.) * ( TMath::Power(DataErrorX,2.) + TMath::Power(AltDataErrorX,2.) ) +
//							TMath::Power(DataEntryX - AltDataEntryX,2.) * ( TMath::Power(DataErrorY,2.) + TMath::Power(AltDataErrorY,2.) ) ), 1E-10) ;

						CovFracError = 1E-8;
						CovError = 1E-8;

					} else if (Syst == "Stat" || Syst == "SmEff_Stat") {

						double DataEntryXCV = DataEntryX;
						double DataEntryYCV = DataEntryY;

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

						CovFracEntry = TMath::Max( (DataEntryX / DataEntryXCV) * (DataEntryY / DataEntryYCV) ,1E-8);
						CovEntry = TMath::Max( DataEntryX * DataEntryY,1E-8);

//						CovFracEntry = TMath::Max( ((AltDataEntryX - DataEntryX) / DataEntryXCV) * ( (AltDataEntryY - DataEntryY) / DataEntryYCV),1E-8);
//						CovEntry = TMath::Max( (AltDataEntryX - DataEntryX) * (AltDataEntryY - DataEntryY),1E-8);

//						CovError = TMath::Max( TMath::Sqrt( 
//							TMath::Power(DataEntryY - AltDataEntryY,2.) * ( TMath::Power(DataErrorX,2.) + TMath::Power(AltDataErrorX,2.) ) +
//							TMath::Power(DataEntryX - AltDataEntryX,2.) * ( TMath::Power(DataErrorY,2.) + TMath::Power(AltDataErrorY,2.) ) ), 1E-10) ;

						CovFracError = 1E-8;
						CovError = 1E-8;

					} else if (
						Syst == "LY" || Syst == "TPC" || Syst == "SCERecomb2" || Syst == "XSec" || Syst == "G4" || Syst == "Flux" ||
						Syst == "MC_LY" || Syst == "MC_TPC" || Syst == "MC_SCERecomb2" || Syst == "MC_XSec" || Syst == "MC_G4" || Syst == "MC_Flux" ||
						Syst == "SmEff_LY" || Syst == "SmEff_TPC" || Syst == "SmEff_SCERecomb2" || Syst == "SmEff_XSec" || Syst == "SmEff_G4" || Syst == "SmEff_Flux" ||
						Syst == "MC_Stat"
					) {						

						for (int alt = 0; alt < NAltModels; alt++ ) {

							if ( (Syst == "LY" || Syst == "MC_LY" || Syst == "SmEff_LY") && Runs[WhichRun] == "Run1" && AltModels[alt] == "_LYAttenuation") 
								{ continue;}

							if ( (Syst == "SCERecomb2" || Syst == "MC_SCERecomb2" || Syst == "SmEff_SCERecomb2") 
								&& Runs[WhichRun] == "Run1" && AltModels[alt] == "_CVextra") 
								{ continue;}	
						
							double CurrentFracCovEntry = FracCovariances[WhichRun][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1);
							double CurrentFracCovError = FracCovariances[WhichRun][WhichPlot]->GetBinError(WhichXBin+1,WhichYBin+1);

							double CurrentCovEntry = Covariances[WhichRun][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1);
							double CurrentCovError = Covariances[WhichRun][WhichPlot]->GetBinError(WhichXBin+1,WhichYBin+1);

							AltDataEntryX = AltBeamOnPlots[WhichPlot][alt]->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
							AltDataErrorX = AltBeamOnPlots[WhichPlot][alt]->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

							AltDataEntryY = AltBeamOnPlots[WhichPlot][alt]->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
							AltDataErrorY = AltBeamOnPlots[WhichPlot][alt]->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

//							if (Syst == "Flux" || Syst == "MC_Flux") {

//								AltDataEntryX = AltBeamOnPlots[WhichPlot][alt]->GetBinContent(WhichXBin+1) / (UniFlux[alt] * NTargets) * Units;
//								AltDataErrorX = AltBeamOnPlots[WhichPlot][alt]->GetBinError(WhichXBin+1) / (UniFlux[alt] * NTargets) * Units;

//								AltDataEntryY = AltBeamOnPlots[WhichPlot][alt]->GetBinContent(WhichYBin+1) / (UniFlux[alt] * NTargets) * Units;
//								AltDataErrorY = AltBeamOnPlots[WhichPlot][alt]->GetBinError(WhichYBin+1) / (UniFlux[alt] * NTargets) * Units;

//							}

							double LocalFracCovEntry = TMath::Max( ((AltDataEntryX - DataEntryX) / DataEntryX) * ( (AltDataEntryY - DataEntryY) / DataEntryY),1E-8);
							double LocalCovEntry = TMath::Max( (AltDataEntryX - DataEntryX) * (AltDataEntryY - DataEntryY),1E-8);

//							double LocalCovError = TMath::Max( TMath::Sqrt( 
//								TMath::Power(DataEntryY - AltDataEntryY,2.) * ( TMath::Power(DataErrorX,2.) + TMath::Power(AltDataErrorX,2.) ) +
//								TMath::Power(DataEntryX - AltDataEntryX,2.) * ( TMath::Power(DataErrorY,2.) + TMath::Power(AltDataErrorY,2.) ) ), 1E-10) ;

							double LocalFracCovError = 1E-8;
							double LocalCovError = 1E-8;

							if (AltUniverses[alt] > 2) {

								LocalFracCovEntry = LocalFracCovEntry * 1. / AltUniverses[alt];
								LocalFracCovError = LocalFracCovError * 1. / AltUniverses[alt];

								LocalCovEntry = LocalCovEntry * 1. / AltUniverses[alt];
								LocalCovError = LocalCovError * 1. / AltUniverses[alt];

							}

							CovFracEntry = CurrentFracCovEntry + LocalFracCovEntry;
							CovEntry = CurrentCovEntry + LocalCovEntry;

//							CovError = TMath::Sqrt( TMath::Power(CurrentCovError,2.) + TMath::Power(LocalCovError,2.) );
							CovFracError = 1E-8;
							CovError = 1E-8;

							FracCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovFracEntry);
							FracCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovFracError);

							Covariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovEntry);
							Covariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovError);

							// -------------------------------------------------------------------------------------------------------

//							if (	PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot" && 
//								(Syst == "LY" || Syst == "TPC" || Syst == "SCERecomb2" ||
//								 Syst == "MC_LY" || Syst == "MC_TPC" || Syst == "MC_SCERecomb2")
//							) {

//								TString CleanString = AltModels[alt];
//								CleanString = CleanString.ReplaceAll("_","");
//								cout << CleanString << " " << LocalCovEntry << endl;

//							}

							// -------------------------------------------------------------------------------------------------------

						}

					}

					else { cout << "Don't know how to handle this systematic uncentainty! Abort !"<< endl << endl; return; }

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
			
			// ---------------------------------------------------------------------------------------	

			// Plot the total covariance matrix GENIE CV

			if (BaseMC == "Overlay9") {	

				// ---------------------------------------------------------------------------------------------------------	

				// Store covariance matrices

				StoreCanvas(Covariances[WhichRun][WhichPlot], "", Syst, PlotNames[WhichPlot], BaseMC, Runs[WhichRun],Tune);

				// -------------------------------------------------------------------------------------------	

				// Store fractional covariance matrices

				StoreCanvas(FracCovariances[WhichRun][WhichPlot], "Frac", Syst, PlotNames[WhichPlot], BaseMC, Runs[WhichRun],Tune);

				// -------------------------------------------------------------------------------------------

			}		

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