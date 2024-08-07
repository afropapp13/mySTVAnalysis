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

// --------------------------------------------------------------------------------------------------------------------------------------------

// TString Syst = "SmEff_Stat" "SmEff_POT" "SmEff_NTarget" "SmEff_LY" "SmEff_TPC" "SmEff_XSec" "SmEff_G4" "SmEff_Flux" "SmEff_Dirt"

void WienerSVD_SmearEff_CovarianceMatrices(TString Syst = "None",TString BaseMC = "Overlay9",TString BeamOnSample = "BeamOn9",TString BeamOffSample = "ExtBNB9",TString DirtSample = "OverlayDirt9") {

	// -------------------------------------------------------------------------------------

	if (Syst == "None") { cout << "What the heck are you doing ? Specify the smearing / efficiency systematic uncertainty that you want to obtain!" << endl; return ;}

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

	PlotNames.push_back("CCQEMuonMomentumPlot"); 
	PlotNames.push_back("CCQEMuonCosThetaPlot"); 
	PlotNames.push_back("CCQEProtonMomentumPlot"); 
	PlotNames.push_back("CCQEProtonCosThetaPlot");

	const int N1DPlots = PlotNames.size();
		
	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");			
	Runs.push_back("Combined");

	if (Syst == "LY" || Syst == "TPC" || Syst == "MC_LY" || Syst == "MC_TPC" || Syst == "SmEff_LY" || Syst == "SmEff_TPC" ) {

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

	if (Syst == "LY" || Syst == "MC_LY" || Syst == "SmEff_LY") {

		AltModels.push_back("_LYDown"); Colors.push_back(kRed+1); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_LYRayleigh"); Colors.push_back(kGreen+2); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_LYAttenuation"); Colors.push_back(kOrange+1); Universes.push_back(1); AltUniverses.push_back(1);

	}

	if (Syst == "TPC" || Syst == "MC_TPC" || Syst == "SmEff_TPC") {

		AltModels.push_back("_X"); Colors.push_back(kRed+1); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_YZ"); Colors.push_back(kGreen+2); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_ThetaXZ"); Colors.push_back(kOrange+1); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_ThetaYZ"); Colors.push_back(kBlue-3); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_SCE"); Colors.push_back(kBlue); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_Recombination2"); Colors.push_back(kMagenta); Universes.push_back(1); AltUniverses.push_back(1);
		//AltModels.push_back("_dEdx"); Colors.push_back(kYellow+2); Universes.push_back(1); AltUniverses.push_back(1);

	}

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

	if (Syst == "G4" || Syst == "MC_G4" || Syst == "SmEff_G4") {

		UniAltModels.push_back("_reinteractions"); Universes.push_back(100);
//		UniAltModels.push_back("_reinteractions_piminus_Geant4"); Universes.push_back(100);
//		UniAltModels.push_back("_reinteractions_piplus_Geant4"); Universes.push_back(100);
//		UniAltModels.push_back("_reinteractions_proton_Geant4"); Universes.push_back(100);

		for (int UniAlt = 0; UniAlt < (int)(UniAltModels.size()); UniAlt++ ) {

			for (int Uni = 0; Uni < Universes[UniAlt]; Uni++ ) {

				AltModels.push_back(UniAltModels[UniAlt]+"_"+TString(std::to_string(Uni))); Colors.push_back(kGreen+2);
				AltUniverses.push_back(Universes[UniAlt]);
			
			}

		}

	}

	if (Syst == "Flux" || Syst == "MC_Flux" || Syst == "SmEff_Flux") {

		UniAltModels.push_back("_fluxes"); Universes.push_back(100);

//		UniAltModels.push_back("_horncurrent_FluxUnisim"); Universes.push_back(100);
//		UniAltModels.push_back("_kminus_PrimaryHadronNormalization"); Universes.push_back(100);
//		UniAltModels.push_back("_kplus_PrimaryHadronFeynmanScaling"); Universes.push_back(100);
//		UniAltModels.push_back("_kzero_PrimaryHadronSanfordWang"); Universes.push_back(100);
//		UniAltModels.push_back("_nucleoninexsec_FluxUnisim"); Universes.push_back(100);
//		UniAltModels.push_back("_nucleonqexsec_FluxUnisim"); Universes.push_back(100);
//		UniAltModels.push_back("_nucleontotxsec_FluxUnisim"); Universes.push_back(100);
//		UniAltModels.push_back("_piminus_PrimaryHadronSWCentralSplineVariation"); Universes.push_back(100);
//		UniAltModels.push_back("_pioninexsec_FluxUnisim"); Universes.push_back(100);
//		UniAltModels.push_back("_pionqexsec_FluxUnisim"); Universes.push_back(100);
//		UniAltModels.push_back("_piontotxsec_FluxUnisim"); Universes.push_back(100);
//		UniAltModels.push_back("_piplus_PrimaryHadronSWCentralSplineVariation"); Universes.push_back(100);
//		UniAltModels.push_back("_expskin_FluxUnisim"); Universes.push_back(10);

		int NUniAltModels = (int)(UniAltModels.size());

		for (int UniAlt = 0; UniAlt < NUniAltModels; UniAlt++ ) {

			for (int Uni = 0; Uni < Universes[UniAlt]; Uni++ ) {

				AltModels.push_back(UniAltModels[UniAlt]+"_"+TString(std::to_string(Uni))); Colors.push_back(kGreen+2);
				AltUniverses.push_back(Universes[UniAlt]);

//				TString DublicateOverlaySample = UniAltModels[UniAlt];
//				TString ReducedOverlaySample = DublicateOverlaySample.ReplaceAll("m_","m");
//				if ( !(string(UniAltModels[UniAlt]).find("expskin") != std::string::npos) ) { ReducedOverlaySample = ReducedOverlaySample.ReplaceAll("n_","n"); }
//				ReducedOverlaySample = ReducedOverlaySample.ReplaceAll("g_","g");				
//				for (int i = 0; i < 10;i++) { ReducedOverlaySample.ReplaceAll(TString(std::to_string(i)),""); }
//				TString FluxHistoName = "numu_ms"+ReducedOverlaySample+"/hEnumu"+ReducedOverlaySample+"_ms_"+TString(std::to_string(Uni));

				TString FluxHistoName = "numu_ms_total/hEnumu_ms_"+TString(std::to_string(Uni));
				UniHistoFlux.push_back( (TH1D*)(FluxFile->Get(FluxHistoName) ) );
	
			}

		}

	}

	// -------------------------------------------------------------------------------------------------------------------------------------

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

		// Until Run3 NuWro is produced
		if ( Runs[WhichRun] == "Run3" && (BaseMC == "Overlay9NuWro" || BeamOnSample == "Overlay9NuWro") ) 
			{ continue;}

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------
	
		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);						
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface) * (SoftFidSurface / Nominal_UB_XY_Surface);

		if (Syst == "Flux" || Syst == "MC_Flux" || Syst == "SmEff_Flux") {

			UniFlux.clear();

			for (int UniAlt = 0; UniAlt < NAltModels; UniAlt++ ) {

				double UniFluxPOT = (UniHistoFlux[UniAlt]->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface) * (SoftFidSurface / Nominal_UB_XY_Surface);
				UniFlux.push_back(UniFluxPOT);

			}

		}

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------

		TString FileName = MigrationMatrixPath+"WienerSVD_"+Syst+"_CovarianceMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		if (BeamOnSample != "BeamOn9") { MigrationMatrixPath+BeamOnSample+"WienerSVD_"+Syst+"_CovarianceMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }

		TFile* FileCovarianceMatrices = new TFile(FileName,"recreate");
		
		// Open base files

		TString ExactFileLocation = PathToFiles+CutExtension;
		TString TStringBaseMC = ExactFileLocation+"/STVStudies_"+BaseMC+"_"+Runs[WhichRun]+CutExtension+".root";
		TString TrueTStringBaseMC = PathToFiles+"/TruthSTVAnalysis_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TString ResponseFileName = MigrationMatrixPath+"FileResponseMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";

		if (Syst == "LY" || Syst == "TPC" || Syst == "MC_LY" || Syst == "MC_TPC"  || Syst == "SmEff_LY" || Syst == "SmEff_TPC") { 

			TStringBaseMC = ExactFileLocation+"/STVStudies_"+BaseMC+"_"+Runs[WhichRun]+"_CV"+CutExtension+".root"; 
			TrueTStringBaseMC = PathToFiles+"/TruthSTVAnalysis_"+BaseMC+"_"+Runs[WhichRun]+"_CV_"+UBCodeVersion+".root"; 
			ResponseFileName = MigrationMatrixPath+"FileResponseMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_CV_"+UBCodeVersion+".root";

		}

		TFile* FileResponseMatrices = new TFile(ResponseFileName,"readonly");
		TrueMCFileSample[WhichRun] = TFile::Open(TrueTStringBaseMC,"readonly");
		MCFileSample[WhichRun] = TFile::Open(TStringBaseMC,"readonly");
		BeamOnFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOnSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");
		BeamOffFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOffSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");
		DirtFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+DirtSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");

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

				TString TStringAltBaseMC = ExactFileLocation+"/STVStudies_"+BaseMC+"_"+Runs[WhichRun]+AltModels[alt]+CutExtension+".root";
				AltMCFileSample[WhichRun][alt] = TFile::Open(TStringAltBaseMC,"readonly");

				AltCC1pPlots[WhichPlot][alt] = (TH1D*)(AltMCFileSample[WhichRun][alt]->Get("CC1pReco"+PlotNames[WhichPlot]));
				AltNonCC1pPlots[WhichPlot][alt] = (TH1D*)(AltMCFileSample[WhichRun][alt]->Get("NonCC1pReco"+PlotNames[WhichPlot]));

				TString TStringAltBaseMCResponseMatrix = MigrationMatrixPath+"FileResponseMatrices_"+BaseMC+"_"+Runs[WhichRun]+AltModels[alt]+"_"+UBCodeVersion+".root";
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

			TH1D* DataPlot = (TH1D*)BeamOnPlots[WhichPlot]->Clone();
			TH1D* AltDataPlot = (TH1D*)BeamOnPlots[WhichPlot]->Clone();

//			if (BeamOnSample == "BeamOn9" || BeamOnSample == "Overlay9") { // Beam On Sample, need to subtract ALL backgrounds

//				DataPlot->Add(NonCC1pPlots[WhichPlot],-1.);
//				DataPlot->Add(BeamOffPlots[WhichPlot],-1.);
//				DataPlot->Add(DirtPlots[WhichPlot],-1.);

//				for (int alt = 0; alt < NAltModels; alt++ ) {

//					if ( (Syst == "LY" || Syst == "MC_LY"  || Syst == "SmEff_LY") && Runs[WhichRun] == "Run1" && AltModels[alt] == "_LYAttenuation") 
//						{ continue;}

//					AltBeamOnPlots[WhichPlot][alt] = (TH1D*)BeamOnPlots[WhichPlot]->Clone();
//					AltBeamOnPlots[WhichPlot][alt]->Add(AltNonCC1pPlots[WhichPlot][alt],-1.);
//					AltBeamOnPlots[WhichPlot][alt]->Add(BeamOffPlots[WhichPlot],-1.);
//					AltBeamOnPlots[WhichPlot][alt]->Add(DirtPlots[WhichPlot],-1.);

//				}

//				if (Syst == "Dirt" || Syst == "MC_Dirt" || Syst == "SmEff_Dirt") {

//					AltDataPlot->Add(NonCC1pPlots[WhichPlot],-1.);
//					AltDataPlot->Add(BeamOffPlots[WhichPlot],-1.);

//					// Scale the dirt sample by 25% down
//					TH1D* DirtClone = (TH1D*)(DirtPlots[WhichPlot]->Clone());
//					DirtClone->Scale(0.75);
//					AltDataPlot->Add(DirtClone,-1.);

//				}

//			} else { // Fake Data Studies, working with MC CC1p signal events as BeamOn

				DataPlot = CC1pPlots[WhichPlot];

				for (int alt = 0; alt < NAltModels; alt++ ) {

					if ( (Syst == "LY" || Syst == "MC_LY" || Syst == "SmEff_LY") && Runs[WhichRun] == "Run1" && AltModels[alt] == "_LYAttenuation") 
						{ continue;}

//					AltBeamOnPlots[WhichPlot][alt] = AltCC1pPlots[WhichPlot][alt];
					AltBeamOnPlots[WhichPlot][alt] = AltForwardFoldedCC1pPlots[WhichPlot][alt];

				}

			//}

			// -------------------------------------------------------------------------------------------------------
			// -------------------------------------------------------------------------------------------------------

			double BeamOnMax = 1.25*DataPlot->GetMaximum();

			// GENIE CV LY / TPC overlays

			if (BaseMC == "Overlay9" && (Syst == "LY" || Syst == "TPC" || Syst == "SmEff_LY" || Syst == "SmEff_TPC") ) {

				TString EventRatePlotCanvasName = Runs[WhichRun]+"_"+PlotNames[WhichPlot]+"_"+Syst;
				TCanvas* EventRatePlotCanvas = new TCanvas(EventRatePlotCanvasName,EventRatePlotCanvasName,205,34,1024,768);
				EventRatePlotCanvas->SetBottomMargin(0.16);
				EventRatePlotCanvas->SetLeftMargin(0.15);

				TLegend* leg = new TLegend(0.12,0.91,0.9,0.99);
				leg->SetBorderSize(0);
				leg->SetTextFont(FontStyle);
				leg->SetTextSize(TextSize-0.02);
				if (Syst == "TPC" || Syst == "MC_TPC" || Syst == "SmEff_TPC") { leg->SetTextSize(TextSize-0.03); }
				leg->SetNColumns(3);

				DataPlot->GetXaxis()->CenterTitle();
				DataPlot->GetXaxis()->SetTitleFont(FontStyle);
				DataPlot->GetXaxis()->SetTitleSize(TextSize);
				DataPlot->GetXaxis()->SetLabelFont(FontStyle);
				DataPlot->GetXaxis()->SetLabelSize(TextSize);
				DataPlot->GetXaxis()->SetNdivisions(6);

				DataPlot->GetYaxis()->SetTitle("# Events / "+ToString(DataPOT));
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

				// Closure Test 
				ClosureTestCC1pPlots[WhichPlot]->SetMarkerColor(kMagenta);
				ClosureTestCC1pPlots[WhichPlot]->Draw("p0 hist same");

				leg->Draw();

				TString EventRateCanvasName = "/EventRate_WienerSVD_"+Syst+"_CovarianceMatrices_"+PlotNames[WhichPlot]+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf";
				EventRatePlotCanvas->SaveAs(PlotPath+BaseMC+EventRateCanvasName);
				delete EventRatePlotCanvas;

			} // End of the cases where we overlay the event rates for GENIE CV  LY / TPC systematics

			// -------------------------------------------------------------------------------------------------------

			// GENIE CV multisim overlays

			if (BaseMC == "Overlay9" && (Syst == "XSec" || Syst == "G4" || Syst == "Flux" || Syst == "SmEff_XSec" || Syst == "SmEff_G4" || Syst == "SmEff_Flux") ) {

				for (int unialt = 0; unialt < (int)(UniAltModels.size()); unialt++ ) {

					TString EventRatePlotCanvasName = Runs[WhichRun]+"_"+PlotNames[WhichPlot]+"_"+Syst;
					TCanvas* EventRatePlotCanvas = new TCanvas(EventRatePlotCanvasName,EventRatePlotCanvasName,205,34,1024,768);
					EventRatePlotCanvas->SetBottomMargin(0.16);
					EventRatePlotCanvas->SetLeftMargin(0.15);

					TLegend* leg = new TLegend(0.12,0.91,0.9,0.99);
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

					DataPlot->GetYaxis()->SetTitle("# Events / "+ToString(DataPOT));
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
					ClosureTestCC1pPlots[WhichPlot]->Draw("p0 hist same");

					leg->Draw();

					TString EventRateCanvasName = "/EventRate_WienerSVD_"+Syst+UniAltModels[unialt]+"_CovarianceMatrices_"+PlotNames[WhichPlot]+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf";
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

			Covariances[WhichRun][WhichPlot] = new TH2D(Syst+"_Covariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";i bin "+XTitle+";j bin "+XTitle,NBins,ArrayBins,NBins,ArrayBins);

			for (int WhichXBin = 0; WhichXBin < NBins; WhichXBin++) { 

				for (int WhichYBin = 0; WhichYBin < NBins; WhichYBin++) {

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

					double CovEntry = -99.;
					double CovError = -99.;

					// Based on the type of systematic unceratainty, choose how to handle it

					if (Syst == "NTarget" || Syst == "MC_NTarget"  || Syst == "SmEff_NTarget") {

						AltDataEntryX = (1+NTargetUncertainty) * DataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorX = (1+NTargetUncertainty) * DataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

						AltDataEntryY = (1+NTargetUncertainty) * DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorY = (1+NTargetUncertainty) * DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

						CovEntry = TMath::Max( ( (AltDataEntryX - DataEntryX) / DataEntryX) * ( (AltDataEntryY - DataEntryY) / DataEntryY),1E-8);
						// CovError = TMath::Max( TMath::Sqrt( 
						// 	TMath::Power(DataEntryY - AltDataEntryY,2.) * ( TMath::Power(DataErrorX,2.) + TMath::Power(AltDataErrorX,2.) ) +
						// 	TMath::Power(DataEntryX - AltDataEntryX,2.) * ( TMath::Power(DataErrorY,2.) + TMath::Power(AltDataErrorY,2.) ) ), 1E-10) ;
						CovError = 1E-8;


					} else if (Syst == "POT" || Syst == "MC_POT" || Syst == "SmEff_POT") {

						AltDataEntryX = (1+POTUncertainty) * DataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorX = (1+POTUncertainty) * DataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

						AltDataEntryY = (1+POTUncertainty) * DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorY = (1+POTUncertainty) * DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

						CovEntry = TMath::Max( ((AltDataEntryX - DataEntryX) / DataEntryX) * ( (AltDataEntryY - DataEntryY) / DataEntryY ),1E-8);
//						CovError = TMath::Max( TMath::Sqrt( 
//							TMath::Power(DataEntryY - AltDataEntryY,2.) * ( TMath::Power(DataErrorX,2.) + TMath::Power(AltDataErrorX,2.) ) +
//							TMath::Power(DataEntryX - AltDataEntryX,2.) * ( TMath::Power(DataErrorY,2.) + TMath::Power(AltDataErrorY,2.) ) ), 1E-10) ;
						CovError = 1E-8;


					} else if (Syst == "Dirt" || Syst == "MC_Dirt" || Syst == "SmEff_Dirt") {

						AltDataEntryX = AltDataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorX = AltDataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

						AltDataEntryY = AltDataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorY = AltDataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

						CovEntry = TMath::Max( ((AltDataEntryX - DataEntryX) / DataEntryX ) * ( (AltDataEntryY - DataEntryY) / DataEntryY ),1E-8);
//						CovError = TMath::Max( TMath::Sqrt( 
//							TMath::Power(DataEntryY - AltDataEntryY,2.) * ( TMath::Power(DataErrorX,2.) + TMath::Power(AltDataErrorX,2.) ) +
//							TMath::Power(DataEntryX - AltDataEntryX,2.) * ( TMath::Power(DataErrorY,2.) + TMath::Power(AltDataErrorY,2.) ) ), 1E-10) ;

						CovError = 1E-8;

					} else if (Syst == "Stat" || Syst == "MC_Stat" || Syst == "SmEff_Stat") {

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

						CovEntry = TMath::Max( ((AltDataEntryX - DataEntryX) / DataEntryXCV) * ( (AltDataEntryY - DataEntryY) / DataEntryYCV),1E-8);
//						CovError = TMath::Max( TMath::Sqrt( 
//							TMath::Power(DataEntryY - AltDataEntryY,2.) * ( TMath::Power(DataErrorX,2.) + TMath::Power(AltDataErrorX,2.) ) +
//							TMath::Power(DataEntryX - AltDataEntryX,2.) * ( TMath::Power(DataErrorY,2.) + TMath::Power(AltDataErrorY,2.) ) ), 1E-10) ;

						CovError = 1E-8;

					} else if (
						Syst == "LY" || Syst == "TPC" || Syst == "XSec" || Syst == "G4" || Syst == "Flux" ||
						Syst == "MC_LY" || Syst == "MC_TPC" || Syst == "MC_XSec" || Syst == "MC_G4" || Syst == "MC_Flux" ||
						Syst == "SmEff_LY" || Syst == "SmEff_TPC" || Syst == "SmEff_XSec" || Syst == "SmEff_G4" || Syst == "SmEff_Flux"
					) {						

						for (int alt = 0; alt < NAltModels; alt++ ) {

							if ( (Syst == "LY" || Syst == "MC_LY" || Syst == "SmEff_LY") && Runs[WhichRun] == "Run1" && AltModels[alt] == "_LYAttenuation") 
								{ continue;}							

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

							double LocalCovEntry = TMath::Max( ((AltDataEntryX - DataEntryX) / DataEntryX) * ( (AltDataEntryY - DataEntryY) / DataEntryY),1E-8);

//							double LocalCovError = TMath::Max( TMath::Sqrt( 
//								TMath::Power(DataEntryY - AltDataEntryY,2.) * ( TMath::Power(DataErrorX,2.) + TMath::Power(AltDataErrorX,2.) ) +
//								TMath::Power(DataEntryX - AltDataEntryX,2.) * ( TMath::Power(DataErrorY,2.) + TMath::Power(AltDataErrorY,2.) ) ), 1E-10) ;

							double LocalCovError = 1E-8;

							if (AltUniverses[alt] > 2) {

								LocalCovEntry = LocalCovEntry * 1. / AltUniverses[alt];
								LocalCovError = LocalCovError * 1. / AltUniverses[alt];

							}

							CovEntry = CurrentCovEntry + LocalCovEntry;
//							CovError = TMath::Sqrt( TMath::Power(CurrentCovError,2.) + TMath::Power(LocalCovError,2.) );
							CovError = 1E-8;

							Covariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovEntry);
							Covariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovError);

						}

					}

					else { cout << "Don't know how to handle this systematic uncentainty! Abort !"<< endl << endl; return; }

					// -------------------------------------------------------------------------------------------------------

					// Setting the elements of the Cov Matrix

					Covariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovEntry);
					Covariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovError);

					// -------------------------------------------------------------------------------------------------------

				}  // end of the loop over bin Y

			} // end of the loop over bin X
			
			FileCovarianceMatrices->cd();
			Covariances[WhichRun][WhichPlot]->Write();
			
			// ---------------------------------------------------------------------------------------	

			// Plot the total covariance matrix GENIE CV

			if (BaseMC == "Overlay9") {		
		
				TString CanvasName = Syst+"_"+PlotNames[WhichPlot]+BaseMC+"_"+Runs[WhichRun];
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

				Covariances[WhichRun][WhichPlot]->SetTitle(Runs[WhichRun] + " " + Syst);	

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
				
				PlotCanvas->SaveAs(PlotPath+BaseMC+"/WienerSVD_"+Syst+"_CovarianceMatrices_"+PlotNames[WhichPlot]+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
				
				delete PlotCanvas;	

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
