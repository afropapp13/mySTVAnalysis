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

void WienerSVD_IndividualCovarianceMatrices(TString Syst = "None",TString Run = "Combined", TString BaseMC = "Overlay9",TString BeamOnSample = "BeamOn9",TString BeamOffSample = "ExtBNB9",TString DirtSample = "OverlayDirt9", TString Tune = "") {

	// -------------------------------------------------------------------------------------

	if (Syst == "None") { cout << "What the heck are you doing ? Specify the smearing / efficiency systematic uncertainty that you want to obtain!" << endl; return ;}

	// -------------------------------------------------------------------------------------

	GlobalSettings();
	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	TGaxis::SetMaxDigits(3);

	// -------------------------------------------------------------------------------------

	vector<TString> PlotNames;
	PlotNames.push_back("MuonCosThetaSingleBinPlot");
	PlotNames.push_back("DeltaPTPlot");	 

	const int N1DPlots = PlotNames.size();
		
	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;		
	Runs.push_back(Run);

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
	vector<vector <TH2D*> > CorrMatrices; CorrMatrices.resize(NRuns,vector<TH2D*>(N1DPlots));	
	vector<vector <TH2D*> > Covariances; Covariances.resize(NRuns,vector<TH2D*>(N1DPlots));	

	vector<vector <TH2D*> > SignalFracCovariances; SignalFracCovariances.resize(NRuns,vector<TH2D*>(N1DPlots));
	vector<vector <TH2D*> > SignalCorrMatrices; SignalCorrMatrices.resize(NRuns,vector<TH2D*>(N1DPlots));	
	vector<vector <TH2D*> > SignalCovariances; SignalCovariances.resize(NRuns,vector<TH2D*>(N1DPlots));
	
	vector<vector <TH2D*> > BkgFracCovariances; BkgFracCovariances.resize(NRuns,vector<TH2D*>(N1DPlots));
	vector<vector <TH2D*> > BkgCorrMatrices; BkgCorrMatrices.resize(NRuns,vector<TH2D*>(N1DPlots));	
	vector<vector <TH2D*> > BkgCovariances; BkgCovariances.resize(NRuns,vector<TH2D*>(N1DPlots));	

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

	// Overlay9_Combined

	if (Syst == "All_UBGenie") {

		UniAltModels.push_back("_All_UBGenie"); Universes.push_back(100);

	}

	//----------------------------------------//

	// Overlay9_Combined

	if (Syst == "AxFFCCQEshape_UBGenie") {

		UniAltModels.push_back("_AxFFCCQEshape_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Combined

	if (Syst == "DecayAngMEC_UBGenie") {

		UniAltModels.push_back("_DecayAngMEC_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Combined

	if (Syst == "NormCCCOH_UBGenie") {

		UniAltModels.push_back("_NormCCCOH_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Combined

	if (Syst == "NormNCCOH_UBGenie") {

		UniAltModels.push_back("_NormNCCOH_UBGenie"); Universes.push_back(2);

	}		

	//----------------------------------------//

	// Overlay9_Combined

	if (Syst == "RPA_CCQE_UBGenie") {

		UniAltModels.push_back("_RPA_CCQE_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Combined

	if (Syst == "ThetaDelta2NRad_UBGenie") {

		UniAltModels.push_back("_ThetaDelta2NRad_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Combined

	if (Syst == "Theta_Delta2Npi_UBGenie") {

		UniAltModels.push_back("_Theta_Delta2Npi_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Combined

	if (Syst == "VecFFCCQEshape_UBGenie") {

		UniAltModels.push_back("_VecFFCCQEshape_UBGenie"); Universes.push_back(2);

	}			

	//----------------------------------------//

	// Overlay9_Combined

	if (Syst == "XSecShape_CCMEC_UBGenie") {

		UniAltModels.push_back("_XSecShape_CCMEC_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "AGKYpT1pi_UBGenie") {

		UniAltModels.push_back("_AGKYpT1pi_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "AGKYxF1pi_UBGenie") {

		UniAltModels.push_back("_AGKYxF1pi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "AhtBY_UBGenie") {

		UniAltModels.push_back("_AhtBY_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "BhtBY_UBGenie") {

		UniAltModels.push_back("_BhtBY_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "CV1uBY_UBGenie") {

		UniAltModels.push_back("_CV1uBY_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "CV2uBY_UBGenie") {

		UniAltModels.push_back("_CV2uBY_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "EtaNCEL_UBGenie") {

		UniAltModels.push_back("_EtaNCEL_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "FrAbs_N_UBGenie") {

		UniAltModels.push_back("_FrAbs_N_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "FrAbs_pi_UBGenie") {

		UniAltModels.push_back("_FrAbs_pi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "FrCEx_N_UBGenie") {

		UniAltModels.push_back("_FrCEx_N_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "FrCEx_pi_UBGenie") {

		UniAltModels.push_back("_FrCEx_pi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "FrInel_N_UBGenie") {

		UniAltModels.push_back("_FrInel_N_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "FrInel_pi_UBGenie") {

		UniAltModels.push_back("_FrInel_pi_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "FrPiProd_N_UBGenie") {

		UniAltModels.push_back("_FrPiProd_N_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "FrPiProd_pi_UBGenie") {

		UniAltModels.push_back("_FrPiProd_pi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "FracDelta_CCMEC_UBGenie") {

		UniAltModels.push_back("_FracDelta_CCMEC_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "FracPN_CCMEC_UBGenie") {

		UniAltModels.push_back("_FracPN_CCMEC_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "MFP_N_UBGenie") {

		UniAltModels.push_back("_MFP_N_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "MFP_pi_UBGenie") {

		UniAltModels.push_back("_MFP_pi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "MaCCQE_UBGenie") {

		UniAltModels.push_back("_MaCCQE_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "MaCCRES_UBGenie") {

		UniAltModels.push_back("_MaCCRES_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "MaNCEL_UBGenie") {

		UniAltModels.push_back("_MaNCEL_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "MaNCRES_UBGenie") {

		UniAltModels.push_back("_MaNCRES_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "MvCCRES_UBGenie") {

		UniAltModels.push_back("_MvCCRES_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "MvNCRES_UBGenie") {

		UniAltModels.push_back("_MvNCRES_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvbarnCC1pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvbarnCC1pi_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvbarnCC2pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvbarnCC2pi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvbarnNC1pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvbarnNC1pi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvbarnNC2pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvbarnNC2pi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvbarpCC1pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvbarpCC1pi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvbarpCC2pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvbarpCC2pi_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvbarpNC1pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvbarpNC1pi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvbarpNC2pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvbarpNC2pi_UBGenie"); Universes.push_back(2);

	}																															

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvnCC1pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvnCC1pi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvnCC2pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvnCC2pi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvnNC1pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvnNC1pi_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvnNC2pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvnNC2pi_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvpCC1pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvpCC1pi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvpCC2pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvpCC2pi_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvpNC1pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvpNC1pi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NonRESBGvpNC2pi_UBGenie") {

		UniAltModels.push_back("_NonRESBGvpNC2pi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NormCCMEC_UBGenie") {

		UniAltModels.push_back("_NormCCMEC_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "NormNCMEC_UBGenie") {

		UniAltModels.push_back("_NormNCMEC_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "RDecBR1eta_UBGenie") {

		UniAltModels.push_back("_RDecBR1eta_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "RDecBR1gamma_UBGenie") {

		UniAltModels.push_back("_RDecBR1gamma_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "UnShortAxFFCCQEshape_UBGenie") {

		UniAltModels.push_back("_UnShortAxFFCCQEshape_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "UnShortDecayAngMEC_UBGenie") {

		UniAltModels.push_back("_UnShortDecayAngMEC_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "UnShortRPA_CCQE_UBGenie") {

		UniAltModels.push_back("_UnShortRPA_CCQE_UBGenie"); Universes.push_back(2);

	}	

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "UnShortTheta_Delta2Npi_UBGenie") {

		UniAltModels.push_back("_UnShortTheta_Delta2Npi_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "UnShortVecFFCCQEshape_UBGenie") {

		UniAltModels.push_back("_UnShortVecFFCCQEshape_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//

	// Overlay9_Run1_DecompXSecUnc

	if (Syst == "UnShortXSecShape_CCMEC_UBGenie") {

		UniAltModels.push_back("_UnShortXSecShape_CCMEC_UBGenie"); Universes.push_back(2);

	}

	//----------------------------------------//


	for (int UniAlt = 0; UniAlt < (int)(UniAltModels.size()); UniAlt++ ) {

		for (int Uni = 0; Uni < Universes[UniAlt]; Uni++ ) {

			AltModels.push_back(UniAltModels[UniAlt]+"_"+TString(std::to_string(Uni))); Colors.push_back(kGreen+2);
			AltUniverses.push_back(Universes[UniAlt]);

		}

	}			

	//----------------------------------------//

	int NAltModels = AltModels.size();

	vector< vector<TFile*> > AltMCFileSample; AltMCFileSample.resize(NRuns,vector<TFile*>(NAltModels));
	vector< vector<TFile*> > AltMCFileSampleResponseMatrix; AltMCFileSampleResponseMatrix.resize(NRuns,vector<TFile*>(NAltModels));

	//----------------------------------------//

	// Alternative Plots (also split into signal and background)

	vector <vector <TH1D*> > AltForwardFoldedCC1pPlots; AltForwardFoldedCC1pPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));
	vector <vector <TH1D*> > AltCC1pPlots; AltCC1pPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));
	vector <vector <TH1D*> > AltNonCC1pPlots; AltNonCC1pPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));
	vector <vector <TH1D*> > AltBeamOnPlots; AltBeamOnPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));
	vector <vector <TH1D*> > AltBeamOffPlots; AltBeamOffPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));

	vector <vector <TH1D*> > SignalAltBeamOnPlots; SignalAltBeamOnPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));

	vector <vector <TH1D*> > BkgAltBeamOnPlots; BkgAltBeamOnPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));	

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

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

		TString FileName = MigrationMatrixPath+Tune+"IndividualWienerSVD_"+Syst+"_CovarianceMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		if (BeamOnSample != "BeamOn9") { FileName = MigrationMatrixPath+BeamOnSample+"IndividualWienerSVD_"+Syst+"_CovarianceMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }

		TFile* FileCovarianceMatrices = new TFile(FileName,"recreate");
		
		// Open base files

		TString ExactFileLocation = PathToFiles+CutExtension;
		TString TStringBaseMC = ExactFileLocation+"/"+Tune+"STVStudies_"+BaseMC+"_"+Runs[WhichRun]+CutExtension+".root";
		TString TrueTStringBaseMC = PathToFiles+"/"+Tune+"TruthSTVAnalysis_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TString ResponseFileName = MigrationMatrixPath+Tune+"FileResponseMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";

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

			if (BeamOnSample == "BeamOn9" || BeamOnSample == "Overlay9") { // Beam On Sample, need to subtract ALL backgrounds

				for (int alt = 0; alt < NAltModels; alt++ ) {

					AltBeamOnPlots[WhichPlot][alt] = (TH1D*)AltForwardFoldedCC1pPlots[WhichPlot][alt]->Clone();
					AltBeamOnPlots[WhichPlot][alt]->Add(AltNonCC1pPlots[WhichPlot][alt]);	

					SignalAltBeamOnPlots[WhichPlot][alt] = (TH1D*)AltForwardFoldedCC1pPlots[WhichPlot][alt]->Clone();
					SignalAltBeamOnPlots[WhichPlot][alt]->Add(NonCC1pPlots[WhichPlot]);

					BkgAltBeamOnPlots[WhichPlot][alt] = (TH1D*)CC1pPlots[WhichPlot]->Clone();
					BkgAltBeamOnPlots[WhichPlot][alt]->Add(AltNonCC1pPlots[WhichPlot][alt]);								

				}		

			} 

			// -------------------------------------------------------------------------------------------------------
			// -------------------------------------------------------------------------------------------------------

			double BeamOnMax = 1.25*DataPlot->GetMaximum();

			// GENIE CV multisim overlays

			for (int unialt = 0; unialt < (int)(UniAltModels.size()); unialt++ ) {

				TString EventRatePlotCanvasName = Runs[WhichRun]+"_"+PlotNames[WhichPlot]+"_"+Syst;
				TCanvas* EventRatePlotCanvas = new TCanvas(EventRatePlotCanvasName,EventRatePlotCanvasName,205,34,1024,768);
				EventRatePlotCanvas->SetBottomMargin(0.17);
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

				TString EventRateCanvasName = "/"+Tune+"IndividualEventRate_WienerSVD_"+Syst+UniAltModels[unialt]+"_CovarianceMatrices_"+PlotNames[WhichPlot]+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf";
				EventRatePlotCanvas->SaveAs(PlotPath+BaseMC+EventRateCanvasName);
				delete EventRatePlotCanvas;

			} // End of the loop over a specific universe name

			// -------------------------------------------------------------------------------------------------------
			// -------------------------------------------------------------------------------------------------------
			
			int NBins = BeamOnPlots[WhichPlot]->GetXaxis()->GetNbins();
			const double* ArrayBins = BeamOnPlots[WhichPlot]->GetXaxis()->GetXbins()->GetArray();
			TString XTitle = BeamOnPlots[WhichPlot]->GetXaxis()->GetTitle();

			// -------------------------------------------------------------------------------------------------------

			// Declare the matrix & initialize the entries to 0

			FracCovariances[WhichRun][WhichPlot] = new TH2D(Syst+"_FracCovariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";i bin "+XTitle+";j bin "+XTitle,NBins,ArrayBins,NBins,ArrayBins);
			Covariances[WhichRun][WhichPlot] = new TH2D(Syst+"_Covariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";i bin "+XTitle+";j bin "+XTitle,NBins,ArrayBins,NBins,ArrayBins);
			
			SignalFracCovariances[WhichRun][WhichPlot] = new TH2D(Syst+"_SignalFracCovariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";i bin "+XTitle+";j bin "+XTitle,NBins,ArrayBins,NBins,ArrayBins);
			SignalCovariances[WhichRun][WhichPlot] = new TH2D(Syst+"_SignalCovariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";i bin "+XTitle+";j bin "+XTitle,NBins,ArrayBins,NBins,ArrayBins);
			
			BkgFracCovariances[WhichRun][WhichPlot] = new TH2D(Syst+"_BkgFracCovariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";i bin "+XTitle+";j bin "+XTitle,NBins,ArrayBins,NBins,ArrayBins);
			BkgCovariances[WhichRun][WhichPlot] = new TH2D(Syst+"_BkgCovariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";i bin "+XTitle+";j bin "+XTitle,NBins,ArrayBins,NBins,ArrayBins);						

			for (int WhichXBin = 0; WhichXBin < NBins; WhichXBin++) { 

				for (int WhichYBin = 0; WhichYBin < NBins; WhichYBin++) {

					FracCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,0.);
					FracCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,0.);

					Covariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,0.);
					Covariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,0.);	

					SignalFracCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,0.);
					SignalFracCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,0.);

					SignalCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,0.);
					SignalCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,0.);
					
					BkgFracCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,0.);
					BkgFracCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,0.);

					BkgCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,0.);
					BkgCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,0.);													

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

					double SignalAltDataEntryX = DataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
					double SignalAltDataErrorX = DataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

					double BkgAltDataEntryX = DataPlot->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
					double BkgAltDataErrorX = DataPlot->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;													

					// Y Bin entry / error

					double DataEntryY = DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
					double DataErrorY = DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

					double AltDataEntryY = DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
					double AltDataErrorY = DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

					double SignalAltDataEntryY = DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
					double SignalAltDataErrorY = DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

					double BkgAltDataEntryY = DataPlot->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
					double BkgAltDataErrorY = DataPlot->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;					

					//----------------------------------------//

					double CovFracEntry = -99.;
					double CovFracError = -99.;

					double CovEntry = -99.;
					double CovError = -99.;

					double SignalCovFracEntry = -99.;
					double SignalCovFracError = -99.;

					double SignalCovEntry = -99.;
					double SignalCovError = -99.;

					double BkgCovFracEntry = -99.;
					double BkgCovFracError = -99.;

					double BkgCovEntry = -99.;
					double BkgCovError = -99.;															

					for (int alt = 0; alt < NAltModels; alt++ ) {
								
						//----------------------------------------//	
							
						// Total covariance matrix							
						
						double CurrentFracCovEntry = FracCovariances[WhichRun][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1);
						double CurrentFracCovError = FracCovariances[WhichRun][WhichPlot]->GetBinError(WhichXBin+1,WhichYBin+1);

						double CurrentCovEntry = Covariances[WhichRun][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1);
						double CurrentCovError = Covariances[WhichRun][WhichPlot]->GetBinError(WhichXBin+1,WhichYBin+1);

						AltDataEntryX = AltBeamOnPlots[WhichPlot][alt]->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorX = AltBeamOnPlots[WhichPlot][alt]->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

						AltDataEntryY = AltBeamOnPlots[WhichPlot][alt]->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
						AltDataErrorY = AltBeamOnPlots[WhichPlot][alt]->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

						double LocalFracCovEntry = TMath::Max( ((AltDataEntryX - DataEntryX) / DataEntryX) * ( (AltDataEntryY - DataEntryY) / DataEntryY),1E-8);
						double LocalCovEntry = TMath::Max( (AltDataEntryX - DataEntryX) * (AltDataEntryY - DataEntryY),1E-8);

						double LocalFracCovError = 1E-8;
						double LocalCovError = 1E-8;

						if (AltUniverses[alt] > 2 || Run == "Run1_DecompXSecUnc") {						

							LocalFracCovEntry = LocalFracCovEntry * 1. / AltUniverses[alt];
							LocalFracCovError = LocalFracCovError * 1. / AltUniverses[alt];

							LocalCovEntry = LocalCovEntry * 1. / AltUniverses[alt];
							LocalCovError = LocalCovError * 1. / AltUniverses[alt];

						}

						CovFracEntry = CurrentFracCovEntry + LocalFracCovEntry;
						CovEntry = CurrentCovEntry + LocalCovEntry;

						CovFracError = 1E-8;
						CovError = 1E-8;

						FracCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovFracEntry);
						FracCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovFracError);

						Covariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovEntry);
						Covariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovError);

						//----------------------------------------//
							
						// Signal covariances
							
						double SignalCurrentFracCovEntry = SignalFracCovariances[WhichRun][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1);
						double SignalCurrentFracCovError = SignalFracCovariances[WhichRun][WhichPlot]->GetBinError(WhichXBin+1,WhichYBin+1);

						double SignalCurrentCovEntry = SignalCovariances[WhichRun][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1);
						double SignalCurrentCovError = SignalCovariances[WhichRun][WhichPlot]->GetBinError(WhichXBin+1,WhichYBin+1);

						SignalAltDataEntryX = SignalAltBeamOnPlots[WhichPlot][alt]->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
						SignalAltDataErrorX = SignalAltBeamOnPlots[WhichPlot][alt]->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

						SignalAltDataEntryY = SignalAltBeamOnPlots[WhichPlot][alt]->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
						SignalAltDataErrorY = SignalAltBeamOnPlots[WhichPlot][alt]->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

						double SignalLocalFracCovEntry = TMath::Max( ((SignalAltDataEntryX - DataEntryX) / DataEntryX) * ( (SignalAltDataEntryY - DataEntryY) / DataEntryY),1E-8);
						double SignalLocalCovEntry = TMath::Max( (SignalAltDataEntryX - DataEntryX) * (SignalAltDataEntryY - DataEntryY),1E-8);

						double SignalLocalFracCovError = 1E-8;
						double SignalLocalCovError = 1E-8;

						if (AltUniverses[alt] > 2 || Run == "Run1_DecompXSecUnc") {

							SignalLocalFracCovEntry = SignalLocalFracCovEntry * 1. / AltUniverses[alt];
							SignalLocalFracCovError = SignalLocalFracCovError * 1. / AltUniverses[alt];

							SignalLocalCovEntry = SignalLocalCovEntry * 1. / AltUniverses[alt];
							SignalLocalCovError = SignalLocalCovError * 1. / AltUniverses[alt];

						}

						SignalCovFracEntry = SignalCurrentFracCovEntry + SignalLocalFracCovEntry;
						SignalCovEntry = SignalCurrentCovEntry + SignalLocalCovEntry;

						SignalCovFracError = 1E-8;
						SignalCovError = 1E-8;

						SignalFracCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,SignalCovFracEntry);
						SignalFracCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,SignalCovFracError);

						SignalCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,SignalCovEntry);
						SignalCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,SignalCovError);
							
						//----------------------------------------//
							
						// Bkg covariances							
							
						double BkgCurrentFracCovEntry = BkgFracCovariances[WhichRun][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1);
						double BkgCurrentFracCovError = BkgFracCovariances[WhichRun][WhichPlot]->GetBinError(WhichXBin+1,WhichYBin+1);

						double BkgCurrentCovEntry = BkgCovariances[WhichRun][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1);
						double BkgCurrentCovError = BkgCovariances[WhichRun][WhichPlot]->GetBinError(WhichXBin+1,WhichYBin+1);

						BkgAltDataEntryX = BkgAltBeamOnPlots[WhichPlot][alt]->GetBinContent(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;
						BkgAltDataErrorX = BkgAltBeamOnPlots[WhichPlot][alt]->GetBinError(WhichXBin+1) / (IntegratedFlux * NTargets) * Units;

						BkgAltDataEntryY = BkgAltBeamOnPlots[WhichPlot][alt]->GetBinContent(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;
						BkgAltDataErrorY = BkgAltBeamOnPlots[WhichPlot][alt]->GetBinError(WhichYBin+1) / (IntegratedFlux * NTargets) * Units;

						double BkgLocalFracCovEntry = TMath::Max( ((BkgAltDataEntryX - DataEntryX) / DataEntryX) * ( (BkgAltDataEntryY - DataEntryY) / DataEntryY),1E-8);
						double BkgLocalCovEntry = TMath::Max( (BkgAltDataEntryX - DataEntryX) * (BkgAltDataEntryY - DataEntryY),1E-8);

						double BkgLocalFracCovError = 1E-8;
						double BkgLocalCovError = 1E-8;

						if (AltUniverses[alt] > 2 || Run == "Run1_DecompXSecUnc") {

							BkgLocalFracCovEntry = BkgLocalFracCovEntry * 1. / AltUniverses[alt];
							BkgLocalFracCovError = BkgLocalFracCovError * 1. / AltUniverses[alt];

							BkgLocalCovEntry = BkgLocalCovEntry * 1. / AltUniverses[alt];
							BkgLocalCovError = BkgLocalCovError * 1. / AltUniverses[alt];

						}

						BkgCovFracEntry = BkgCurrentFracCovEntry + BkgLocalFracCovEntry;
						BkgCovEntry = BkgCurrentCovEntry + BkgLocalCovEntry;

						BkgCovFracError = 1E-8;
						BkgCovError = 1E-8;

						BkgFracCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,BkgCovFracEntry);
						BkgFracCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,BkgCovFracError);

						BkgCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,BkgCovEntry);
						BkgCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,BkgCovError);																				
	
						//----------------------------------------//

					} // End of the loop over the alternative universes

					//----------------------------------------//

					// Setting the elements of the Cov Matrix

					FracCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovFracEntry);
					FracCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovFracError);

					Covariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovEntry);
					Covariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovError);

					//----------------------------------------//

					// Setting the elements of the signal Cov Matrix	

					SignalFracCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,SignalCovFracEntry);
					SignalFracCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,SignalCovFracError);

					SignalCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,SignalCovEntry);
					SignalCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,SignalCovError);

					//----------------------------------------//

					// Setting the elements of the bkg Cov Matrix	

					BkgFracCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,BkgCovFracEntry);
					BkgFracCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,BkgCovFracError);

					BkgCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,BkgCovEntry);
					BkgCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,BkgCovError);					

					//----------------------------------------//

				}  // end of the loop over bin Y

			} // end of the loop over bin X
			
			FileCovarianceMatrices->cd();

			FracCovariances[WhichRun][WhichPlot]->Write();			
			Covariances[WhichRun][WhichPlot]->Write();	

			SignalFracCovariances[WhichRun][WhichPlot]->Write();			
			SignalCovariances[WhichRun][WhichPlot]->Write();

			BkgFracCovariances[WhichRun][WhichPlot]->Write();			
			BkgCovariances[WhichRun][WhichPlot]->Write();								

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
