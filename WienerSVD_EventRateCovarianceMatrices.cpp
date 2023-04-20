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

// ------------------------------------------------------------------------------------------------------------------------------

// TString Syst = "Stat" "POT" "NTarget" "LY" "TPC" "SCERecomb2" "XSec" "DetailedXSec" "G4" "Flux" "Dirt" "MC_Stat" "NuWro"

void WienerSVD_EventRateCovarianceMatrices(TString Syst = "None",TString BaseMC = "Overlay9",TString BeamOnSample = "BeamOn9",TString BeamOffSample = "ExtBNB9",TString DirtSample = "OverlayDirt9", TString Tune = "") {

	// -------------------------------------------------------------------------------------

	if (Syst == "None") { cout << "What the heck are you doing ? Specify the smearing / efficiency systematic uncertainty that you want to obtain!" << endl; return ;}

	// -------------------------------------------------------------------------------------

	GlobalSettings();
	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	TGaxis::SetMaxDigits(3);

	TString CutExtension = "_NoCuts";

	// -------------------------------------------------------------------------------------

	vector<TString> PlotNames;
	PlotNames.push_back("LLRPIDPlot");
	PlotNames.push_back("MuonLLRPIDPlot");
	PlotNames.push_back("ProtonLLRPIDPlot");

	const int N1DPlots = PlotNames.size();
		
	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;		
	Runs.push_back("Combined");

	if (Syst == "LY" || Syst == "TPC" || Syst == "SCERecomb2" || Syst == "MC_LY" || Syst == "MC_TPC" || Syst == "MC_SCERecomb2" 
	|| Syst == "SmEff_LY" || Syst == "SmEff_TPC" || Syst == "SmEff_SCERecomb2" ) {

		Runs.clear();
		Runs.push_back("Run3");

	}

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
	
	vector<vector <TH2D*> > FracCovariances; FracCovariances.resize(NRuns,vector<TH2D*>(N1DPlots));
	vector<vector <TH2D*> > Covariances; Covariances.resize(NRuns,vector<TH2D*>(N1DPlots));	

	// -------------------------------------------------------------------------------------------------------------------------------------

	// Alternative MC Models

	std::vector<TString> AltModels;
	std::vector<int> AltUniverses;
	std::vector<int> Universes;
	std::vector<TString> UniAltModels;

	//----------------------------------------//

	if (Syst == "LY" || Syst == "MC_LY" || Syst == "SmEff_LY") {

		AltModels.push_back("_LYDown"); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_LYRayleigh"); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_LYAttenuation"); Universes.push_back(1); AltUniverses.push_back(1);

	}

	//----------------------------------------//	

	if (Syst == "TPC" || Syst == "MC_TPC" || Syst == "SmEff_TPC") {

		AltModels.push_back("_X"); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_YZ"); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_ThetaXZ"); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_ThetaYZ"); Universes.push_back(1); AltUniverses.push_back(1);

	}

	//----------------------------------------//	

	if (Syst == "SCERecomb2" || Syst == "MC_SCERecomb2" || Syst == "SmEff_SCERecomb2") {

		AltModels.push_back("_SCE"); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_Recombination2"); Universes.push_back(1); AltUniverses.push_back(1);

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

				AltModels.push_back(UniAltModels[UniAlt]+"_"+TString(std::to_string(Uni)));
	 			AltUniverses.push_back(Universes[UniAlt]);

			}

		}

	}	

	//----------------------------------------//	

	if (Syst == "G4" || Syst == "MC_G4" || Syst == "SmEff_G4") {

		UniAltModels.push_back("_reinteractions"); Universes.push_back(100);

		for (int UniAlt = 0; UniAlt < (int)(UniAltModels.size()); UniAlt++ ) {

			for (int Uni = 0; Uni < Universes[UniAlt]; Uni++ ) {

				AltModels.push_back(UniAltModels[UniAlt]+"_"+TString(std::to_string(Uni)));
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

				AltModels.push_back(UniAltModels[UniAlt]+"_"+TString(std::to_string(Uni)));
				AltUniverses.push_back(Universes[UniAlt]);
	
			}

		}

	}

	//----------------------------------------//	

	if (Syst == "MC_Stat") {

		UniAltModels.push_back("_MC_Stat"); Universes.push_back(100);

		for (int UniAlt = 0; UniAlt < (int)(UniAltModels.size()); UniAlt++ ) {

			for (int Uni = 0; Uni < Universes[UniAlt]; Uni++ ) {

				AltModels.push_back(UniAltModels[UniAlt]+"_"+TString(std::to_string(Uni)));
				AltUniverses.push_back(Universes[UniAlt]);
			
			}

		}

	}	

	//----------------------------------------//	

	if (Syst == "NuWro") {

		UniAltModels.push_back("_NuWro"); Universes.push_back(1);

		for (int UniAlt = 0; UniAlt < (int)(UniAltModels.size()); UniAlt++ ) {

			for (int Uni = 0; Uni < Universes[UniAlt]; Uni++ ) {

				AltModels.push_back(UniAltModels[UniAlt]+"_"+TString(std::to_string(Uni)));
				AltUniverses.push_back(Universes[UniAlt]);
			
			}

		}

	}	

	//----------------------------------------//

	int NAltModels = AltModels.size();

	vector< vector<TFile*> > AltMCFileSample; AltMCFileSample.resize(NRuns,vector<TFile*>(NAltModels));

	//----------------------------------------//

	// Alternative Plots (also split into signal and background)

	vector <vector <TH1D*> > AltCC1pPlots; AltCC1pPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));
	vector <vector <TH1D*> > AltNonCC1pPlots; AltNonCC1pPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));
	vector <vector <TH1D*> > AltBeamOnPlots; AltBeamOnPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));
	vector <vector <TH1D*> > AltBeamOffPlots; AltBeamOffPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------

		TString FileName = MigrationMatrixPath+Tune+"WienerSVD_"+Syst+"_EventRateCovarianceMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		if (BeamOnSample != "BeamOn9") { FileName = MigrationMatrixPath+BeamOnSample+"WienerSVD_"+Syst+"_EventRateCovarianceMatrices_"+BaseMC+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }

		TFile* FileCovarianceMatrices = new TFile(FileName,"recreate");
		
		// Open base files

		TString ExactFileLocation = PathToFiles+CutExtension;
		TString TStringBaseMC = ExactFileLocation+"/"+Tune+"STVStudies_"+BaseMC+"_"+Runs[WhichRun]+CutExtension+".root";

		if (Syst == "LY" || Syst == "TPC" || Syst == "MC_LY" || Syst == "MC_TPC"  || Syst == "SmEff_LY" || Syst == "SmEff_TPC") { 

			TStringBaseMC = ExactFileLocation+"/"+Tune+"STVStudies_"+BaseMC+"_"+Runs[WhichRun]+"_CV"+CutExtension+".root"; 

		}

		if (Syst == "SCERecomb2" || Syst == "MC_SCERecomb2"  || Syst == "SmEff_SCERecomb2") { 

			TStringBaseMC = ExactFileLocation+"/"+Tune+"STVStudies_"+BaseMC+"_"+Runs[WhichRun]+"_CVextra"+CutExtension+".root"; 

		}

		MCFileSample[WhichRun] = TFile::Open(TStringBaseMC,"readonly");
		BeamOnFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOnSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");
		BeamOffFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOffSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");
		DirtFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/"+Tune+"STVStudies_"+DirtSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");

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

			// Grab alternative plots

			for (int alt = 0; alt < NAltModels; alt++ ) {

				if ( (Syst == "LY" || Syst == "MC_LY" || Syst == "SmEff_LY") && Runs[WhichRun] == "Run1" && AltModels[alt] == "_LYAttenuation") 
					{ continue;}

				// Open Alternative MC files

				TString TStringAltBaseMC = ExactFileLocation+"/"+Tune+"STVStudies_"+BaseMC+"_"+Runs[WhichRun]+AltModels[alt]+CutExtension+".root";
				if (Syst == "NuWro") { TStringAltBaseMC = ExactFileLocation+"/"+Tune+"STVStudies_Overlay9NuWro_"+Runs[WhichRun]+CutExtension+".root"; }
				AltMCFileSample[WhichRun][alt] = TFile::Open(TStringAltBaseMC,"readonly");

				AltCC1pPlots[WhichPlot][alt] = (TH1D*)(AltMCFileSample[WhichRun][alt]->Get("CC1pReco"+PlotNames[WhichPlot]));
				AltNonCC1pPlots[WhichPlot][alt] = (TH1D*)(AltMCFileSample[WhichRun][alt]->Get("NonCC1pReco"+PlotNames[WhichPlot]));			

				AltCC1pPlots[WhichPlot][alt]->SetDirectory(0); // to decouple it from the open file directory
				AltNonCC1pPlots[WhichPlot][alt]->SetDirectory(0); // to decouple it from the open file directory

				AltMCFileSample[WhichRun][alt]->Close();

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

					AltBeamOnPlots[WhichPlot][alt] = (TH1D*)AltCC1pPlots[WhichPlot][alt]->Clone();
					AltBeamOnPlots[WhichPlot][alt]->Add(AltNonCC1pPlots[WhichPlot][alt]);

					if (Syst == "MC_Stat") {

						AltBeamOnPlots[WhichPlot][alt] = (TH1D*)AltCC1pPlots[WhichPlot][alt]->Clone();

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

				}				

			} 

			// -------------------------------------------------------------------------------------------------------
			
			int NBins = BeamOnPlots[WhichPlot]->GetXaxis()->GetNbins();
			const double* ArrayBins = BeamOnPlots[WhichPlot]->GetXaxis()->GetXbins()->GetArray();
			TString XTitle = BeamOnPlots[WhichPlot]->GetXaxis()->GetTitle();

			// -------------------------------------------------------------------------------------------------------

			// Declare the matrix & initialize the entries to 0

			FracCovariances[WhichRun][WhichPlot] = new TH2D(Syst+"_FracCovariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";i bin "+XTitle+";j bin "+XTitle,NBins,ArrayBins,NBins,ArrayBins);
			Covariances[WhichRun][WhichPlot] = new TH2D(Syst+"_Covariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";i bin "+XTitle+";j bin "+XTitle,NBins,ArrayBins,NBins,ArrayBins);

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

					double DataEntryX = DataPlot->GetBinContent(WhichXBin+1);
					double DataErrorX = DataPlot->GetBinError(WhichXBin+1);

					double AltDataEntryX = DataPlot->GetBinContent(WhichXBin+1);
					double AltDataErrorX = DataPlot->GetBinError(WhichXBin+1);								

					// Y Bin entry / error

					double DataEntryY = DataPlot->GetBinContent(WhichYBin+1);
					double DataErrorY = DataPlot->GetBinError(WhichYBin+1);

					double AltDataEntryY = DataPlot->GetBinContent(WhichYBin+1);
					double AltDataErrorY = DataPlot->GetBinError(WhichYBin+1);			

					// -------------------------------------------------------------------------------------------------------
					// -------------------------------------------------------------------------------------------------------

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

					// Based on the type of systematic unceratainty, choose how to handle it

					if (Syst == "NTarget" || Syst == "MC_NTarget"  || Syst == "SmEff_NTarget") {

						AltDataEntryX = (1+NTargetUncertainty) * DataPlot->GetBinContent(WhichXBin+1);
						AltDataErrorX = (1+NTargetUncertainty) * DataPlot->GetBinError(WhichXBin+1);

						AltDataEntryY = (1+NTargetUncertainty) * DataPlot->GetBinContent(WhichYBin+1);
						AltDataErrorY = (1+NTargetUncertainty) * DataPlot->GetBinError(WhichYBin+1);

						CovFracEntry = TMath::Max( ( (AltDataEntryX - DataEntryX) / DataEntryX) * ( (AltDataEntryY - DataEntryY) / DataEntryY),1E-8);
						CovEntry = TMath::Max( (AltDataEntryX - DataEntryX) * (AltDataEntryY - DataEntryY),1E-8);

						CovFracError = 1E-8;
						CovError = 1E-8;

					} else if (Syst == "POT" || Syst == "MC_POT" || Syst == "SmEff_POT") {

						AltDataEntryX = (1+POTUncertainty) * DataPlot->GetBinContent(WhichXBin+1);
						AltDataErrorX = (1+POTUncertainty) * DataPlot->GetBinError(WhichXBin+1);

						AltDataEntryY = (1+POTUncertainty) * DataPlot->GetBinContent(WhichYBin+1);
						AltDataErrorY = (1+POTUncertainty) * DataPlot->GetBinError(WhichYBin+1);

						CovFracEntry = TMath::Max( ((AltDataEntryX - DataEntryX) / DataEntryX) * ( (AltDataEntryY - DataEntryY) / DataEntryY ),1E-8);
						CovEntry = TMath::Max( (AltDataEntryX - DataEntryX) * (AltDataEntryY - DataEntryY),1E-8);

						CovFracError = 1E-8;
						CovError = 1E-8;

					} else if (Syst == "Dirt" || Syst == "MC_Dirt" || Syst == "SmEff_Dirt") {

						AltDataEntryX = AltDataPlot->GetBinContent(WhichXBin+1);
						AltDataErrorX = AltDataPlot->GetBinError(WhichXBin+1);

						AltDataEntryY = AltDataPlot->GetBinContent(WhichYBin+1);
						AltDataErrorY = AltDataPlot->GetBinError(WhichYBin+1);

						CovFracEntry = TMath::Max( ((AltDataEntryX - DataEntryX) / DataEntryX ) * ( (AltDataEntryY - DataEntryY) / DataEntryY ),1E-8);
						CovEntry = TMath::Max( (AltDataEntryX - DataEntryX) * (AltDataEntryY - DataEntryY),1E-8);

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

						CovFracError = 1E-8;
						CovError = 1E-8;

					} else if (
						Syst == "LY" || Syst == "TPC" || Syst == "SCERecomb2" || Syst == "XSec" || Syst == "DetailedXSec" || 
						Syst == "G4" || Syst == "Flux" ||
						Syst == "MC_LY" || Syst == "MC_TPC" || Syst == "MC_SCERecomb2" || Syst == "MC_XSec" || Syst == "MC_G4" || Syst == "MC_Flux" ||
						Syst == "SmEff_LY" || Syst == "SmEff_TPC" || Syst == "SmEff_SCERecomb2" || Syst == "SmEff_XSec" || Syst == "SmEff_G4" || Syst == "SmEff_Flux" ||
						Syst == "MC_Stat" || Syst == "NuWro" 
					) {						

						for (int alt = 0; alt < NAltModels; alt++ ) {
						
							//----------------------------------------//

							if ( (Syst == "LY" || Syst == "MC_LY" || Syst == "SmEff_LY") && Runs[WhichRun] == "Run1" && AltModels[alt] == "_LYAttenuation") 
								{ continue;}

							if ( (Syst == "SCERecomb2" || Syst == "MC_SCERecomb2" || Syst == "SmEff_SCERecomb2") 
								&& Runs[WhichRun] == "Run1" && AltModels[alt] == "_CVextra") 
								{ continue;}	
								
							//----------------------------------------//	
							
							// Total covariance matrix							
						
							double CurrentFracCovEntry = FracCovariances[WhichRun][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1);
							double CurrentFracCovError = FracCovariances[WhichRun][WhichPlot]->GetBinError(WhichXBin+1,WhichYBin+1);

							double CurrentCovEntry = Covariances[WhichRun][WhichPlot]->GetBinContent(WhichXBin+1,WhichYBin+1);
							double CurrentCovError = Covariances[WhichRun][WhichPlot]->GetBinError(WhichXBin+1,WhichYBin+1);

							AltDataEntryX = AltBeamOnPlots[WhichPlot][alt]->GetBinContent(WhichXBin+1);
							AltDataErrorX = AltBeamOnPlots[WhichPlot][alt]->GetBinError(WhichXBin+1);

							AltDataEntryY = AltBeamOnPlots[WhichPlot][alt]->GetBinContent(WhichYBin+1);
							AltDataErrorY = AltBeamOnPlots[WhichPlot][alt]->GetBinError(WhichYBin+1);

							double LocalFracCovEntry = TMath::Max( ((AltDataEntryX - DataEntryX) / DataEntryX) * ( (AltDataEntryY - DataEntryY) / DataEntryY),1E-8);
							double LocalCovEntry = TMath::Max( (AltDataEntryX - DataEntryX) * (AltDataEntryY - DataEntryY),1E-8);

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

							CovFracError = 1E-8;
							CovError = 1E-8;

							FracCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovFracEntry);
							FracCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovFracError);

							Covariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovEntry);
							Covariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovError);
							
							//----------------------------------------//

						}

					}

					else { cout << endl << "Don't know how to handle this systematic uncentainty! Abort !"<< endl << endl; return; }

					//----------------------------------------//

					// Setting the elements of the Cov Matrix

					FracCovariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovFracEntry);
					FracCovariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovFracError);

					Covariances[WhichRun][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovEntry);
					Covariances[WhichRun][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovError);

					//----------------------------------------//

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
