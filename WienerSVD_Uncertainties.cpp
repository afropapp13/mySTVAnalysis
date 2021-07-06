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

#include <iomanip>
#include <iostream>

#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "ubana/AnalysisCode/Secondary_Code/GlobalSettings.cpp"
#include "ubana/AnalysisCode/Secondary_Code/myFunctions.cpp"

// TString Syst = "NTarget" "POT" "Stat" "LY" "TPC" "XSec" "G4" "Flux" "Dirt" "MC_Stat" "MC_POT" "MC_NTarget" "MC_LY" "MC_TPC" "MC_XSec" "MC_G4" "MC_Flux" "MC_Dirt"

void WienerSVD_Uncertainties(TString Syst = "None",TString BaseMC = "Overlay9",TString BeamOnSample = "BeamOn9",TString BeamOffSample = "ExtBNB9",TString DirtSample = "OverlayDirt9") {

	// -------------------------------------------------------------------------------------

	if (Syst == "None") { cout << "What the heck are you doing ? Specify the systematic uncertainty that you want to obtain!" << endl; return ;}

	// -------------------------------------------------------------------------------------

	GlobalSettings();
	TGaxis::SetMaxDigits(3);
	std::cout << std::fixed << std::setprecision(2);

	// -------------------------------------------------------------------------------------

	vector<TString> PlotNames;
//	PlotNames.push_back("DeltaPTPlot"); 
//	PlotNames.push_back("DeltaAlphaTPlot"); 
//	PlotNames.push_back("DeltaPhiTPlot"); 
//	PlotNames.push_back("MuonMomentumPlot"); 
//	PlotNames.push_back("MuonPhiPlot"); 
	PlotNames.push_back("MuonCosThetaPlot");
//	PlotNames.push_back("ProtonMomentumPlot"); 
//	PlotNames.push_back("ProtonPhiPlot"); 
//	PlotNames.push_back("ProtonCosThetaPlot");

	const int N1DPlots = PlotNames.size();
		
	// -------------------------------------------------------------------------------------------------------------------------------------

	//vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
	//Runs.push_back("Run3");
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

	if (Syst == "LY" || Syst == "MC_LY") {

		AltModels.push_back("_LYDown"); Colors.push_back(kRed+1); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_LYRayleigh"); Colors.push_back(kGreen+2); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_LYAttenuation"); Colors.push_back(kOrange+1); Universes.push_back(1); AltUniverses.push_back(1);

	}

	if (Syst == "TPC" || Syst == "MC_TPC") {

		AltModels.push_back("_WireModX"); Colors.push_back(kRed+1); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_WireModYZ"); Colors.push_back(kGreen+2); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_WireModThetaXZ"); Colors.push_back(kOrange+1); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_WireModThetaYZ"); Colors.push_back(kBlue-3); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_SCE"); Colors.push_back(kBlue); Universes.push_back(1); AltUniverses.push_back(1);
		AltModels.push_back("_Recombination2"); Colors.push_back(kMagenta); Universes.push_back(1); AltUniverses.push_back(1);
		//AltModels.push_back("_dEdx"); Colors.push_back(kYellow+2); Universes.push_back(1); AltUniverses.push_back(1);

	}

	if (Syst == "XSec" || Syst == "MC_XSec") {

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

	if (Syst == "G4" || Syst == "MC_G4") {

		UniAltModels.push_back("_reinteractions_piminus_Geant4"); Universes.push_back(100);
		UniAltModels.push_back("_reinteractions_piplus_Geant4"); Universes.push_back(100);
		UniAltModels.push_back("_reinteractions_proton_Geant4"); Universes.push_back(100);

		for (int UniAlt = 0; UniAlt < (int)(UniAltModels.size()); UniAlt++ ) {

			for (int Uni = 0; Uni < Universes[UniAlt]; Uni++ ) {

				AltModels.push_back(UniAltModels[UniAlt]+"_"+TString(std::to_string(Uni))); Colors.push_back(kGreen+2);
				AltUniverses.push_back(Universes[UniAlt]);
			
			}

		}

	}

	if (Syst == "Flux" || Syst == "MC_Flux") {

		UniAltModels.push_back("_horncurrent_FluxUnisim"); Universes.push_back(100);
		UniAltModels.push_back("_kminus_PrimaryHadronNormalization"); Universes.push_back(100);
		UniAltModels.push_back("_kplus_PrimaryHadronFeynmanScaling"); Universes.push_back(100);
		UniAltModels.push_back("_kzero_PrimaryHadronSanfordWang"); Universes.push_back(100);
		UniAltModels.push_back("_nucleoninexsec_FluxUnisim"); Universes.push_back(100);
		UniAltModels.push_back("_nucleonqexsec_FluxUnisim"); Universes.push_back(100);
		UniAltModels.push_back("_nucleontotxsec_FluxUnisim"); Universes.push_back(100);
		UniAltModels.push_back("_piminus_PrimaryHadronSWCentralSplineVariation"); Universes.push_back(100);
		UniAltModels.push_back("_pioninexsec_FluxUnisim"); Universes.push_back(100);
		UniAltModels.push_back("_pionqexsec_FluxUnisim"); Universes.push_back(100);
		UniAltModels.push_back("_piontotxsec_FluxUnisim"); Universes.push_back(100);
		UniAltModels.push_back("_piplus_PrimaryHadronSWCentralSplineVariation"); Universes.push_back(100);
		UniAltModels.push_back("_expskin_FluxUnisim"); Universes.push_back(10);

		int NUniAltModels = (int)(UniAltModels.size());

		for (int UniAlt = 0; UniAlt < NUniAltModels; UniAlt++ ) {

			for (int Uni = 0; Uni < Universes[UniAlt]; Uni++ ) {

				AltModels.push_back(UniAltModels[UniAlt]+"_"+TString(std::to_string(Uni))); Colors.push_back(kGreen+2);
				AltUniverses.push_back(Universes[UniAlt]);

				TString DublicateOverlaySample = UniAltModels[UniAlt];
				TString ReducedOverlaySample = DublicateOverlaySample.ReplaceAll("m_","m");
				if ( !(string(UniAltModels[UniAlt]).find("expskin") != std::string::npos) ) { ReducedOverlaySample = ReducedOverlaySample.ReplaceAll("n_","n"); }
				ReducedOverlaySample = ReducedOverlaySample.ReplaceAll("g_","g");				
				
				for (int i = 0; i < 10;i++) { ReducedOverlaySample.ReplaceAll(TString(std::to_string(i)),""); }
				
				TString FluxHistoName = "numu_ms"+ReducedOverlaySample+"/hEnumu"+ReducedOverlaySample+"_ms_"+TString(std::to_string(Uni));
				UniHistoFlux.push_back( (TH1D*)(FluxFile->Get(FluxHistoName) ) );
	
			}

		}

	}

	// -------------------------------------------------------------------------------------------------------------------------------------

	int NAltModels = AltModels.size();

	vector< vector<TFile*> > AltMCFileSample; AltMCFileSample.resize(NRuns,vector<TFile*>(NAltModels));

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	// Alternative Plots

	vector <vector <TH1D*> > AltCC1pPlots; AltCC1pPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));
	vector <vector <TH1D*> > AltNonCC1pPlots; AltNonCC1pPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));
	vector <vector <TH1D*> > AltBeamOnPlots; AltBeamOnPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));
	vector <vector <TH1D*> > AltBeamOffPlots; AltBeamOffPlots.resize(N1DPlots,vector<TH1D*>(NAltModels));

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------
	
		double DataPOT = ReturnBeamOnRunPOT(Runs[WhichRun]);						
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface) * (SoftFidSurface / Nominal_UB_XY_Surface);

		if (Syst == "Flux" || Syst == "MC_Flux") {

			UniFlux.clear();

			for (int UniAlt = 0; UniAlt < NAltModels; UniAlt++ ) {

				double UniFluxPOT = (UniHistoFlux[UniAlt]->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface) * (SoftFidSurface / Nominal_UB_XY_Surface);
				UniFlux.push_back(UniFluxPOT);

			}

		}

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------
		
		// Open base files

		TString ExactFileLocation = PathToFiles+CutExtension;
		TString TStringBaseMC = ExactFileLocation+"/STVStudies_"+BaseMC+"_"+Runs[WhichRun]+CutExtension+".root";
		if (Syst == "LY" || Syst == "TPC" || Syst == "MC_LY" || Syst == "MC_TPC") { TStringBaseMC = ExactFileLocation+"/STVStudies_"+BaseMC+"_"+Runs[WhichRun]+"_CV"+CutExtension+".root"; }

		MCFileSample[WhichRun] = TFile::Open(TStringBaseMC,"readonly");
		BeamOnFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOnSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");
		BeamOffFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOffSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");
		DirtFileSample[WhichRun] = TFile::Open(ExactFileLocation+"/STVStudies_"+DirtSample+"_"+Runs[WhichRun]+CutExtension+".root","readonly");

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

				if ( (Syst == "LY" || Syst == "MC_LY") && Runs[WhichRun] == "Run1" && AltModels[alt] == "_LYAttenuation") 
					{ continue;}

				// Open Alternative MC files

				TString TStringAltBaseMC = ExactFileLocation+"/STVStudies_"+BaseMC+"_"+Runs[WhichRun]+AltModels[alt]+CutExtension+".root";
				AltMCFileSample[WhichRun][alt] = TFile::Open(TStringAltBaseMC,"readonly");

				AltCC1pPlots[WhichPlot][alt] = (TH1D*)(AltMCFileSample[WhichRun][alt]->Get("CC1pReco"+PlotNames[WhichPlot]));
				AltNonCC1pPlots[WhichPlot][alt] = (TH1D*)(AltMCFileSample[WhichRun][alt]->Get("NonCC1pReco"+PlotNames[WhichPlot]));

				AltCC1pPlots[WhichPlot][alt]->SetDirectory(0); // to decouple it from the open file directory
				AltNonCC1pPlots[WhichPlot][alt]->SetDirectory(0); // to decouple it from the open file directory

				AltMCFileSample[WhichRun][alt]->Close();

			}

			// -------------------------------------------------------------------------------------------------------

			TH1D* DataPlot = (TH1D*)BeamOnPlots[WhichPlot]->Clone();
			TH1D* AltDataPlot = (TH1D*)BeamOnPlots[WhichPlot]->Clone();

			if (BeamOnSample == "BeamOn9" || BeamOnSample == "Overlay9") { // Beam On Sample, need to subtract ALL backgrounds

				DataPlot->Add(NonCC1pPlots[WhichPlot],-1.);
				DataPlot->Add(BeamOffPlots[WhichPlot],-1.);
				DataPlot->Add(DirtPlots[WhichPlot],-1.);

				for (int alt = 0; alt < NAltModels; alt++ ) {

					if ( (Syst == "LY" || Syst == "MC_LY") && Runs[WhichRun] == "Run1" && AltModels[alt] == "_LYAttenuation") 
						{ continue;}

					AltBeamOnPlots[WhichPlot][alt] = (TH1D*)BeamOnPlots[WhichPlot]->Clone();
					AltBeamOnPlots[WhichPlot][alt]->Add(AltNonCC1pPlots[WhichPlot][alt],-1.);
					AltBeamOnPlots[WhichPlot][alt]->Add(BeamOffPlots[WhichPlot],-1.);
					AltBeamOnPlots[WhichPlot][alt]->Add(DirtPlots[WhichPlot],-1.);

				}

				if (Syst == "Dirt" || Syst == "MC_Dirt") {

					AltDataPlot->Add(NonCC1pPlots[WhichPlot],-1.);
					AltDataPlot->Add(BeamOffPlots[WhichPlot],-1.);

					// Scale the dirt sample by 25% down
					TH1D* DirtClone = (TH1D*)(DirtPlots[WhichPlot]->Clone());
					DirtClone->Scale(0.75);
					AltDataPlot->Add(DirtClone,-1.);

				}

			} else { // Fake Data Studies, working with MC CC1p signal events as BeamOn

				DataPlot = CC1pPlots[WhichPlot];

				for (int alt = 0; alt < NAltModels; alt++ ) {

					if ( (Syst == "LY" || Syst == "MC_LY") && Runs[WhichRun] == "Run1" && AltModels[alt] == "_LYAttenuation") 
						{ continue;}

					AltBeamOnPlots[WhichPlot][alt] = AltCC1pPlots[WhichPlot][alt];

				}

			}	

			// -------------------------------------------------------------------------------------------------------
			// -------------------------------------------------------------------------------------------------------
			
			int NBins = BeamOnPlots[WhichPlot]->GetXaxis()->GetNbins();

			// -------------------------------------------------------------------------------------------------------

			cout << endl;

			DataPlot->Scale( 1. / (IntegratedFlux * NTargets) * Units );
			double DataIntXSec = DataPlot->Integral("width");

			double Total = 0.;

			// -------------------------------------------------------------------------------------------------------
			// -------------------------------------------------------------------------------------------------------

			// Based on the type of systematic unceratainty, choose how to handle it

			if (Syst == "NTarget" || Syst == "MC_NTarget") {

				AltDataPlot = (TH1D*)(DataPlot->Clone());
				AltDataPlot->Scale( (1+NTargetUncertainty) );
				double AltDataIntXSec = AltDataPlot->Integral("width");

				double FracContr = (AltDataIntXSec - DataIntXSec) / DataIntXSec * 100.;  
				Total = FracContr;

				cout << Syst << " & " << FracContr << " \\\\" << endl;

			} else if (Syst == "POT" || Syst == "MC_POT") {

				AltDataPlot = (TH1D*)(DataPlot->Clone());
				AltDataPlot->Scale( (1+POTUncertainty) );
				double AltDataIntXSec = AltDataPlot->Integral("width");

				double FracContr = (AltDataIntXSec - DataIntXSec) / DataIntXSec * 100.;  
				Total = FracContr;

				cout << Syst << " & " << FracContr << " \\\\" << endl;

			} else if (Syst == "Dirt" || Syst == "MC_Dirt") {

				AltDataPlot->Scale( 1. / (IntegratedFlux * NTargets) * Units );
				double AltDataIntXSec = AltDataPlot->Integral("width");

				double FracContr = (AltDataIntXSec - DataIntXSec) / DataIntXSec * 100.;  
				Total = FracContr;

				cout << Syst << " & " << FracContr << " \\\\" << endl;

			} else if (Syst == "Stat" || Syst == "MC_Stat") {

				double StatUnc = 0.;

				for (int WhichBin = 0; WhichBin < NBins; WhichBin++) {

					double BinError = DataPlot->GetBinError(WhichBin+1);
					double BinWidth = DataPlot->GetBinWidth(WhichBin+1);

					StatUnc += BinError * BinWidth;

				}

				double FracContr = StatUnc / DataIntXSec * 100.;
				Total = FracContr;

				cout << Syst << " & " << FracContr << " \\\\" << endl;

			} else if (
				Syst == "LY" || Syst == "TPC" || Syst == "XSec" || Syst == "G4" || Syst == "Flux" ||
				Syst == "MC_LY" || Syst == "MC_TPC" || Syst == "MC_XSec" || Syst == "MC_G4" || Syst == "MC_Flux"
			) {						

				for (int alt = 0; alt < NAltModels; alt++ ) {

					if ( (Syst == "LY" || Syst == "MC_LY") && Runs[WhichRun] == "Run1" && AltModels[alt] == "_LYAttenuation") 
						{ continue;}	

					if (Syst == "Flux" || Syst == "MC_Flux") {
						AltBeamOnPlots[WhichPlot][alt]->Scale( 1. / (UniFlux[alt] * NTargets) * Units );
					} else {
						AltBeamOnPlots[WhichPlot][alt]->Scale( 1. / (IntegratedFlux * NTargets) * Units );
					}
			
					double AltDataIntXSec = AltBeamOnPlots[WhichPlot][alt]->Integral("width");

					double FracContr = (AltDataIntXSec - DataIntXSec) / DataIntXSec * 100.;  
					if (AltUniverses[alt] > 2) { FracContr = FracContr * 1. / TMath::Sqrt(AltUniverses[alt]); }

					Total = TMath::Sqrt( TMath::Power(Total,2.) + TMath::Power(FracContr,2.) );

					TString copy = AltModels[alt];
					TString label = copy.ReplaceAll("_","");

					cout << label << " & " << FracContr << " \\\\" << endl;						

				}

			}

			else { cout << "Don't know how to handle this systematic uncentainty! Abort !"<< endl << endl; return; }

			// -------------------------------------------------------------------------------------------------------

			cout << "\\hline \\hline" << endl;
			cout << Syst << " & " << Total << " \\\\" << endl;		

		} // End of the loop over the plots

		// ---------------------------------------------------------------------------------------	

		// Close base files

		MCFileSample[WhichRun]->Close();
		BeamOnFileSample[WhichRun]->Close();
		BeamOffFileSample[WhichRun]->Close();
		DirtFileSample[WhichRun]->Close();

		// ---------------------------------------------------------------------------------------	

		cout << endl; 

	} // End of the loop over the runs	

} // End of the program
