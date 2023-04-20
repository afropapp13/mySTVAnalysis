#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TEfficiency.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLine.h>
#include <TGaxis.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "ubana/myClasses/Constants.h"
#include "ubana/myClasses/Util.h"
#include "ubana/myClasses/WienerSVD.h"

using namespace std;
using namespace Constants;

#include "ubana/AnalysisCode/Secondary_Code/GlobalSettings.cpp"
#include "ubana/AnalysisCode/Secondary_Code/myFunctions.cpp"

// -----------------------------------------------------------------------------------------------

void WienerSVD_Merge_EventRateCovariances(TString OverlaySample = "Overlay9", TString BeamOn9 = "", TString Tune = "") {

	// -------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();

	gStyle->SetOptStat(0);	

	// -------------------------------------------------------------------------------------

	GlobalSettings();
	TGaxis::SetMaxDigits(3);			

	// -------------------------------------------------------------------------------------

	vector<TString> PlotNames;
	PlotNames.push_back("LLRPIDPlot");
	PlotNames.push_back("MuonLLRPIDPlot");
	PlotNames.push_back("ProtonLLRPIDPlot");

	const int NPlots = PlotNames.size();

	// -----------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;			
	Runs.push_back("Combined");				

	int NRuns = (int)(Runs.size());

	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> UncSources;

	if (BeamOn9 != "") { // Fake data study, we need only the Stat/MC_Stat & XSec covariances

		UncSources.push_back("Stat");
		UncSources.push_back("XSec");
		UncSources.push_back("MC_Stat");
		UncSources.push_back("NuWro");			

	} else {

		UncSources.push_back("Stat");
		UncSources.push_back("LY");
		UncSources.push_back("TPC");
		UncSources.push_back("SCERecomb2");
		UncSources.push_back("XSec");
		UncSources.push_back("G4");
		UncSources.push_back("Flux");
		UncSources.push_back("Dirt");
		UncSources.push_back("POT"); 
		UncSources.push_back("NTarget");
		UncSources.push_back("MC_Stat");
		UncSources.push_back("NuWro");		

	}

	int NSamples = UncSources.size();

	vector<TFile*> CovFiles;
	CovFiles.resize(NSamples);

	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TH2D*> FracCovariances; FracCovariances.resize(NPlots);
	vector<TH2D*> FracStatCovariances; FracStatCovariances.resize(NPlots);
	vector<TH2D*> FracSystCovariances; FracSystCovariances.resize(NPlots);

	vector<TH2D*> Covariances; Covariances.resize(NPlots);
	vector<TH2D*> StatCovariances; StatCovariances.resize(NPlots);
	vector<TH2D*> MCStatCovariances; MCStatCovariances.resize(NPlots);	
	vector<TH2D*> SystCovariances; SystCovariances.resize(NPlots);	
	vector<TH2D*> LYCovariances; LYCovariances.resize(NPlots);	
	vector<TH2D*> TPCCovariances; TPCCovariances.resize(NPlots);
	vector<TH2D*> SCERecomb2Covariances; SCERecomb2Covariances.resize(NPlots);
	vector<TH2D*> XSecCovariances; XSecCovariances.resize(NPlots);
	vector<TH2D*> G4Covariances; G4Covariances.resize(NPlots);
	vector<TH2D*> FluxCovariances; FluxCovariances.resize(NPlots);
	vector<TH2D*> DirtCovariances; DirtCovariances.resize(NPlots);	
	vector<TH2D*> POTCovariances; POTCovariances.resize(NPlots);
	vector<TH2D*> NuWroCovariances; NuWroCovariances.resize(NPlots);
	vector<TH2D*> NTargetCovariances; NTargetCovariances.resize(NPlots);							

	// -------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {					

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString TotalFileCovarianceSpecName = Tune + BeamOn9 + "WienerSVD_Total_EventRateCovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		if (BeamOn9 != "") { TotalFileCovarianceSpecName = Tune + "WienerSVD_Total_EventRateCovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }

		TString TotalFileCovarianceName = MigrationMatrixPath + TotalFileCovarianceSpecName;
		TFile* TotalFileCovarianceMatrices = new TFile(TotalFileCovarianceName,"recreate");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

			TCanvas* MCERPlotCanvas = nullptr;
			TLegend* legMC = nullptr;
			TString MCERCanvasName = Tune + "MCERSyst_"+PlotNames[WhichPlot]+OverlaySample+"_"+Runs[WhichRun];

			// ----------------------------------------------------------------------------------------------------

			// Loop over the samples / sources of uncertainty

			int Counter = 0; // Counter for colors

			TH2D* DetFracCovMatrix = nullptr; // matrix for all detector variations

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				// -----------------------------------------------------------------------------------------------------------------------------------------

				// Opening the file containing the covariance matrices for each one of the systematics

//				TString FileCovarianceSpecName = BeamOn9+"WienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
				TString FileCovarianceSpecName = "WienerSVD_" + UncSources[WhichSample] + "_EventRateCovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
				if (BeamOn9 != "" && UncSources[WhichSample] == "Stat") { FileCovarianceSpecName = Tune + "WienerSVD_" + UncSources[WhichSample] + "_EventRateCovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }
				if (BeamOn9 != "" && UncSources[WhichSample] == "NuWro" && ( Tune == "GENIEv2" || Tune == "NoTune" || Tune == "TwiceMEC" ) ) 
					{ FileCovarianceSpecName = "WienerSVD_" + UncSources[WhichSample] + "_EventRateCovarianceMatrices_"+OverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }


				// For the detector variation, we follow the PeLEE recipe
				// Only Run3 and propagate across all runs

				if (UncSources[WhichSample] == "LY" || UncSources[WhichSample] == "TPC"  || UncSources[WhichSample] == "SCERecomb2" || UncSources[WhichSample] == "MC_LY" 
				|| UncSources[WhichSample] == "MC_TPC" || UncSources[WhichSample] == "MC_SCERecomb2" || UncSources[WhichSample] == "SmEff_LY" 
				|| UncSources[WhichSample] == "SmEff_TPC" || UncSources[WhichSample] == "SmEff_SCERecomb2" ) {

//					FileCovarianceSpecName = BeamOn9+"WienerSVD_" + UncSources[WhichSample] + "_CovarianceMatrices_"+OverlaySample+"_Run3_"+UBCodeVersion+".root";
					FileCovarianceSpecName = Tune + "WienerSVD_" + UncSources[WhichSample] + "_EventRateCovarianceMatrices_"+OverlaySample+"_Run3_"+UBCodeVersion+".root";

				}

				TString FileCovarianceName = MigrationMatrixPath + FileCovarianceSpecName;
				CovFiles[WhichSample] = new TFile(FileCovarianceName,"readonly");

				// -----------------------------------------------------------------------------------------------------------------------------------------

				TString LocalCovMatrixName = UncSources[WhichSample]+"_Covariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun];
				TString LocalFracCovMatrixName = UncSources[WhichSample]+"_FracCovariance_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun];			

				// For the detector variation, we follow the PeLEE recipe
				// Only Run3 and propagate across all runs

				if (UncSources[WhichSample] == "LY" || UncSources[WhichSample] == "TPC" || UncSources[WhichSample] == "SCERecomb2" || UncSources[WhichSample] == "MC_LY" 
				|| UncSources[WhichSample] == "MC_TPC" || UncSources[WhichSample] == "MC_SCERecomb2" || UncSources[WhichSample] == "SmEff_LY" 
				|| UncSources[WhichSample] == "SmEff_TPC" || UncSources[WhichSample] == "SmEff_SCERecomb2" ) {

					LocalCovMatrixName = UncSources[WhichSample]+"_Covariance_"+PlotNames[WhichPlot]+"_Run3";
					LocalFracCovMatrixName = UncSources[WhichSample]+"_FracCovariance_"+PlotNames[WhichPlot]+"_Run3";					

				}

				TH2D* LocalCovMatrix = (TH2D*)( CovFiles[WhichSample]->Get(LocalCovMatrixName) );
				TH2D* LocalCovMatrixClone = (TH2D*)(LocalCovMatrix->Clone() );
				TH2D* LocalFracCovMatrix = (TH2D*)( CovFiles[WhichSample]->Get(LocalFracCovMatrixName) );			

				LocalCovMatrix->SetDirectory(0);
				LocalFracCovMatrix->SetDirectory(0);				
				LocalCovMatrixClone->SetDirectory(0);				
				CovFiles[WhichSample]->Close();

				// -------------------------------------------------------------------------------------------------

				// Stat
				if (WhichSample == 0) { 

					StatCovariances[WhichPlot] = LocalCovMatrix; 
					FracStatCovariances[WhichPlot] = LocalFracCovMatrix;					
				
				} else {

					if (WhichSample == 1) { 

						if ( StatCovariances[WhichPlot]->GetNbinsX() != LocalCovMatrixClone->GetNbinsX() )
							{ cout << UncSources[WhichSample] << " " << PlotNames[WhichPlot] << " covariance matrix with different number of bins" << endl; }
						
						SystCovariances[WhichPlot] = LocalCovMatrixClone; 
						FracSystCovariances[WhichPlot] = LocalFracCovMatrix;

					} else { 

						if ( StatCovariances[WhichPlot]->GetNbinsX() != LocalCovMatrixClone->GetNbinsX() )
							{ cout << UncSources[WhichSample] << " " << PlotNames[WhichPlot] << " covariance matrix with different number of bins" << endl; }						
						
						SystCovariances[WhichPlot]->Add(LocalCovMatrixClone); 
						FracSystCovariances[WhichPlot]->Add(LocalFracCovMatrix);						 
						
					}					

					if (UncSources[WhichSample] == "MC_Stat") { MCStatCovariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "LY") { LYCovariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "TPC") { TPCCovariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "SCERecomb2") { SCERecomb2Covariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "XSec") { XSecCovariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "G4") { G4Covariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "Flux") { FluxCovariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "Dirt") { DirtCovariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "POT") { POTCovariances[WhichPlot] = LocalCovMatrix; }
					if (UncSources[WhichSample] == "NuWro") { NuWroCovariances[WhichPlot] = LocalCovMatrix; }					
					if (UncSources[WhichSample] == "NTarget") { NTargetCovariances[WhichPlot] = LocalCovMatrix; }																				

				}

			} // End of the loop over the samples

			// ------------------------------------------------------------------			

			// Storing Stat/Syst/Total Covariance Matrices in root file

			TotalFileCovarianceMatrices->cd();

			MCStatCovariances[WhichPlot]->Write("MCStatCovariance_"+PlotNames[WhichPlot]);

			if (BeamOn9 == "") {

				LYCovariances[WhichPlot]->Write("LYCovariance_"+PlotNames[WhichPlot]);
				TPCCovariances[WhichPlot]->Write("TPCCovariance_"+PlotNames[WhichPlot]);			
				SCERecomb2Covariances[WhichPlot]->Write("SCERecomb2Covariance_"+PlotNames[WhichPlot]);			
				XSecCovariances[WhichPlot]->Write("XSecCovariance_"+PlotNames[WhichPlot]);
				G4Covariances[WhichPlot]->Write("G4Covariance_"+PlotNames[WhichPlot]);			
				FluxCovariances[WhichPlot]->Write("FluxCovariance_"+PlotNames[WhichPlot]);
				DirtCovariances[WhichPlot]->Write("DirtCovariance_"+PlotNames[WhichPlot]);
				POTCovariances[WhichPlot]->Write("POTCovariance_"+PlotNames[WhichPlot]);
				NuWroCovariances[WhichPlot]->Write("NuWroCovariance_"+PlotNames[WhichPlot]);							
				NTargetCovariances[WhichPlot]->Write("NTargetCovariance_"+PlotNames[WhichPlot]);

			}			

			StatCovariances[WhichPlot]->Write("StatCovariance_"+PlotNames[WhichPlot]);
			FracStatCovariances[WhichPlot]->Write("FracStatCovariance_"+PlotNames[WhichPlot]);		

			SystCovariances[WhichPlot]->Write("SystCovariance_"+PlotNames[WhichPlot]);
			FracSystCovariances[WhichPlot]->Write("FracSystCovariance_"+PlotNames[WhichPlot]);		

			TH2D* CloneCovariances = (TH2D*)(StatCovariances[WhichPlot]->Clone());
			if ( SystCovariances[WhichPlot]->GetNbinsX() != CloneCovariances->GetNbinsX() )
				{ cout << PlotNames[WhichPlot] << " covariance matrix with different number of bins" << endl; }			
			CloneCovariances->Add(SystCovariances[WhichPlot]);
			CloneCovariances->Write("TotalCovariance_"+PlotNames[WhichPlot]);

			TH2D* CloneFracCovariances = (TH2D*)(FracStatCovariances[WhichPlot]->Clone());
			if ( FracSystCovariances[WhichPlot]->GetNbinsX() != CloneFracCovariances->GetNbinsX() )
				{ cout << PlotNames[WhichPlot] << " covariance matrix with different number of bins" << endl; }			
			CloneFracCovariances->Add(FracSystCovariances[WhichPlot]);
			CloneFracCovariances->Write("FracTotalCovariance_"+PlotNames[WhichPlot]);		
		
		} // End of the loop over the plots

		TotalFileCovarianceMatrices->Close();
		cout << endl << "Covariance matrix file " << TotalFileCovarianceName << " has been created" << endl << endl;

		cout << "Merging of covariance matrices for run " << Runs[WhichRun] << " completed!" << endl;

	} // End of the loop over the runs	

} // End of the program 