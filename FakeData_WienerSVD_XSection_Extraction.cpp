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

#include "ubana/AnalysisCode/Secondary_Code/myFunctions.cpp"

// -------------------------------------------------------------------------------------------------------------------------------------

TString ToStringPOT(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

// -------------------------------------------------------------------------------------------------------------------------------------

void ReweightXSec(TH1D* h, double SF = 1.) {

	int NBins = h->GetXaxis()->GetNbins();

	for (int i = 0; i < NBins; i++) {

		double CurrentEntry = h->GetBinContent(i+1);
		double NewEntry = CurrentEntry * SF / h->GetBinWidth(i+1);

		double CurrentError = h->GetBinError(i+1);
		double NewError = CurrentError * SF / h->GetBinWidth(i+1);

		h->SetBinContent(i+1,NewEntry); 
		h->SetBinError(i+1,NewError); 
//		h->SetBinError(i+1,0.000001); 

	}

}

// -------------------------------------------------------------------------------------------------------------------------------------

void FakeData_WienerSVD_XSection_Extraction(TString OverlaySample = "Overlay9", TString BeamOnSample = "BeamOn9") {

	// -------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();	
	gStyle->SetOptStat(0);
	gStyle->SetEndErrorSize(4);	

	TString Subtract = "";

	int DecimalAccuracy = 2;

	// -------------------------------------------------------------------------------------

	int NEventsPassingSelectionCuts = 0;
	TString CutExtension = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	// v52
	VectorCuts.push_back("");
	VectorCuts.push_back("_PID");
	VectorCuts.push_back("_NuScore");

	int NCuts = (int)(VectorCuts.size());	

	for (int i = 0; i < NCuts; i++) { CutExtension = CutExtension + VectorCuts[i]; }

	// -------------------------------------------------------------------------------------

//	vector<TString> PlotNames;
//	PlotNames.push_back("DeltaPTPlot"); 
//	PlotNames.push_back("DeltaAlphaTPlot"); 
//	PlotNames.push_back("DeltaPhiTPlot");
//	PlotNames.push_back("MuonMomentumPlot"); 
//	PlotNames.push_back("MuonCosThetaPlot"); 
//	PlotNames.push_back("MuonPhiPlot");
//	PlotNames.push_back("ProtonMomentumPlot"); 
//	PlotNames.push_back("ProtonCosThetaPlot");
//	PlotNames.push_back("ProtonPhiPlot");
//	PlotNames.push_back("ECalPlot");
//	PlotNames.push_back("EQEPlot"); 
//	PlotNames.push_back("Q2Plot");

//	PlotNames.push_back("CCQEMuonMomentumPlot"); 
//	PlotNames.push_back("CCQEMuonCosThetaPlot"); 
//	PlotNames.push_back("CCQEProtonMomentumPlot"); 
//	PlotNames.push_back("CCQEProtonCosThetaPlot");

	const int N1DPlots = PlotNames.size();
	//cout << "Number of 1D Plots = " << N1DPlots << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");

	vector<TString> LabelsOfSamples;
	vector<TString> NameOfSamples;

	NameOfSamples.push_back(OverlaySample);
	NameOfSamples.push_back(BeamOnSample);
	NameOfSamples.push_back("ExtBNB9"); 
	NameOfSamples.push_back("OverlayDirt9");
	NameOfSamples.push_back("GenieOverlay");	
	NameOfSamples.push_back("AltEventGen");	 

	//----------------------------------------//

	vector<TString> Runs;
//	Runs.push_back("Run1");
////	Runs.push_back("Run2");
//	Runs.push_back("Run3");
////	Runs.push_back("Run4");
////	Runs.push_back("Run5");	
	Runs.push_back("Combined");			

	int NRuns = (int)(Runs.size());
	//cout << "Number of Runs = " << NRuns << endl;

	//----------------------------------------//

	// CV Flux File

	TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));

	//----------------------------------------//	

	// open the file that contains the unfolding uncertainty

	//TString NameUnfUnc = PathToExtractedXSec+"WienerSVD_UnfoldingUnc_Combined_"+UBCodeVersion+Subtract+".root";
	//TFile* FileUnfUnc = TFile::Open(NameUnfUnc,"readonly");			

	//----------------------------------------//

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		//----------------------------------------//
	
		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);					
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface);	
			
		//----------------------------------------//

		vector<TCanvas*> PlotCanvas; PlotCanvas.clear();

		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
		vector<vector<TH1D*> > QEPlotsTrue; QEPlotsTrue.clear();
		vector<vector<TH1D*> > MECPlotsTrue; MECPlotsTrue.clear();
		vector<vector<TH1D*> > RESPlotsTrue; RESPlotsTrue.clear();
		vector<vector<TH1D*> > DISPlotsTrue; DISPlotsTrue.clear();
		vector<vector<TH1D*> > COHPlotsTrue; COHPlotsTrue.clear();		
		vector<vector<TH1D*> > PlotsBkgReco; PlotsBkgReco.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();

		//vector<TH1D*> UnfUnc; UnfUnc.clear();
		vector<TH2D*> ResponseMatrices; ResponseMatrices.clear();
		vector<TH2D*> CovarianceMatrices; CovarianceMatrices.clear();
		vector<TH2D*> MCStatCovarianceMatrices; MCStatCovarianceMatrices.clear();		

		//----------------------------------------//

		TString FileResponseName = MigrationMatrixPath+"FileResponseMatrices_"+NameOfSamples[0]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";

		if (OverlaySample == "NoTuneOverlay9") { 
			FileResponseName = MigrationMatrixPath+"NoTuneFileResponseMatrices_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		}

		if (OverlaySample == "TwiceMECOverlay9") { 
			FileResponseName = MigrationMatrixPath+"TwiceMECFileResponseMatrices_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		}

		cout << "File Responses = " << FileResponseName << endl;
		TFile* FileResponseMatrices = new TFile(FileResponseName,"readonly");

		//----------------------------------------//		

		TString FileCovarianceName = MigrationMatrixPath+"WienerSVD_Total_CovarianceMatrices_"+NameOfSamples[0]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";

		if (OverlaySample == "NoTuneOverlay9") { 
			FileCovarianceName = MigrationMatrixPath+"NoTuneWienerSVD_Total_CovarianceMatrices_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; 
		}

		if (OverlaySample == "TwiceMECOverlay9") { 
			FileCovarianceName = MigrationMatrixPath+"TwiceMECWienerSVD_Total_CovarianceMatrices_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; 
		}

		// For the fake data studies with the default overlay MC and alternative fake data
		// we need only the stat, mc stat, xsec and unfolding uncertainties

		if (OverlaySample == "Overlay9" && BeamOnSample == "Overlay9NuWro") { 
			FileCovarianceName = MigrationMatrixPath+"Overlay9NuWroWienerSVD_Total_CovarianceMatrices_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; 
		}		

		if (OverlaySample == "Overlay9" && BeamOnSample == "NoTuneOverlay9") { 
			FileCovarianceName = MigrationMatrixPath+"NoTuneWienerSVD_Total_CovarianceMatrices_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; 
		}

		if (OverlaySample == "Overlay9" && BeamOnSample == "TwiceMECOverlay9") { 
			FileCovarianceName = MigrationMatrixPath+"TwiceMECWienerSVD_Total_CovarianceMatrices_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; 
		}		

		cout << "File Covariances = " << FileCovarianceName << endl;			
		TFile* FileCovarianceMatrices = new TFile(FileCovarianceName,"readonly");

		//----------------------------------------//

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();
		
		TString PathToFilesUBCodeExtension = PathToFiles+CutExtension;

		// Store the extracted xsections & associated files in dedicated file

		TString NameExtractedXSec = PathToExtractedXSec+BeamOnSample+"WienerSVD_ExtractedXSec_"+NameOfSamples[0]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+Subtract+".root";
		TFile* ExtractedXSec = TFile::Open(NameExtractedXSec,"recreate");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Loop over the samples

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			if (
			NameOfSamples[WhichSample] == "BeamOn9" || NameOfSamples[WhichSample] == "Overlay9NuWro" ||
			NameOfSamples[WhichSample] == "ExtBNB9" || 
			NameOfSamples[WhichSample] == "OverlayDirt9"
			) { 
			
			TString FileName = "STVStudies_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+CutExtension+".root";
			FileSample.push_back(TFile::Open(PathToFilesUBCodeExtension+"/"+FileName)); 
			}
			
			if (NameOfSamples[WhichSample] == "Overlay9") { 
			
			TString FileName = "STVStudies_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+CutExtension+".root";
			FileSample.push_back(TFile::Open(PathToFilesUBCodeExtension+"/"+FileName)); 
			
			}

			if (NameOfSamples[WhichSample] == "NoTuneOverlay9") { 
			
			TString FileName = "NoTuneSTVStudies_Overlay9_"+Runs[WhichRun]+CutExtension+".root";
			FileSample.push_back(TFile::Open(PathToFilesUBCodeExtension+"/"+FileName)); 
			
			}

			if (NameOfSamples[WhichSample] == "TwiceMECOverlay9") { 
			
			TString FileName = "TwiceMECSTVStudies_Overlay9_"+Runs[WhichRun]+CutExtension+".root";
			FileSample.push_back(TFile::Open(PathToFilesUBCodeExtension+"/"+FileName)); 
			
			}						

			if (NameOfSamples[WhichSample] == "GenieOverlay") { 
			
				TString FileName = "TruthSTVAnalysis_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";					
				FileSample.push_back(TFile::Open(PathToFiles+FileName));  
			
			}

			if (NameOfSamples[WhichSample] == "AltEventGen") { 
			
				TString FileName = "TruthSTVAnalysis_"+BeamOnSample+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
				if (BeamOnSample == "NoTuneOverlay9") { FileName = "NoTuneTruthSTVAnalysis_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }
				if (BeamOnSample == "TwiceMECOverlay9") { FileName = "TwiceMECTruthSTVAnalysis_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }		
				if (BeamOnSample == "Overlay9NuWro") { FileName = "TruthSTVAnalysis_Overlay9NuWro_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }	
				if (BeamOnSample == "BeamOn9" && OverlaySample == "NoTuneOverlay9") 
					{ FileName = "NoTuneTruthSTVAnalysis_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }
				if (BeamOnSample == "BeamOn9" && OverlaySample == "TwiceMECOverlay9") 
					{ FileName = "TwiceMECTruthSTVAnalysis_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; }						
						
				FileSample.push_back(TFile::Open(PathToFiles+FileName));  
			
			}			

			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();
			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();
			vector<TH1D*> QECurrentPlotsTrue; QECurrentPlotsTrue.clear();
			vector<TH1D*> MECCurrentPlotsTrue; MECCurrentPlotsTrue.clear();
			vector<TH1D*> RESCurrentPlotsTrue; RESCurrentPlotsTrue.clear();
			vector<TH1D*> DISCurrentPlotsTrue; DISCurrentPlotsTrue.clear();
			vector<TH1D*> COHCurrentPlotsTrue; COHCurrentPlotsTrue.clear();			
			vector<TH1D*> CurrentPlotsBkgReco; CurrentPlotsBkgReco.clear();
			vector<TH1D*> CurrentPlotsCC1pReco; CurrentPlotsCC1pReco.clear();

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
				CurrentPlotsReco.push_back(histReco);

				TH1D* histBkgReco = (TH1D*)(FileSample[WhichSample]->Get("NonCC1pReco"+PlotNames[WhichPlot]));
				CurrentPlotsBkgReco.push_back(histBkgReco);

				TH1D* histCC1pReco = (TH1D*)(FileSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]));
				CurrentPlotsCC1pReco.push_back(histCC1pReco);

				TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
				CurrentPlotsTrue.push_back(histTrue);

				TH1D* QEhistTrue = (TH1D*)(FileSample[WhichSample]->Get("QETrue"+PlotNames[WhichPlot]));
				QECurrentPlotsTrue.push_back(QEhistTrue);

				TH1D* MEChistTrue = (TH1D*)(FileSample[WhichSample]->Get("MECTrue"+PlotNames[WhichPlot]));
				MECCurrentPlotsTrue.push_back(MEChistTrue);	

				TH1D* REShistTrue = (TH1D*)(FileSample[WhichSample]->Get("RESTrue"+PlotNames[WhichPlot]));
				RESCurrentPlotsTrue.push_back(REShistTrue);

				TH1D* DIShistTrue = (TH1D*)(FileSample[WhichSample]->Get("DISTrue"+PlotNames[WhichPlot]));
				DISCurrentPlotsTrue.push_back(DIShistTrue);

				TH1D* COHhistTrue = (TH1D*)(FileSample[WhichSample]->Get("COHTrue"+PlotNames[WhichPlot]));
				COHCurrentPlotsTrue.push_back(COHhistTrue);					
		
			} // End of the loop over the plots

			PlotsReco.push_back(CurrentPlotsReco);		
			PlotsTrue.push_back(CurrentPlotsTrue);	
			QEPlotsTrue.push_back(QECurrentPlotsTrue);
			MECPlotsTrue.push_back(MECCurrentPlotsTrue);
			RESPlotsTrue.push_back(RESCurrentPlotsTrue);
			DISPlotsTrue.push_back(DISCurrentPlotsTrue);
			COHPlotsTrue.push_back(COHCurrentPlotsTrue);				
			PlotsBkgReco.push_back(CurrentPlotsBkgReco);
			PlotsCC1pReco.push_back(CurrentPlotsCC1pReco);

		} // End of the loop over the samples

		//----------------------------------------//

		// File to store the unfolding model uncertainties

//		TFile* fUnc = nullptr;	
//		TString fUncName = "";
//		if (BeamOnSample == "Overlay9NuWro" && OverlaySample == "Overlay9") { 
			
//			fUncName = PathToExtractedXSec+"WienerSVD_UnfoldingUnc_Combined_"+UBCodeVersion+".root";
//			fUnc = new TFile(fUncName,"recreate"); 
			
//		}	

		// ----------------------------------------------------------------------------------------------------------------------------------

		// Loop over the event rate plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

			//----------------------------------------//			

			//UnfUnc.push_back((TH1D*)FileUnfUnc->Get("UnfUnc_"+PlotNames[WhichPlot]));			

			//----------------------------------------//

			ResponseMatrices.push_back((TH2D*)FileResponseMatrices->Get("POTScaledCC1pReco"+PlotNames[WhichPlot]+"2D"));

			// Already flux-averaged rates
			CovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("TotalCovariance_"+PlotNames[WhichPlot]));
			MCStatCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("MCStatCovariance_"+PlotNames[WhichPlot]));			

			//----------------------------------------//

			// True CC1p Signal MC // No detector/ reconstruction / smearing effects

			int n = PlotsTrue[4][WhichPlot]->GetNbinsX();
			double Nuedges[n+1];
			    
			for (int i = 0; i < n+1; i++) { Nuedges[i] = PlotsTrue[4][WhichPlot]->GetBinLowEdge(i+1); }

			// -------------------------------------------------------------------------------------------------

			// BeamOn = Alternative model for fake data studies in this case
			// Thus we grab the CC1p part, without any subtractions

			TH1D* DataPlot = (TH1D*)(PlotsCC1pReco[1][WhichPlot]->Clone());
//			TH1D* DataPlot = (TH1D*)(PlotsReco[1][WhichPlot]->Clone());
//			DataPlot->Add(PlotsBkgReco[0][WhichPlot],-1);

			int m = DataPlot->GetNbinsX();			
			TString XTitle = DataPlot->GetXaxis()->GetTitle();
			TString YTitle = DataPlot->GetYaxis()->GetTitle();	

			// Flux-averaged event rates
			DataPlot->Scale(Units/(IntegratedFlux*NTargets));
			PlotsTrue[4][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets)); // CV
			QEPlotsTrue[4][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
			MECPlotsTrue[4][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
			RESPlotsTrue[4][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
			DISPlotsTrue[4][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
			COHPlotsTrue[4][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));

			PlotsTrue[5][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets)); // Alt CV
			QEPlotsTrue[5][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
			MECPlotsTrue[5][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
			RESPlotsTrue[5][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
			DISPlotsTrue[5][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
			COHPlotsTrue[5][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));			

			// -------------------------------------------------------------------------------------------

			// Construct vectors (for 1D histogram) and matrices (for 2D histogram) for input

			TVectorD signal(n);
			TVectorD QEsignal(n);
			TVectorD MECsignal(n);
			TVectorD RESsignal(n);
			TVectorD DISsignal(n);
			TVectorD COHsignal(n);			
			TVectorD altsignal(n);
			TVectorD QEaltsignal(n);
			TVectorD MECaltsignal(n);
			TVectorD RESaltsignal(n);
			TVectorD DISaltsignal(n);
			TVectorD COHaltsignal(n);						
			TVectorD measure(m);
			TMatrixD response(m, n);
			TMatrixD covariance(m, m);
			TMatrixD mcstatcovariance(m, m);			
			TMatrixD statcovariance(m, m);
			TMatrixD systcovariance(m, m);

			// Convert input into mathematical formats, easy and clean to be processed. 
			// Converted defined/implemented in source files, see include/Util.h

			H2V(PlotsTrue[4][WhichPlot], signal);
			H2V(QEPlotsTrue[4][WhichPlot], QEsignal);
			H2V(MECPlotsTrue[4][WhichPlot], MECsignal);
			H2V(RESPlotsTrue[4][WhichPlot], RESsignal);
			H2V(DISPlotsTrue[4][WhichPlot], DISsignal);
			H2V(COHPlotsTrue[4][WhichPlot], COHsignal);

			H2V(PlotsTrue[5][WhichPlot], altsignal);
			H2V(QEPlotsTrue[5][WhichPlot], QEaltsignal);
			H2V(MECPlotsTrue[5][WhichPlot], MECaltsignal);
			H2V(RESPlotsTrue[5][WhichPlot], RESaltsignal);
			H2V(DISPlotsTrue[5][WhichPlot], DISaltsignal);
			H2V(COHPlotsTrue[5][WhichPlot], COHaltsignal);			

			H2V(DataPlot, measure);
			H2M(ResponseMatrices[WhichPlot], response, kFALSE); // X axis: Reco, Y axis: True
			H2M(CovarianceMatrices[WhichPlot], covariance, kTRUE); // X axis: True, Y axis: Reco
			H2M(MCStatCovarianceMatrices[WhichPlot], mcstatcovariance, kTRUE); // X axis: True, Y axis: Reco

			// ------------------------------------------------------------------------------------------

			// Construct to record additinal smearing matrix and wiener filter (diagomal matrix) elements. 
    
			TMatrixD AddSmear(n,n);
			TVectorD WF(n);
			TMatrixD UnfoldCov(n,n);
			TMatrixD CovRotation(n,n);

			// ------------------------------------------------------------------------------------------	

			TH2D* smear = new TH2D("smear_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+XTitle,n,Nuedges,n,Nuedges);
			TH1D* wiener = new TH1D("wiener_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"Wiener Filter Vector",n,0,n);
			TH2D* unfcov = new TH2D("unfcov_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"Unfolded spectrum covariance", n, Nuedges, n, Nuedges);

			// --------------------------------------------------------------------------------------------------

			// Core implementation of Wiener-SVD
			// AddSmear and WF to record the core information in the unfolding.

			TVectorD unfold = WienerSVD(response, signal, measure, covariance, 2, 0.5, AddSmear, WF, UnfoldCov,CovRotation);

			// --------------------------------------------------------------------------------------------------

			// Start plotting
		
			TString CanvasName = PlotNames[WhichPlot]+"_"+Runs[WhichRun];
			TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.16);
			PlotCanvas->SetTopMargin(0.13);			
			PlotCanvas->SetLeftMargin(0.21);			
		
			TH1D* unf = new TH1D("unf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);

			// --------------------------------------------------------------------------------------------------

			V2H(unfold, unf); 	
			M2H(UnfoldCov, unfcov);		

			// --------------------------------------------------------------------------------------------------						

			// Setting the uncertainty using the fake data matrix: Stat + MC Stat + XSec

			for (int i = 1; i <= n;i++ ) { 

				// default / total uncertainty
				// XSec / Stat / MC Stat unc
				double CovUnc = TMath::Sqrt(UnfoldCov(i-1,i-1) );
				unf->SetBinError(i, CovUnc );							

			}

			ReweightXSec(unf);		

			// ------------------------------------------------------------------------------------

			// Start plotting 				

			unf->GetXaxis()->CenterTitle();
			unf->GetXaxis()->SetLabelSize(TextSize);
			unf->GetXaxis()->SetLabelFont(FontStyle);
			unf->GetXaxis()->SetTitleSize(TextSize);
			unf->GetXaxis()->SetTitleFont(FontStyle);			
			unf->GetXaxis()->SetNdivisions(6);	

			unf->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);
			unf->GetYaxis()->CenterTitle();
			unf->GetYaxis()->SetLabelSize(TextSize);
			unf->GetYaxis()->SetLabelFont(FontStyle);
			unf->GetYaxis()->SetTitleSize(TextSize);
			unf->GetYaxis()->SetTitleFont(FontStyle);			
			unf->GetYaxis()->SetNdivisions(6);

			unf->SetLineColor(BeamOnColor);
			unf->SetMarkerColor(BeamOnColor);
			unf->SetMarkerStyle(20);
			unf->SetMarkerSize(1.);	

			//----------------------------------------//

			// Unfolded MC Stat covariance
			TMatrixD CovRotation_T (TMatrixD::kTransposed, CovRotation); 
			TMatrixD UnfMCStatCov = CovRotation*mcstatcovariance*CovRotation_T; 				

			TH1D* unfMCStat = (TH1D*)(unf->Clone());	

			for (int i = 1; i <= n;i++ ) { 

				double BinWidth = unf->GetBinWidth(i);			

				// Stat, MC Stat, XSec uncertainties
				double CovUnc = TMath::Sqrt(UnfoldCov(i-1,i-1) ) / BinWidth;					
				unfMCStat->SetBinError(i, CovUnc );												

			}						

			// ----------------------------------------//		

			// The Nominal MC CC1p prediction has to be multiplied by the additional smearing matrix Ac

			TH1D* TrueUnf = new TH1D("TrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TVectorD AcTrueUnfold = AddSmear * signal;
			V2H(AcTrueUnfold, TrueUnf);

			TH1D* QETrueUnf = new TH1D("QETrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TVectorD QEAcTrueUnfold = AddSmear * QEsignal;
			V2H(QEAcTrueUnfold, QETrueUnf);

			TH1D* MECTrueUnf = new TH1D("MECTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TVectorD MECAcTrueUnfold = AddSmear * MECsignal;
			V2H(MECAcTrueUnfold, MECTrueUnf);

			TH1D* RESTrueUnf = new TH1D("RESTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TVectorD RESAcTrueUnfold = AddSmear * RESsignal;
			V2H(RESAcTrueUnfold, RESTrueUnf);

			TH1D* DISTrueUnf = new TH1D("DISTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TVectorD DISAcTrueUnfold = AddSmear * DISsignal;
			V2H(DISAcTrueUnfold, DISTrueUnf);

			TH1D* COHTrueUnf = new TH1D("COHTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TVectorD COHAcTrueUnfold = AddSmear * COHsignal;
			V2H(COHAcTrueUnfold, COHTrueUnf);			

			ReweightXSec(TrueUnf);
			ReweightXSec(QETrueUnf);
			ReweightXSec(MECTrueUnf);			
			ReweightXSec(RESTrueUnf);
			ReweightXSec(DISTrueUnf);			
			ReweightXSec(COHTrueUnf);
			
			TrueUnf->SetLineColor(OverlayColor);
			TrueUnf->SetMarkerColor(OverlayColor);	

			// ----------------------------------------//

			// Same for the alternative MC CC1p prediction (NuWro et al)

			TH1D* NoSmearAltTrueUnf = new TH1D("NoSmearAltTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* AltTrueUnf = new TH1D("AltTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TVectorD NoSmearAltTrueUnfold = altsignal;
			TVectorD AltAcTrueUnfold = AddSmear * altsignal;
			V2H(NoSmearAltTrueUnfold, NoSmearAltTrueUnf);
			V2H(AltAcTrueUnfold, AltTrueUnf);

			TH1D* QENoSmearAltTrueUnf = new TH1D("QENoSmearAltTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* QEAltTrueUnf = new TH1D("QEAltTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TVectorD QENoSmearAltTrueUnfold = QEaltsignal;
			TVectorD QEAltAcTrueUnfold = AddSmear * QEaltsignal;
			V2H(QENoSmearAltTrueUnfold, QENoSmearAltTrueUnf);
			V2H(QEAltAcTrueUnfold, QEAltTrueUnf);

			TH1D* MECNoSmearAltTrueUnf = new TH1D("MECNoSmearAltTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* MECAltTrueUnf = new TH1D("MECAltTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TVectorD MECNoSmearAltTrueUnfold = MECaltsignal;
			TVectorD MECAltAcTrueUnfold = AddSmear * MECaltsignal;
			V2H(MECNoSmearAltTrueUnfold, MECNoSmearAltTrueUnf);
			V2H(MECAltAcTrueUnfold, MECAltTrueUnf);

			TH1D* RESNoSmearAltTrueUnf = new TH1D("RESNoSmearAltTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* RESAltTrueUnf = new TH1D("RESAltTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TVectorD RESNoSmearAltTrueUnfold = RESaltsignal;
			TVectorD RESAltAcTrueUnfold = AddSmear * RESaltsignal;
			V2H(RESNoSmearAltTrueUnfold, RESNoSmearAltTrueUnf);
			V2H(RESAltAcTrueUnfold, RESAltTrueUnf);

			TH1D* DISNoSmearAltTrueUnf = new TH1D("DISNoSmearAltTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* DISAltTrueUnf = new TH1D("DISAltTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TVectorD DISNoSmearAltTrueUnfold = DISaltsignal;
			TVectorD DISAltAcTrueUnfold = AddSmear * DISaltsignal;
			V2H(DISNoSmearAltTrueUnfold, DISNoSmearAltTrueUnf);
			V2H(DISAltAcTrueUnfold, DISAltTrueUnf);

			TH1D* COHNoSmearAltTrueUnf = new TH1D("COHNoSmearAltTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* COHAltTrueUnf = new TH1D("COHAltTrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TVectorD COHNoSmearAltTrueUnfold = COHaltsignal;
			TVectorD COHAltAcTrueUnfold = AddSmear * COHaltsignal;
			V2H(COHNoSmearAltTrueUnfold, COHNoSmearAltTrueUnf);
			V2H(COHAltAcTrueUnfold, COHAltTrueUnf);															

			ReweightXSec(NoSmearAltTrueUnf);
			ReweightXSec(AltTrueUnf);
			ReweightXSec(QENoSmearAltTrueUnf);
			ReweightXSec(QEAltTrueUnf);
			ReweightXSec(MECNoSmearAltTrueUnf);
			ReweightXSec(MECAltTrueUnf);
			ReweightXSec(RESNoSmearAltTrueUnf);
			ReweightXSec(RESAltTrueUnf);
			ReweightXSec(DISNoSmearAltTrueUnf);
			ReweightXSec(DISAltTrueUnf);
			ReweightXSec(COHNoSmearAltTrueUnf);
			ReweightXSec(COHAltTrueUnf);															

			AltTrueUnf->SetLineColor(kOrange+7);
			AltTrueUnf->SetMarkerColor(kOrange+7);	

			// ----------------------------------------//	

			double max = 1.3 * TMath::Max(unf->GetMaximum(), AltTrueUnf->GetMaximum());		
//			unf->GetYaxis()->SetRangeUser(XSecRange[PlotNames[WhichPlot]].first,max);
			unfMCStat->GetYaxis()->SetRangeUser(XSecRange[PlotNames[WhichPlot]].first,max);	

			//------------------------------//

			// The N-dimensional analysis has been developed based on the bin number, not the actual range

			if (string(PlotNames[WhichPlot]).find("Serial") != std::string::npos) {	

				TString XaxisTitle = unfMCStat->GetXaxis()->GetTitle();
				XaxisTitle.ReplaceAll("deg","bin #");
				XaxisTitle.ReplaceAll("GeV/c","bin #");
				XaxisTitle.ReplaceAll("GeV","bin #");				
				unfMCStat->GetXaxis()->SetTitle(XaxisTitle);

				TString YaxisTitle = VarLabel[PlotNames[WhichPlot]];
				YaxisTitle.ReplaceAll("deg","");
				YaxisTitle.ReplaceAll("GeV/c","");
				YaxisTitle.ReplaceAll("GeV","");
				YaxisTitle.ReplaceAll("/c","");
				YaxisTitle.ReplaceAll("^{2}","");												
				unfMCStat->GetYaxis()->SetTitle(YaxisTitle);				

			}			

			//------------------------------//					

			// Draw the data points first to get the beautiful canvas 
			PlotCanvas->cd();
			//unf->Draw("e1x0"); // Full unc : XSec + Stat + MC Stat
			unfMCStat->Draw("e1x0");
			TrueUnf->Draw("hist same");
			AltTrueUnf->Draw("hist same");

			// Plotting again so that the data points are on top 
			PlotCanvas->cd();
			//unf->Draw("e1x0 same"); // Full unc : XSec + Stat + MC Stat
			unfMCStat->Draw("e1x0 same"); // Only MC Stat

			//------------------------------//
			/*
			if (BeamOnSample == "Overlay9NuWro" && OverlaySample == "Overlay9") {			

				// Only for the NuWro fake data study
				// use any residual (GENIE - NuWro)/ sqrt(12) diffs
				// as an unfolding model uncertainty

				TH1D* UnfUnc = (TH1D*)(TrueUnf->Clone());

				for (int ibin = 1; ibin <= (int)(TrueUnf->GetXaxis()->GetNbins()); ibin++) {

					double UnfPoint = unfMCStat-> GetBinContent(ibin);
					double UnfError = unfMCStat-> GetBinError(ibin);
					double TruePoint = AltTrueUnf-> GetBinContent(ibin);								

					if ( TMath::Abs(UnfPoint - TruePoint) <  UnfError) {

						UnfUnc->SetBinContent(ibin,0.);

					} else {

						double Uncertainty = TMath::Abs( TMath::Abs(UnfPoint - TruePoint) - UnfError ) / TMath::Sqrt(12);
						UnfUnc->SetBinContent(ibin,Uncertainty);

					}

				}

				UnfUnc->SetLineColor(kRed+1);
				UnfUnc->SetFillColor(kRed+1);			
				UnfUnc->Draw("e2 same");

				fUnc->cd();
				UnfUnc->Write("UnfUnc_"+PlotNames[WhichPlot]);

			}
			*/

			//------------------------------//

			// Legend & POT Normalization

			double tor860_wcut = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);
			TString Label = ToStringPOT(tor860_wcut)+" POT";

			TLegend* legData = new TLegend(0.23,0.89,0.95,0.98);
			legData->SetBorderSize(0);
			legData->SetTextSize(0.05);
			legData->SetTextFont(FontStyle);
			legData->SetNColumns(2);
			legData->SetMargin(0.1);

			//------------------------------//

			// Chi2 calculation

			double FDChi2, FDpval; int FDNdof;
			double CVChi2, CVpval; int CVNdof;	

			TH2D* CovClone = (TH2D*)(unfcov->Clone());

			for (int ix = 1; ix <= n; ix++) {

				for (int iy = 1; iy <= n; iy++) {

					double WidthX = unfcov->GetXaxis()->GetBinWidth(ix);
					double WidthY = unfcov->GetYaxis()->GetBinWidth(iy);

					double TwoDWidth = WidthX * WidthY;

					double BinContent = unfcov->GetBinContent(ix,iy);	
					double NewBinContent = BinContent/TwoDWidth;																	

					// Only for the diagonal elements
					// Add the unfolding uncertainty
					// On top of everything else
					// That is done both for the final xsec result and for the unfolded covariance
					if (ix == iy) { 
						// unfolded covariance matrix
//						double UnfUncBin = UnfUnc[WhichPlot]->GetBinContent(ix);
						double UnfUncBin = 0.;

						NewBinContent = NewBinContent + TMath::Power(UnfUncBin,2.) ;					 

						// xsec uncertainty
						double CurrentUnc = PlotsReco[0][WhichPlot]->GetBinError(ix);
						double NewError = TMath::Sqrt( TMath::Power(CurrentUnc,2.) + TMath::Power(UnfUncBin,2.) ) ;
						PlotsReco[0][WhichPlot]->SetBinError(ix,NewError);
						
					}

					CovClone->SetBinContent(ix,iy,NewBinContent);										

				}	

			}		

			CalcChiSquared(TrueUnf,unfMCStat,CovClone,CVChi2,CVNdof,CVpval);	
			TString CVChi2NdofAlt = "(" + to_string_with_precision(CVChi2,1) + "/" + TString(std::to_string(CVNdof)) +")";							
			CalcChiSquared(AltTrueUnf,unfMCStat,CovClone,FDChi2,FDNdof,FDpval);	
			TString FDChi2NdofAlt = "(" + to_string_with_precision(FDChi2,1) + "/" + TString(std::to_string(FDNdof)) +")";			

			//------------------------------//				

			TString ReducedLabel = 	BeamOnSample;
			TString PrintLabel = ReducedLabel.ReplaceAll("Overlay9","");	

			legData->AddEntry(unf,"Unfolded " + PrintLabel,"ep");

			TLegendEntry* lMC = legData->AddEntry(TrueUnf,"GENIE Tune "+ CVChi2NdofAlt,"l");
			lMC->SetTextColor(OverlayColor);

			TLegendEntry* lAltMC = legData->AddEntry(AltTrueUnf,"True " + PrintLabel + " " + FDChi2NdofAlt,"l");
			lAltMC->SetTextColor(kOrange+7);				

			legData->Draw();
			
			TString CanvasPath = PlotPath+"Overlay9";
			TString FullCanvasName = "/"+BeamOnSample+"WienerSVD_XSections_"+CanvasName+"_"+UBCodeVersion+Subtract+".pdf";
			if (OverlaySample != "Overlay9") { FullCanvasName = "/"+BeamOnSample+"WienerSVD_XSections_"+CanvasName+"_"+UBCodeVersion+Subtract+"_"+OverlaySample+".pdf"; }
			PlotCanvas->SaveAs(CanvasPath+FullCanvasName);	
			delete PlotCanvas;			

			//------------------------------//

			TH1D* diff = new TH1D("diff_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"Fractional difference of unf and signal model",n, Nuedges);
			
			for(int i=1; i <= n; i++) {
			
				double s2s = 0;
				double u = unf->GetBinContent(i);
				double t = signal(i-1);
				if (t != 0) { s2s = u/t - 1; }
				else { s2s = 1.; } // attention to this t = 0
				diff->SetBinContent(i, s2s); // in percentage 

			}	

			//------------------------------//
/*
			// intrinsic bias (Ac-I) * s_bar formula

			TH1D* bias = new TH1D("bias_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"intrinsic bias w.r.t. model",n, Nuedges);
			TH1D* bias2 = new TH1D("bias2_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"intrinsic bias2 w.r.t. unfolded result",n, Nuedges);
			TMatrixD unit(n,n);
			unit.UnitMatrix();
			TVectorD intrinsicbias = (AddSmear - unit)*signal;
			TVectorD intrinsicbias2 = (AddSmear - unit)*unfold;

			for (int i = 0; i < n; i++) {

			    if (signal(i) != 0) { intrinsicbias(i) = intrinsicbias(i)/signal(i); }
			    else { intrinsicbias(i) = 0.; }
			    if (unfold(i) != 0) { intrinsicbias2(i) = intrinsicbias2(i)/unfold(i); }
			    else { intrinsicbias2(i) = 0.; }
			}	

			V2H(intrinsicbias, bias);
			V2H(intrinsicbias2, bias2);				

			//------------------------------//

			// Diagonal Uncertainty

			TH1D* fracError = new TH1D("fracError_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun], "Fractional uncertainty", n, Nuedges);
			TH1D* absError = new TH1D("absError_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun], "absolute uncertainty", n, Nuedges);

			for (int i = 1; i <= n; i++) {

			    fracError->SetBinContent(i, TMath::Sqrt(UnfoldCov(i-1, i-1))/unfold(i-1));
			    absError->SetBinContent(i, TMath::Sqrt(UnfoldCov(i-1, i-1)));
			}
    
			/// MSE
			TH1D* MSE = new TH1D("MSE_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun], "Mean Square Error: variance+bias^2", n, Nuedges);
			TH1D* MSE2 = new TH1D("MSE2_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun], "Mean Square Error: variance", n, Nuedges);
    
			for (int i = 0; i < n; i++) {

			    MSE->SetBinContent(i+1, TMath::Power(intrinsicbias2(i)*unfold(i),2)+UnfoldCov(i,i));
			    MSE2->SetBinContent(i+1, UnfoldCov(i,i));

			}

			// convert matrix/vector to histogram and save

			M2H(AddSmear, smear);
			V2H(WF, wiener);
*/
			//------------------------------//

			ExtractedXSec->cd();

			unf->Write("Reco"+PlotNames[WhichPlot]); // "Data" xsec
			unfMCStat->Write("MCStatReco"+PlotNames[WhichPlot]); // "Data" only with MC stat uncertainties		
			TrueUnf->Write("True"+PlotNames[WhichPlot]); // MC after multiplication by Ac
			QETrueUnf->Write("QETrue"+PlotNames[WhichPlot]);
			MECTrueUnf->Write("MECTrue"+PlotNames[WhichPlot]);
			RESTrueUnf->Write("RESTrue"+PlotNames[WhichPlot]);
			DISTrueUnf->Write("DISTrue"+PlotNames[WhichPlot]);				
			COHTrueUnf->Write("COHTrue"+PlotNames[WhichPlot]);			
			NoSmearAltTrueUnf->Write("NoSmearAltTrue"+PlotNames[WhichPlot]); // alternative MC w/o multiplication by Ac
			QENoSmearAltTrueUnf->Write("QENoSmearAltTrue"+PlotNames[WhichPlot]);	
			MECNoSmearAltTrueUnf->Write("MECNoSmearAltTrue"+PlotNames[WhichPlot]);
			RESNoSmearAltTrueUnf->Write("RESNoSmearAltTrue"+PlotNames[WhichPlot]);
			DISNoSmearAltTrueUnf->Write("DISNoSmearAltTrue"+PlotNames[WhichPlot]);
			COHNoSmearAltTrueUnf->Write("COHNoSmearAltTrue"+PlotNames[WhichPlot]);									
			AltTrueUnf->Write("AltTrue"+PlotNames[WhichPlot]); // alternative MC with multiplication by Ac
			QEAltTrueUnf->Write("QEAltTrue"+PlotNames[WhichPlot]);
			MECAltTrueUnf->Write("MECAltTrue"+PlotNames[WhichPlot]);
			RESAltTrueUnf->Write("RESAltTrue"+PlotNames[WhichPlot]);
			DISAltTrueUnf->Write("DISAltTrue"+PlotNames[WhichPlot]);
			COHAltTrueUnf->Write("COHAltTrue"+PlotNames[WhichPlot]);												
			smear->Write("Ac"+PlotNames[WhichPlot]); // additional smearing matrix
			unfcov->Write("UnfCov"+PlotNames[WhichPlot]); // unfolded covariance matrix
			CovarianceMatrices[WhichPlot]->Write("Cov"+PlotNames[WhichPlot]); // covariance matrix		
			//wiener->Write("Wiener"+PlotNames[WhichPlot]);
			//diff->Write("Diff"+PlotNames[WhichPlot]);
			//bias->Write("Bias"+PlotNames[WhichPlot]);
			//bias2->Write("Bias2"+PlotNames[WhichPlot]);
			//fracError->Write("FracErr"+PlotNames[WhichPlot]);
			//absError->Write("AbsErr"+PlotNames[WhichPlot]);
			//MSE->Write("MSE"+PlotNames[WhichPlot]);
			//MSE2->Write("MSE2"+PlotNames[WhichPlot]);
			ResponseMatrices[WhichPlot]->Write("Response"+PlotNames[WhichPlot]); // response matrix

			// ---------------------------------------------------------------------------------------------------------------------------
/*
			// Make the additional smearing matrix pretty

			TString SmearCanvasName = "Smear_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun];
			TCanvas* SmearPlotCanvas = new TCanvas(SmearCanvasName,SmearCanvasName,205,34,1024,768);
			SmearPlotCanvas->cd();
			SmearPlotCanvas->SetBottomMargin(0.16);
			SmearPlotCanvas->SetTopMargin(0.13);			
			SmearPlotCanvas->SetLeftMargin(0.2);
			SmearPlotCanvas->SetRightMargin(0.2);

			smear->GetXaxis()->CenterTitle();
			smear->GetXaxis()->SetLabelFont(FontStyle);
			smear->GetXaxis()->SetTitleFont(FontStyle);
			smear->GetXaxis()->SetLabelSize(TextSize);
			smear->GetXaxis()->SetTitleSize(TextSize);
			smear->GetXaxis()->SetNdivisions(6);

			smear->GetYaxis()->CenterTitle();
			smear->GetYaxis()->SetLabelFont(FontStyle);
			smear->GetYaxis()->SetTitleFont(FontStyle);
			smear->GetYaxis()->SetLabelSize(TextSize);
			smear->GetYaxis()->SetTitleSize(TextSize);
			smear->GetYaxis()->SetNdivisions(6);	

			smear->GetZaxis()->SetLabelFont(FontStyle);
			smear->GetZaxis()->SetTitleFont(FontStyle);
			smear->GetZaxis()->SetLabelSize(TextSize);
			smear->GetZaxis()->SetTitleSize(TextSize);

			smear->Draw("coltz");

			TLatex *text = new TLatex();
			text->SetTextFont(FontStyle);
			text->SetTextSize(0.06);
			text->DrawTextNDC(0.47, 0.92, Runs[WhichRun]);

			TString SmearCanvas = "/"+BeamOnSample+"Smear_WienerSVD_XSections_"+CanvasName+"_"+UBCodeVersion+Subtract+".pdf";
			if (OverlaySample != "Overlay9") { SmearCanvas = "/"+BeamOnSample+"Smear_WienerSVD_XSections_"+CanvasName+"_"+UBCodeVersion+Subtract+"_"+OverlaySample+".pdf"; }
			SmearPlotCanvas->SaveAs(CanvasPath+SmearCanvas);
			delete SmearPlotCanvas;
*/
			// ---------------------------------------------------------------------------------------------------------------------------

		} // End of the loop over the plots

		//if (BeamOnSample == "Overlay9NuWro" && OverlaySample == "Overlay9") {

		//	cout << endl << "File " << fUncName << " created" << endl;
		//	fUnc->Close();

		//}

		ExtractedXSec->Close();
		std::cout << std::endl << "File " << NameExtractedXSec << " created" << std::endl << std::endl;

	} // End of the loop over the runs	

} // End of the program 
