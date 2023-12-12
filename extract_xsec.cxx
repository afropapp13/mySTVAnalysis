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

#include "../myClasses/Tools.h"
#include "../myClasses/Constants.h"
#include "../myClasses/Util.h"
#include "../myClasses/WienerSVD.h"

using namespace std;
using namespace Constants;

#include "../myClasses/myFunctions.cpp"

//----------------------------------------//

// For multi dimentional development
// Don't forget to divide by the slice range
// That is NOT done in this module !!!

//----------------------------------------//

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

void extract_xsec(TString OverlaySample = "", bool ClosureTest = false, TString BeamOnSample = "", TString Tune = "") {

	// -------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();	
	gStyle->SetOptStat(0);
	gStyle->SetEndErrorSize(6);	

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
	VectorCuts.push_back("_CRT");

	int NCuts = (int)(VectorCuts.size());	

	for (int i = 0; i < NCuts; i++) { CutExtension = CutExtension + VectorCuts[i]; }

	// -------------------------------------------------------------------------------------

//	vector<TString> PlotNames;
//	PlotNames.push_back("DeltaPTPlot"); 

	const int N1DPlots = PlotNames.size();
	//cout << "Number of 1D Plots = " << N1DPlots << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");

	vector<TString> LabelsOfSamples;
	vector<TString> NameOfSamples;

	NameOfSamples.push_back("Overlay9"); if (OverlaySample != "") { NameOfSamples[0] = OverlaySample; }
	NameOfSamples.push_back("BeamOn9");  if (BeamOnSample != "") { NameOfSamples[1] = BeamOnSample; }
	NameOfSamples.push_back("ExtBNB9"); 
	NameOfSamples.push_back("OverlayDirt9");
	NameOfSamples.push_back("GenieOverlay");	 
	
	int DataIndex = -1.;

	// -----------------------------------------------------------------------------------------------------------------------------------

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

	// For the systematics decomposition into stat /shape / norm 

	Tools tools;

	// -------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------
	
		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);						
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface);
		//cout << Runs[WhichRun] << " Integrated flux = " << IntegratedFlux << endl;
				
		// -------------------------------------------------------------------------------------		

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

		vector<TH2D*> ResponseMatrices; ResponseMatrices.clear();

		vector<TH2D*> CovarianceMatrices; CovarianceMatrices.clear();
		vector<TH2D*> StatCovarianceMatrices; StatCovarianceMatrices.clear();
		vector<TH2D*> MCStatCovarianceMatrices; MCStatCovarianceMatrices.clear();			
		vector<TH2D*> SystCovarianceMatrices; SystCovarianceMatrices.clear();	
		vector<TH2D*> LYCovarianceMatrices; LYCovarianceMatrices.clear();
		vector<TH2D*> TPCCovarianceMatrices; TPCCovarianceMatrices.clear();
		vector<TH2D*> SCERecomb2CovarianceMatrices; SCERecomb2CovarianceMatrices.clear();
		vector<TH2D*> XSecCovarianceMatrices; XSecCovarianceMatrices.clear();
		vector<TH2D*> FluxCovarianceMatrices; FluxCovarianceMatrices.clear();
		vector<TH2D*> G4CovarianceMatrices; G4CovarianceMatrices.clear();
		vector<TH2D*> DirtCovarianceMatrices; DirtCovarianceMatrices.clear();
		vector<TH2D*> POTCovarianceMatrices; POTCovarianceMatrices.clear();
		vector<TH2D*> NTargetCovarianceMatrices; NTargetCovarianceMatrices.clear();	
		//vector<TH2D*> NuWroCovarianceMatrices; NuWroCovarianceMatrices.clear();

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString FileResponseName = MigrationMatrixPath+Tune+"FileResponseMatrices_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
		TFile* FileResponseMatrices = new TFile(FileResponseName,"readonly");

		TString FileCovarianceName = MigrationMatrixPath+Tune+"WienerSVD_Total_CovarianceMatrices_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
		TFile* FileCovarianceMatrices = new TFile(FileCovarianceName,"readonly");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();
		
		TString PathToFilesUBCodeExtension = PathToFiles+CutExtension;

		// Store the extracted xsections & associated files in dedicated file

		TFile* ExtractedXSec = nullptr;
		TString NameExtractedXSec = "";

		// File to store xsecs for data release

		TString XSecTxtName = PathToExtractedXSec+"TxtXSec_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+Subtract+".txt";
		ofstream myXSecTxtFile;

		if (ClosureTest == false) {

			NameExtractedXSec = PathToExtractedXSec+Tune+"WienerSVD_ExtractedXSec_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+Subtract+".root";
			if (Tune != "") 
				{ NameExtractedXSec = PathToExtractedXSec+"AltMC"+Tune+"WienerSVD_ExtractedXSec_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+Subtract+".root"; }
			ExtractedXSec = TFile::Open(NameExtractedXSec,"recreate");

			// ---------------------------------------------------------------------------------------------------------------------------------------

			//File to store xsecs for data release

			myXSecTxtFile.open(XSecTxtName);
			myXSecTxtFile << std::fixed << std::setprecision(2);
			myXSecTxtFile << Runs[WhichRun] << endl << endl;
			myXSecTxtFile << "Bin #; Bin Low Edge; Bin High Edge; Bin Entry [10^{-38} cm^{2}]; Bin Error [10^{-38} cm^{2}]" << endl << endl;

		}		

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Loop over the samples

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {
		
			if (NameOfSamples[WhichSample] == "BeamOn9") { DataIndex = WhichSample; }

			if (
				NameOfSamples[WhichSample] == "BeamOn9" || 
				NameOfSamples[WhichSample] == "ExtBNB9" || 
				NameOfSamples[WhichSample] == "OverlayDirt9"
			) { 
			
				TString FileName = "STVStudies_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+CutExtension+".root";
				FileSample.push_back(TFile::Open(PathToFilesUBCodeExtension+"/"+FileName)); 
			}
			
			if (NameOfSamples[WhichSample] == "Overlay9") { 
			
				TString FileName = Tune+"STVStudies_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+CutExtension+".root";
				FileSample.push_back(TFile::Open(PathToFilesUBCodeExtension+"/"+FileName)); 
				
			}

			if (NameOfSamples[WhichSample] == "GenieOverlay") { 
			
				TString FileName = Tune+"TruthSTVAnalysis_Overlay9_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
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

		// ----------------------------------------------------------------------------------------------------------------------------------

		// Loop over the event rate plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

			// ------------------------------------------------------------------------------------------------

			ResponseMatrices.push_back((TH2D*)FileResponseMatrices->Get("POTScaledCC1pReco"+PlotNames[WhichPlot]+"2D"));

			// Already flux-averaged rates
			CovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("TotalCovariance_"+PlotNames[WhichPlot]));
			SystCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("SystCovariance_"+PlotNames[WhichPlot]));
			StatCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("StatCovariance_"+PlotNames[WhichPlot]));
			MCStatCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("MCStatCovariance_"+PlotNames[WhichPlot]));

			LYCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("LYCovariance_"+PlotNames[WhichPlot]));
			TPCCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("TPCCovariance_"+PlotNames[WhichPlot]));
			SCERecomb2CovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("SCERecomb2Covariance_"+PlotNames[WhichPlot]));
			XSecCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("XSecCovariance_"+PlotNames[WhichPlot]));
			FluxCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("FluxCovariance_"+PlotNames[WhichPlot]));
			G4CovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("G4Covariance_"+PlotNames[WhichPlot]));
			DirtCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("DirtCovariance_"+PlotNames[WhichPlot]));
			POTCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("POTCovariance_"+PlotNames[WhichPlot]));
			//NuWroCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("NuWroCovariance_"+PlotNames[WhichPlot]));			
			NTargetCovarianceMatrices.push_back((TH2D*)FileCovarianceMatrices->Get("NTargetCovariance_"+PlotNames[WhichPlot]));

			// -----------------------------------------------------------------------------------------------------

			// True CC1p Signal MC // No detector/ reconstruction / smearing effects

			int n = PlotsTrue[4][WhichPlot]->GetNbinsX();
			double Nuedges[n+1];
			    
			for (int i = 0; i < n+1; i++) { Nuedges[i] = PlotsTrue[4][WhichPlot]->GetBinLowEdge(i+1); }

			// -------------------------------------------------------------------------------------------------

			// BeamOn

			TH1D* DataPlot = (TH1D*)(PlotsReco[1][WhichPlot]->Clone());

			DataPlot->Add(PlotsReco[2][WhichPlot],-1); // Subtract ExtBNB
			DataPlot->Add(PlotsReco[3][WhichPlot],-1); // Subtract Dirt
			DataPlot->Add(PlotsBkgReco[0][WhichPlot],-1); // Subtract NonCC1p Beam Related Background

			// If performing a closure test, use the CC1p0pi MC part and no bkg subtraction

			if (ClosureTest == true) { DataPlot = PlotsCC1pReco[0][WhichPlot]; }

			int m = DataPlot->GetNbinsX();			
			TString XTitle = DataPlot->GetXaxis()->GetTitle();
			TString YTitle = DataPlot->GetYaxis()->GetTitle();	

			// Flux-averaged event rates 
			// both for the reco and for the true level spectrum
			DataPlot->Scale(Units/(IntegratedFlux*NTargets));		
			PlotsTrue[4][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
			QEPlotsTrue[4][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
			MECPlotsTrue[4][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
			RESPlotsTrue[4][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
			DISPlotsTrue[4][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
			COHPlotsTrue[4][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));																	 

			// -------------------------------------------------------------------------------------------

			// Construct vectors (for 1D histogram) and matrices (for 2D histogram) for input

			TVectorD signal(n);
			TVectorD QEsignal(n);
			TVectorD MECsignal(n);
			TVectorD RESsignal(n);
			TVectorD DISsignal(n);
			TVectorD COHsignal(n);															
			TVectorD measure(m);
			TMatrixD response(m, n);

			TMatrixD covariance(m, m);
			TMatrixD mcstatcovariance(m, m);
			TMatrixD statcovariance(m, m);
			TMatrixD systcovariance(m, m);
			TMatrixD lycovariance(m, m);	
			TMatrixD tpccovariance(m, m);
			TMatrixD scerecomb2covariance(m, m);
			TMatrixD fluxcovariance(m, m);	
			TMatrixD xseccovariance(m, m);
			TMatrixD g4covariance(m, m);
			TMatrixD dirtcovariance(m, m);	
			TMatrixD potcovariance(m, m);
			TMatrixD nuwrocovariance(m, m);			
			TMatrixD ntargetcovariance(m, m);

			TMatrixD signalcovariance(m, m);
			TMatrixD signalmcstatcovariance(m, m);
			TMatrixD signalstatcovariance(m, m);
			TMatrixD signalsystcovariance(m, m);
			TMatrixD signallycovariance(m, m);	
			TMatrixD signaltpccovariance(m, m);
			TMatrixD signalscerecomb2covariance(m, m);
			TMatrixD signalfluxcovariance(m, m);	
			TMatrixD signalxseccovariance(m, m);
			TMatrixD signalg4covariance(m, m);
			TMatrixD signaldirtcovariance(m, m);	
			TMatrixD signalpotcovariance(m, m);
			TMatrixD signalnuwrocovariance(m, m);			
			TMatrixD signalntargetcovariance(m, m);

			TMatrixD bkgcovariance(m, m);
			TMatrixD bkgmcstatcovariance(m, m);
			TMatrixD bkgstatcovariance(m, m);
			TMatrixD bkgsystcovariance(m, m);
			TMatrixD bkglycovariance(m, m);	
			TMatrixD bkgtpccovariance(m, m);
			TMatrixD bkgscerecomb2covariance(m, m);
			TMatrixD bkgfluxcovariance(m, m);	
			TMatrixD bkgxseccovariance(m, m);
			TMatrixD bkgg4covariance(m, m);
			TMatrixD bkgdirtcovariance(m, m);	
			TMatrixD bkgpotcovariance(m, m);
			TMatrixD bkgnuwrocovariance(m, m);			
			TMatrixD bkgntargetcovariance(m, m);															

			// Convert input into mathematical formats, easy and clean to be processed. 
			// Converted defined/implemented in source files, see include/Util.h

			H2V(PlotsTrue[4][WhichPlot], signal);
			H2V(QEPlotsTrue[4][WhichPlot], QEsignal);
			H2V(MECPlotsTrue[4][WhichPlot], MECsignal);
			H2V(RESPlotsTrue[4][WhichPlot], RESsignal);
			H2V(DISPlotsTrue[4][WhichPlot], DISsignal);
			H2V(COHPlotsTrue[4][WhichPlot], COHsignal);												

			H2V(DataPlot, measure);
			H2M(ResponseMatrices[WhichPlot], response, kFALSE); // X axis: Reco, Y axis: True

			H2M(CovarianceMatrices[WhichPlot], covariance, kTRUE); // X axis: True, Y axis: Reco
			H2M(MCStatCovarianceMatrices[WhichPlot], mcstatcovariance, kTRUE); // X axis: True, Y axis: Reco
			H2M(StatCovarianceMatrices[WhichPlot], statcovariance, kTRUE); // X axis: True, Y axis: Reco
			H2M(SystCovarianceMatrices[WhichPlot], systcovariance, kTRUE); // X axis: True, Y axis: Reco

			H2M(LYCovarianceMatrices[WhichPlot], lycovariance, kTRUE); // X axis: True, Y axis: Reco
			H2M(TPCCovarianceMatrices[WhichPlot], tpccovariance, kTRUE); // X axis: True, Y axis: Reco
			H2M(SCERecomb2CovarianceMatrices[WhichPlot], scerecomb2covariance, kTRUE); // X axis: True, Y axis: Reco
			H2M(FluxCovarianceMatrices[WhichPlot], fluxcovariance, kTRUE); // X axis: True, Y axis: Reco
			H2M(XSecCovarianceMatrices[WhichPlot], xseccovariance, kTRUE); // X axis: True, Y axis: Reco
			H2M(G4CovarianceMatrices[WhichPlot], g4covariance, kTRUE); // X axis: True, Y axis: Reco
			H2M(DirtCovarianceMatrices[WhichPlot], dirtcovariance, kTRUE); // X axis: True, Y axis: Reco
			H2M(POTCovarianceMatrices[WhichPlot], potcovariance, kTRUE); // X axis: True, Y axis: Reco
			//H2M(NuWroCovarianceMatrices[WhichPlot], nuwrocovariance, kTRUE); // X axis: True, Y axis: Reco			
			H2M(NTargetCovarianceMatrices[WhichPlot], ntargetcovariance, kTRUE); // X axis: True, Y axis: Reco

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
			TH2D* normunfcov = new TH2D("normunfcov_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"Norm Unfolded spectrum covariance", n, Nuedges, n, Nuedges);
			TH2D* shapeunfcov = new TH2D("shapeunfcov_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"Shape Unfolded spectrum covariance", n, Nuedges, n, Nuedges);						
			TH2D* covrot = new TH2D("covrot_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"Covariance rotation matrix", n, Nuedges, n, Nuedges);
			TH2D* transpcovrot = new TH2D("transpcovrot_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"Transpose covariance rotation matrix", n, Nuedges, n, Nuedges);			

			// --------------------------------------------------------------------------------------------------

			// Core implementation of Wiener-SVD
			// AddSmear and WF to record the core information in the unfolding.

			TVectorD unfold = WienerSVD(response, signal, measure, covariance, 2, 0.5, AddSmear, WF, UnfoldCov, CovRotation);

			// --------------------------------------------------------------------------------------------------

			TMatrixD CovRotation_T (TMatrixD::kTransposed, CovRotation); 

			TMatrixD UnfMCStatCov = CovRotation*mcstatcovariance*CovRotation_T;
			TMatrixD UnfStatCov = CovRotation*statcovariance*CovRotation_T; 
			TMatrixD UnfSystCov = CovRotation*systcovariance*CovRotation_T; 
			TMatrixD UnfLYCov = CovRotation*lycovariance*CovRotation_T;
			TMatrixD UnfTPCCov = CovRotation*tpccovariance*CovRotation_T;
			TMatrixD UnfSCERecomb2Cov = CovRotation*scerecomb2covariance*CovRotation_T;
			TMatrixD UnfFluxCov = CovRotation*fluxcovariance*CovRotation_T;
			TMatrixD UnfXSecCov = CovRotation*xseccovariance*CovRotation_T;
			TMatrixD UnfG4Cov = CovRotation*g4covariance*CovRotation_T;
			TMatrixD UnfDirtCov = CovRotation*dirtcovariance*CovRotation_T;
			TMatrixD UnfPOTCov = CovRotation*potcovariance*CovRotation_T;
			//TMatrixD UnfNuWroCov = CovRotation*nuwrocovariance*CovRotation_T;			
			TMatrixD UnfNTargetCov = CovRotation*ntargetcovariance*CovRotation_T;	

			// Decomposition of systematic uncertainties into shape / normalization uncertainty

			std::vector<TMatrixD> NormShapeVector = tools.MatrixDecomp(n,measure,UnfSystCov);

			// --------------------------------------------------------------------------------------------------

			// Start plotting
		
			TString CanvasName = PlotNames[WhichPlot]+"_"+Runs[WhichRun];
			TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.16);
			PlotCanvas->SetTopMargin(0.13);			
			PlotCanvas->SetLeftMargin(0.21);			
		
			TH1D* unf = new TH1D("unf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* unfFullUnc = new TH1D("unfFullUnc_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* unfShapeOnly = new TH1D("unfShapeOnly_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* unfNormOnly = new TH1D("unfNormOnly_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);			

			TH1D* unfMCStat = new TH1D("unfMCStat_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* unfStat = new TH1D("unfStat_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* unfLY = new TH1D("unfLY_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* unfTPC = new TH1D("unfTPC_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* unfSCERecomb2 = new TH1D("unfSCERecomb2_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* unfFlux = new TH1D("unfFlux_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* unfXSec = new TH1D("unfXSec_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* unfG4 = new TH1D("unfG4_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* unfDirt = new TH1D("unfDirt_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* unfPOT = new TH1D("unfPOT_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			//TH1D* unfNuWro = new TH1D("unfNuWro_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);			
			TH1D* unfNTarget = new TH1D("unfNTarget_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);	

			TH1D* signalunfMCStat = new TH1D("signalunfMCStat_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* signalunfStat = new TH1D("signalunfStat_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* signalunfLY = new TH1D("signalunfLY_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* signalunfTPC = new TH1D("signalunfTPC_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* signalunfSCERecomb2 = new TH1D("signalunfSCERecomb2_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* signalunfFlux = new TH1D("signalunfFlux_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* signalunfXSec = new TH1D("signalunfXSec_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* signalunfG4 = new TH1D("signalunfG4_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* signalunfDirt = new TH1D("signalunfDirt_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* signalunfPOT = new TH1D("signalunfPOT_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			//TH1D* signalunfNuWro = new TH1D("signalunfNuWro_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);			
			TH1D* signalunfNTarget = new TH1D("signalunfNTarget_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);

			TH1D* bkgunfMCStat = new TH1D("bkgunfMCStat_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* bkgunfStat = new TH1D("bkgunfStat_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* bkgunfLY = new TH1D("bkgunfLY_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* bkgunfTPC = new TH1D("bkgunfTPC_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* bkgunfSCERecomb2 = new TH1D("bkgunfSCERecomb2_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* bkgunfFlux = new TH1D("bkgunfFlux_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* bkgunfXSec = new TH1D("bkgunfXSec_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* bkgunfG4 = new TH1D("bkgunfG4_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* bkgunfDirt = new TH1D("bkgunfDirt_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TH1D* bkgunfPOT = new TH1D("bkgunfPOT_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			//TH1D* bkgunfNuWro = new TH1D("bkgunfNuWro_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);			
			TH1D* bkgunfNTarget = new TH1D("bkgunfNTarget_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);																	

			// --------------------------------------------------------------------------------------------------

			V2H(unfold, unf);

			// --------------------------------------------------------------------------------------------------					

			// Scaling by correct units & bin width division

			for (int i = 1; i <= n;i++ ) { 

				double CVInBin = unf->GetBinContent(i);

				// default / total uncertainty
				double CovUnc = TMath::Sqrt(UnfoldCov(i-1,i-1) );

				unf->SetBinError(i, CovUnc );

			}

			ReweightXSec(unf);
			ReweightXSec(PlotsTrue[4][WhichPlot]);
			ReweightXSec(QEPlotsTrue[4][WhichPlot]);
			ReweightXSec(MECPlotsTrue[4][WhichPlot]);
			ReweightXSec(RESPlotsTrue[4][WhichPlot]);
			ReweightXSec(DISPlotsTrue[4][WhichPlot]);
			ReweightXSec(COHPlotsTrue[4][WhichPlot]);																		

			unfFullUnc = (TH1D*)(unf->Clone());			
			unfShapeOnly = (TH1D*)(unf->Clone());
			unfNormOnly = (TH1D*)(unf->Clone());

			unfStat = (TH1D*)(unf->Clone());
			unfMCStat = (TH1D*)(unf->Clone());
			unfLY = (TH1D*)(unf->Clone());
			unfTPC = (TH1D*)(unf->Clone());
			unfSCERecomb2 = (TH1D*)(unf->Clone());
			unfFlux = (TH1D*)(unf->Clone());
			unfXSec = (TH1D*)(unf->Clone());
			unfG4 = (TH1D*)(unf->Clone());
			unfDirt = (TH1D*)(unf->Clone());
			unfPOT = (TH1D*)(unf->Clone());
			//unfNuWro = (TH1D*)(unf->Clone());			
			unfNTarget = (TH1D*)(unf->Clone());

			signalunfStat = (TH1D*)(unf->Clone());
			signalunfMCStat = (TH1D*)(unf->Clone());
			signalunfLY = (TH1D*)(unf->Clone());
			signalunfTPC = (TH1D*)(unf->Clone());
			signalunfSCERecomb2 = (TH1D*)(unf->Clone());
			signalunfFlux = (TH1D*)(unf->Clone());
			signalunfXSec = (TH1D*)(unf->Clone());
			signalunfG4 = (TH1D*)(unf->Clone());
			signalunfDirt = (TH1D*)(unf->Clone());
			signalunfPOT = (TH1D*)(unf->Clone());
			//signalunfNuWro = (TH1D*)(unf->Clone());			
			signalunfNTarget = (TH1D*)(unf->Clone());

			bkgunfStat = (TH1D*)(unf->Clone());
			bkgunfMCStat = (TH1D*)(unf->Clone());
			bkgunfLY = (TH1D*)(unf->Clone());
			bkgunfTPC = (TH1D*)(unf->Clone());
			bkgunfSCERecomb2 = (TH1D*)(unf->Clone());
			bkgunfFlux = (TH1D*)(unf->Clone());
			bkgunfXSec = (TH1D*)(unf->Clone());
			bkgunfG4 = (TH1D*)(unf->Clone());
			bkgunfDirt = (TH1D*)(unf->Clone());
			bkgunfPOT = (TH1D*)(unf->Clone());
			//bkgunfNuWro = (TH1D*)(unf->Clone());			
			bkgunfNTarget = (TH1D*)(unf->Clone());																		

			myXSecTxtFile << PlotNames[WhichPlot] << endl << endl;			

			for (int i = 1; i <= n;i++ ) {

				double CV = unf->GetBinContent(i);
				double CVError = unf->GetBinError(i);
				double LowEdge = unf->GetBinLowEdge(i);
				double Width = unf->GetBinWidth(i);
				double HighEdge = LowEdge + Width;

				double StatError = TMath::Sqrt( UnfStatCov(i-1,i-1) ) / Width;
				unfStat->SetBinError(i, StatError);

				double MCStatError = TMath::Sqrt( UnfMCStatCov(i-1,i-1) ) / Width;
				unfMCStat->SetBinError(i, MCStatError);				

				// Set unc = stat + shape syst
				double ShapeStatUnc = NormShapeVector[1](i-1,i-1) + UnfStatCov(i-1,i-1);
				unf->SetBinError(i, TMath::Sqrt( ShapeStatUnc ) / Width );	

				// Keep track of the total unc as well
				unfFullUnc->SetBinError(i, TMath::Sqrt(UnfoldCov(i-1,i-1)) / Width );				

				unfShapeOnly->SetBinError(i,TMath::Sqrt( TMath::Abs( NormShapeVector[1](i-1,i-1) ) ) / Width );	
				unfNormOnly->SetBinContent(i,0.);	
				unfNormOnly->SetBinError(i,TMath::Sqrt( TMath::Abs( NormShapeVector[0](i-1,i-1) ) ) / Width );

				//----------------------------------------//

				// Sources of uncertainty

				double LYError = TMath::Sqrt( UnfLYCov(i-1,i-1) ) / Width;
				unfLY->SetBinError(i, LYError);

				double TPCError = TMath::Sqrt( UnfTPCCov(i-1,i-1) ) / Width;
				unfTPC->SetBinError(i, TPCError);

				double SCERecomb2Error = TMath::Sqrt( UnfSCERecomb2Cov(i-1,i-1) ) / Width;
				unfSCERecomb2->SetBinError(i, SCERecomb2Error);

				double FluxError = TMath::Sqrt( UnfFluxCov(i-1,i-1) ) / Width;
				unfFlux->SetBinError(i, FluxError);

				double XSecError = TMath::Sqrt( UnfXSecCov(i-1,i-1) ) / Width;
				unfXSec->SetBinError(i, XSecError);

				double G4Error = TMath::Sqrt( UnfG4Cov(i-1,i-1) ) / Width;
				unfG4->SetBinError(i, G4Error);

				double DirtErrorCovEntry = UnfDirtCov(i-1,i-1);
				if (DirtErrorCovEntry < 0.) { DirtErrorCovEntry = 0.; }
				double DirtError = TMath::Sqrt( DirtErrorCovEntry ) / Width;
				unfDirt->SetBinError(i, DirtError);

				double POTError = TMath::Sqrt( UnfPOTCov(i-1,i-1) ) / Width;
				unfPOT->SetBinError(i, POTError);

				//double NuWroError = TMath::Sqrt( UnfNuWroCov(i-1,i-1) ) / Width;
				//unfNuWro->SetBinError(i, NuWroError);				

				double NTargetError = TMath::Sqrt( UnfNTargetCov(i-1,i-1) ) / Width;
				unfNTarget->SetBinError(i, NTargetError);																				

				//----------------------------------------//

				// Data release

				myXSecTxtFile << i << "; " << LowEdge << "; " << HighEdge << "; " << CV << "; " << CVError << endl << endl;			

			}															

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

			unf->GetYaxis()->SetRangeUser(XSecRange[PlotNames[WhichPlot]].first,XSecRange[PlotNames[WhichPlot]].second);
			unf->SetLineColor(BeamOnColor);
			unf->SetMarkerColor(BeamOnColor);
			unf->SetMarkerStyle(20);
			unf->SetMarkerSize(1.);	

			// Draw the data points first to get the beautiful canvas 
			PlotCanvas->cd();
			if (ClosureTest == true) { unf->Draw("p0 hist"); }
			else { unf->Draw("e1x0"); }	

			TLatex *textSlice = new TLatex();
			textSlice->SetTextFont(FontStyle);
			textSlice->SetTextSize(0.06);
			TString PlotNameDuplicate = PlotNames[WhichPlot];
			TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("Reco","") ;
			textSlice->DrawLatexNDC(0.24, 0.8, LatexLabel[ReducedPlotName]);				

			// -------------------------------------- //	

			// The MC CC1p prediction has to be multiplied by the additional smearing matrix Ac

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
			PlotCanvas->cd();					
			TrueUnf->Draw("hist same");		

			// -------------------------------------- //						

			// Plotting again so that the data points are on top 
			PlotCanvas->cd();
			if (ClosureTest == true) { unf->Draw("p0 hist same"); }
			else { 

				// Draw the truth w/o the additional smearing
				//PlotsTrue[4][WhichPlot]->SetLineColor(kMagenta);
				//PlotsTrue[4][WhichPlot]->Draw("hist same");				
				
				unf->Draw("e1x0 same"); 

				unfStat->SetLineColor(kBlack);
				unfStat->Draw("e1x0 same");			

				unfShapeOnly->SetLineColor(kOrange+7);
				//unfShapeOnly->Draw("e1x0 same");

				unfNormOnly->SetFillColorAlpha(kGray+1, 0.45);
				unfNormOnly->SetLineColor(kGray+1);
				//unfNormOnly->SetFillStyle(3000);
				unfNormOnly->Draw("e2 hist same");			

				gPad->RedrawAxis();		
		
			}

			// ------------------------------------------------------------------------------

			// Legend & POT Normalization

			double tor860_wcut = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);
			TString Label = ToStringPOT(tor860_wcut)+" POT";

			TLegend* legData = new TLegend(0.23,0.89,0.95,0.98);
			legData->SetBorderSize(0);
			legData->SetTextSize(0.06);
			legData->SetTextFont(FontStyle);
			legData->SetNColumns(1);
			legData->SetMargin(0.1);			

//			legData->AddEntry(unf,"MicroBooNE Data " + Runs[WhichRun] + " " + Label,"ep");
			legData->AddEntry(unf,"MicroBooNE Data " + Label,"ep");

			TLegendEntry* lMC = legData->AddEntry(TrueUnf,"MC uB Tune","l");
			lMC->SetTextColor(OverlayColor);	

			legData->Draw();
			
			TString CanvasPath = PlotPath+NameOfSamples[0];
			TString FullCanvasName = "/"+Tune+"WienerSVD_XSections_"+CanvasName+"_"+UBCodeVersion+Subtract+".pdf";
			if (ClosureTest == true) { FullCanvasName = "/ClosureTest_"+Tune+"WienerSVD_XSections_"+CanvasName+"_"+UBCodeVersion+Subtract+".pdf"; }
			PlotCanvas->SaveAs(CanvasPath+FullCanvasName);	
			delete PlotCanvas;			

			// ----------------------------------------------------------------------------------------------------------------		

			TH1D* diff = new TH1D("diff_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],"Fractional difference of unf and signal model",n, Nuedges);
			
			for(int i=1; i <= n; i++) {
			
				double s2s = 0;
				double u = unf->GetBinContent(i);
				double t = signal(i-1);
				if (t != 0) { s2s = u/t - 1; }
				else { s2s = 1.; } // attention to this t = 0
				diff->SetBinContent(i, s2s); // in percentage 
				
			}	

			// ---------------------------------------------------------------------------------------------------------------------------

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

			// ---------------------------------------------------------------------------------------------------------------------------

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
			M2H(UnfoldCov, unfcov);
			M2H( NormShapeVector[0] , normunfcov);
			M2H( (NormShapeVector[1] + UnfStatCov), shapeunfcov);	
			M2H(CovRotation, covrot);
			M2H(CovRotation_T, transpcovrot);											

			// ---------------------------------------------------------------------------------------------------------------------------
    
			if (ClosureTest == false) {

				ExtractedXSec->cd();

				unf->Write("Reco"+PlotNames[WhichPlot]);
				unfShapeOnly->Write("ShapeOnlyReco"+PlotNames[WhichPlot]);
				unfNormOnly->Write("NormOnlyReco"+PlotNames[WhichPlot]);

				unfMCStat->Write("MCStatReco"+PlotNames[WhichPlot]);
				unfStat->Write("StatReco"+PlotNames[WhichPlot]);
				unfLY->Write("LYReco"+PlotNames[WhichPlot]);
				unfTPC->Write("TPCReco"+PlotNames[WhichPlot]);
				unfSCERecomb2->Write("SCERecomb2Reco"+PlotNames[WhichPlot]);
				unfFlux->Write("FluxReco"+PlotNames[WhichPlot]);
				unfXSec->Write("XSecReco"+PlotNames[WhichPlot]);
				unfG4->Write("G4Reco"+PlotNames[WhichPlot]);
				unfDirt->Write("DirtReco"+PlotNames[WhichPlot]);
				unfPOT->Write("POTReco"+PlotNames[WhichPlot]);
				//unfNuWro->Write("NuWroReco"+PlotNames[WhichPlot]);				
				unfNTarget->Write("NTargetReco"+PlotNames[WhichPlot]);

				unfFullUnc->Write("RecoFullUnc"+PlotNames[WhichPlot]);				
				TrueUnf->Write("True"+PlotNames[WhichPlot]);
				QETrueUnf->Write("QETrue"+PlotNames[WhichPlot]);
				MECTrueUnf->Write("MECTrue"+PlotNames[WhichPlot]);
				RESTrueUnf->Write("RESTrue"+PlotNames[WhichPlot]);
				DISTrueUnf->Write("DISTrue"+PlotNames[WhichPlot]);				
				COHTrueUnf->Write("COHTrue"+PlotNames[WhichPlot]);				
				PlotsTrue[4][WhichPlot]->Write("NoSmearTrue"+PlotNames[WhichPlot]);
				QEPlotsTrue[4][WhichPlot]->Write("QENoSmearTrue"+PlotNames[WhichPlot]);
				MECPlotsTrue[4][WhichPlot]->Write("MECNoSmearTrue"+PlotNames[WhichPlot]);
				RESPlotsTrue[4][WhichPlot]->Write("RESNoSmearTrue"+PlotNames[WhichPlot]);
				DISPlotsTrue[4][WhichPlot]->Write("DISNoSmearTrue"+PlotNames[WhichPlot]);
				COHPlotsTrue[4][WhichPlot]->Write("COHNoSmearTrue"+PlotNames[WhichPlot]);																			
				smear->Write("Ac"+PlotNames[WhichPlot]);
				unfcov->Write("UnfCov"+PlotNames[WhichPlot]);
				shapeunfcov->Write("ShapeUnfCov"+PlotNames[WhichPlot]); // Shape systematic unc
				normunfcov->Write("NormUnfCov"+PlotNames[WhichPlot]); // Norm systematic unc			
				CovarianceMatrices[WhichPlot]->Write("Cov"+PlotNames[WhichPlot]); // covariance before unfolding			
				//wiener->Write("Wiener"+PlotNames[WhichPlot]);
				//diff->Write("Diff"+PlotNames[WhichPlot]);
				//bias->Write("Bias"+PlotNames[WhichPlot]);
				//bias2->Write("Bias2"+PlotNames[WhichPlot]);
				//fracError->Write("FracErr"+PlotNames[WhichPlot]);
				//absError->Write("AbsErr"+PlotNames[WhichPlot]);
				//MSE->Write("MSE"+PlotNames[WhichPlot]);
				//MSE2->Write("MSE2"+PlotNames[WhichPlot]);
				ResponseMatrices[WhichPlot]->Write("Response"+PlotNames[WhichPlot]);
				covrot->Write("CovRot"+PlotNames[WhichPlot]);
				transpcovrot->Write("TranspCovRot"+PlotNames[WhichPlot]);								

				// ---------------------------------------------------------------------------------------------------------------------------

				// Make the additional smearing matrix pretty

				TString SmearCanvasName = "Smear_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun];
				TCanvas* SmearPlotCanvas = new TCanvas(SmearCanvasName,SmearCanvasName,205,34,1024,768);
				SmearPlotCanvas->cd();
				SmearPlotCanvas->SetBottomMargin(0.17);
				SmearPlotCanvas->SetTopMargin(0.13);			
				SmearPlotCanvas->SetLeftMargin(0.2);
				SmearPlotCanvas->SetRightMargin(0.2);

				smear->GetXaxis()->CenterTitle();
				smear->GetXaxis()->SetLabelFont(FontStyle);
				smear->GetXaxis()->SetTitleFont(FontStyle);
				smear->GetXaxis()->SetLabelSize(TextSize);
				smear->GetXaxis()->SetTitleSize(TextSize);
				smear->GetXaxis()->SetNdivisions(6);
				TString Xtitle = smear->GetXaxis()->GetTitle();
				smear->GetXaxis()->SetTitle("True " + Xtitle);

				smear->GetYaxis()->CenterTitle();
				smear->GetYaxis()->SetLabelFont(FontStyle);
				smear->GetYaxis()->SetTitleFont(FontStyle);
				smear->GetYaxis()->SetLabelSize(TextSize);
				smear->GetYaxis()->SetTitleSize(TextSize);
				smear->GetYaxis()->SetNdivisions(6);
				TString Ytitle = smear->GetYaxis()->GetTitle();
				smear->GetYaxis()->SetTitle("Reco " + Ytitle);					

				smear->GetZaxis()->SetLabelFont(FontStyle);
				smear->GetZaxis()->SetTitleFont(FontStyle);
				smear->GetZaxis()->SetLabelSize(TextSize);
				smear->GetZaxis()->SetTitleSize(TextSize);
				smear->GetZaxis()->SetRangeUser(-0.1,1.);				

				gStyle->SetPaintTextFormat("4.2f");
				smear->SetMarkerColor(kWhite);
				smear->Draw("coltz text");

				TLatex *text = new TLatex();
				text->SetTextFont(FontStyle);
				text->SetTextSize(0.06);
				text->DrawTextNDC(0.3, 0.92, Runs[WhichRun] + ", " + LatexLabel[PlotNames[WhichPlot]]);

				TString SmearCanvas = "/Smear_"+Tune+"WienerSVD_XSections_"+CanvasName+"_"+UBCodeVersion+Subtract+".pdf";
				SmearPlotCanvas->SaveAs(CanvasPath+SmearCanvas);
				delete SmearPlotCanvas;

			}

			// ---------------------------------------------------------------------------------------------------------------------------

		} // End of the loop over the plots

		if (ClosureTest == false ) {

			ExtractedXSec->Close();
			std::cout << std::endl << "File " << NameExtractedXSec << " created" << std::endl << std::endl;

		}

//		for (int i = 0; i < NSamples; i++) { FileSample[i]->Close(); }

	} // End of the loop over the runs	

} // End of the program 
