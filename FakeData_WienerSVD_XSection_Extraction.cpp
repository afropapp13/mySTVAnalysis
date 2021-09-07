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

	// We want the number of events, as if the bin width is 1
	double ExtraFactor = 1.;
	if (NBins == 1) { ExtraFactor = 2.; }

	for (int i = 0; i < NBins; i++) {

		double CurrentEntry = h->GetBinContent(i+1);
		double NewEntry = CurrentEntry * SF / h->GetBinWidth(i+1) * ExtraFactor;

		double CurrentError = h->GetBinError(i+1);
		double NewError = CurrentError * SF / h->GetBinWidth(i+1) * ExtraFactor;

		h->SetBinContent(i+1,NewEntry); 
		h->SetBinError(i+1,NewError); 
//		h->SetBinError(i+1,0.000001); 

	}

}

// -------------------------------------------------------------------------------------------------------------------------------------

void FakeData_WienerSVD_XSection_Extraction(TString OverlaySample = "", TString BeamOnSample = "") {

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

	NameOfSamples.push_back("Overlay9"); if (OverlaySample != "") { NameOfSamples[0] = OverlaySample; }
	NameOfSamples.push_back("BeamOn9");  if (BeamOnSample != "") { NameOfSamples[1] = BeamOnSample; }
	NameOfSamples.push_back("ExtBNB9"); 
	NameOfSamples.push_back("OverlayDirt9");
	NameOfSamples.push_back("GenieOverlay");	 

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

	// -------------------------------------------------------------------------------------------------------------------------------------

	// CV Flux File

	TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));

	// -------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------
	
		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);					
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface);	
			
		// -------------------------------------------------------------------------------------		

		vector<TCanvas*> PlotCanvas; PlotCanvas.clear();

		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
		vector<vector<TH1D*> > PlotsBkgReco; PlotsBkgReco.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();

		vector<TH2D*> ResponseMatrices; ResponseMatrices.clear();
		vector<TH2D*> CovarianceMatrices; CovarianceMatrices.clear();

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString FileResponseName = MigrationMatrixPath+"FileResponseMatrices_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
		TFile* FileResponseMatrices = new TFile(FileResponseName,"readonly");

		TString FileCovarianceName = MigrationMatrixPath+BeamOnSample+"WienerSVD_Total_CovarianceMatrices_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
		TFile* FileCovarianceMatrices = new TFile(FileCovarianceName,"readonly");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();
		
		TString PathToFilesUBCodeExtension = PathToFiles+CutExtension;

		// Store the extracted xsections & associated files in dedicated file

		TString NameExtractedXSec = PathToExtractedXSec+BeamOnSample+"WienerSVD_ExtractedXSec_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+Subtract+".root";
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
			
			TString FileName = "STVStudies_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+CutExtension+".root";
			FileSample.push_back(TFile::Open(PathToFilesUBCodeExtension+"/"+FileName)); 
			
			}

			if (NameOfSamples[WhichSample] == "GenieOverlay") { 
			
			TString FileName = "TruthSTVAnalysis_"+BeamOnSample+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
			FileSample.push_back(TFile::Open(PathToFiles+FileName));  
			
			}

			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();
			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();
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
		
			} // End of the loop over the plots

			PlotsReco.push_back(CurrentPlotsReco);		
			PlotsTrue.push_back(CurrentPlotsTrue);		
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

			// -----------------------------------------------------------------------------------------------------

			// True CC1p Signal MC // No detector/ reconstruction / smearing effects

			int n = PlotsTrue[4][WhichPlot]->GetNbinsX();
			double Nuedges[n+1];
			    
			for (int i = 0; i < n+1; i++) { Nuedges[i] = PlotsTrue[4][WhichPlot]->GetBinLowEdge(i+1); }

			// -------------------------------------------------------------------------------------------------

			// BeamOn = Alternative model for fake data studies in this case
			// Thus we grab the CC1p part, without any subtractions

			TH1D* DataPlot = (TH1D*)(PlotsCC1pReco[1][WhichPlot]->Clone());

			int m = DataPlot->GetNbinsX();			
			TString XTitle = DataPlot->GetXaxis()->GetTitle();
			TString YTitle = DataPlot->GetYaxis()->GetTitle();	

			// Flux-averaged event rates
			DataPlot->Scale(Units/(IntegratedFlux*NTargets));
			PlotsTrue[4][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));		 

			// -------------------------------------------------------------------------------------------

			// Construct vectors (for 1D histogram) and matrices (for 2D histogram) for input

			TVectorD signal(n);
			TVectorD measure(m);
			TMatrixD response(m, n);
			TMatrixD covariance(m, m);
			TMatrixD statcovariance(m, m);
			TMatrixD systcovariance(m, m);

			// Convert input into mathematical formats, easy and clean to be processed. 
			// Converted defined/implemented in source files, see include/Util.h

			H2V(PlotsTrue[4][WhichPlot], signal);
			H2V(DataPlot, measure);
			H2M(ResponseMatrices[WhichPlot], response, kFALSE); // X axis: Reco, Y axis: True
			H2M(CovarianceMatrices[WhichPlot], covariance, kTRUE); // X axis: True, Y axis: Reco

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

			// --------------------------------------------------------------------------------------------------				

			// Scaling by correct units & bin width division

			for (int i = 1; i <= n;i++ ) { 

			double CVInBin = unf->GetBinContent(i);

			// default / total uncertainty
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

			unf->GetYaxis()->SetRangeUser(XSecRange[PlotNames[WhichPlot]].first,XSecRange[PlotNames[WhichPlot]].second);
			unf->SetLineColor(BeamOnColor);
			unf->SetMarkerColor(BeamOnColor);
			unf->SetMarkerStyle(20);
			unf->SetMarkerSize(1.);	

			// Draw the data points first to get the beautiful canvas 
			PlotCanvas->cd();
			unf->Draw("e1x0");			

			// The MC CC1p prediction has to be multiplied by the additional smearing matrix Ac

			TH1D* TrueUnf = new TH1D("TrueUnf_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun],";"+XTitle+";"+YTitle,n,Nuedges);
			TVectorD AcTrueUnfold = AddSmear * signal;
			V2H(AcTrueUnfold, TrueUnf);

			ReweightXSec(TrueUnf);
//			TrueUnf->Scale(Units/(IntegratedFlux*NTargets));
			TrueUnf->SetLineColor(OverlayColor);
			TrueUnf->SetMarkerColor(OverlayColor);	
			PlotCanvas->cd();				
			TrueUnf->Draw("hist same");

			// Plotting again so that the data points are on top 
			PlotCanvas->cd();
			unf->Draw("e1x0 same"); 

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

			TString ReducedLabel = 	BeamOnSample;
			TString PrintLabel = ReducedLabel.ReplaceAll("Overlay9","");	

			legData->AddEntry(unf,"Unfolded " + PrintLabel,"ep");

			TLegendEntry* lMC = legData->AddEntry(TrueUnf,"True " + PrintLabel,"l");
			lMC->SetTextColor(OverlayColor);	

			legData->Draw();
			
			TString CanvasPath = PlotPath+NameOfSamples[0];
			TString FullCanvasName = "/"+BeamOnSample+"WienerSVD_XSections_"+CanvasName+"_"+UBCodeVersion+Subtract+".pdf";
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

			// ---------------------------------------------------------------------------------------------------------------------------

			ExtractedXSec->cd();

			unf->Write("Reco"+PlotNames[WhichPlot]);
			TrueUnf->Write("True"+PlotNames[WhichPlot]);			
			smear->Write("Ac"+PlotNames[WhichPlot]);
			unfcov->Write("UnfCov"+PlotNames[WhichPlot]);	
			CovarianceMatrices[WhichPlot]->Write("Cov"+PlotNames[WhichPlot]);				
			//wiener->Write("Wiener"+PlotNames[WhichPlot]);
			//diff->Write("Diff"+PlotNames[WhichPlot]);
			//bias->Write("Bias"+PlotNames[WhichPlot]);
			//bias2->Write("Bias2"+PlotNames[WhichPlot]);
			//fracError->Write("FracErr"+PlotNames[WhichPlot]);
			//absError->Write("AbsErr"+PlotNames[WhichPlot]);
			//MSE->Write("MSE"+PlotNames[WhichPlot]);
			//MSE2->Write("MSE2"+PlotNames[WhichPlot]);
			ResponseMatrices[WhichPlot]->Write("Response"+PlotNames[WhichPlot]);	

			// ---------------------------------------------------------------------------------------------------------------------------

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
			SmearPlotCanvas->SaveAs(CanvasPath+SmearCanvas);
			delete SmearPlotCanvas;

			// ---------------------------------------------------------------------------------------------------------------------------

		} // End of the loop over the plots

		ExtractedXSec->Close();
		std::cout << std::endl << "File " << NameExtractedXSec << " created" << std::endl << std::endl;

	} // End of the loop over the runs	

} // End of the program 
