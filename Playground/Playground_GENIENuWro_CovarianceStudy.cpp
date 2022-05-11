#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TLatex.h>
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

//----------------------------------------//

TString ToStringPOT(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

//----------------------------------------//

void ReweightXSec(TH1D* h, double SF = 1.) {

	int NBins = h->GetXaxis()->GetNbins();

	for (int i = 0; i < NBins; i++) {

		double CurrentEntry = h->GetBinContent(i+1);
		double NewEntry = CurrentEntry * SF / h->GetBinWidth(i+1);

		double CurrentError = h->GetBinError(i+1);
		double NewError = CurrentError * SF / h->GetBinWidth(i+1);

		h->SetBinContent(i+1,NewEntry); 
		h->SetBinError(i+1,NewError); 

	}

}

//----------------------------------------//

void Playground_GENIENuWro_CovarianceStudy() {

	//----------------------------------------//

	int LineWidth = 2;

	gStyle->SetPalette(55); 
	const Int_t NCont = 999; 
	gStyle->SetNumberContours(NCont); 
	gStyle->SetTitleSize(0.07,"t");

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();	
	gStyle->SetOptStat(0);
	gStyle->SetEndErrorSize(4);	

	TString CutExtension = "_NoCuts_PID_NuScore";

	//----------------------------------------//

	vector<TString> PlotNames;
	PlotNames.push_back("DeltaPTPlot");	
	PlotNames.push_back("DeltaAlphaTPlot");
//	PlotNames.push_back("DeltaPtxPlot");
//	PlotNames.push_back("DeltaPtyPlot");		 

	const int N1DPlots = PlotNames.size();

	//----------------------------------------//

	// CV Flux File

	TFile* FluxFile = TFile::Open("../MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));

	//----------------------------------------//	

	// open the file that contains the unfolding uncertainty

	TString NameUnfUnc = PathToExtractedXSec+"WienerSVD_UnfoldingUnc_Combined_"+UBCodeVersion+".root";
	TFile* FileUnfUnc = TFile::Open(NameUnfUnc,"readonly");			

	//----------------------------------------//
	
	double DataPOT = PeLEE_ReturnBeamOnRunPOT("Combined");					
	double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface);	
			
	//----------------------------------------//

	vector<TCanvas*> PlotCanvas; PlotCanvas.clear();

	vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();
	vector<vector<TH1D*> > PlotsCVTrue; PlotsCVTrue.clear();

	vector<TH1D*> UnfUnc; UnfUnc.clear();
	vector<TH2D*> ResponseMatrices; ResponseMatrices.clear();
	vector<TH2D*> FullCovarianceMatrices; FullCovarianceMatrices.clear();
	vector<TH2D*> PartialCovarianceMatrices; PartialCovarianceMatrices.clear();	
	vector<TH2D*> MCStatCovarianceMatrices; MCStatCovarianceMatrices.clear();
	vector<TH2D*> StatCovarianceMatrices; StatCovarianceMatrices.clear();
	vector<TH2D*> XSecCovarianceMatrices; XSecCovarianceMatrices.clear();						

	//----------------------------------------//

	TString FileResponseName = MigrationMatrixPath+"FileResponseMatrices_Overlay9_Combined_"+UBCodeVersion+".root";
	cout << "File Responses = " << FileResponseName << endl;
	TFile* FileResponseMatrices = new TFile(FileResponseName,"readonly");

	//----------------------------------------//		

	TString CVFileCovarianceName = MigrationMatrixPath+"WienerSVD_Total_CovarianceMatrices_Overlay9_Combined_"+UBCodeVersion+".root";
	cout << "File CV Covariances = " << CVFileCovarianceName << endl;			
	TFile* CVFileCovarianceMatrices = new TFile(CVFileCovarianceName,"readonly");

	//----------------------------------------//		

	TString AltCVFileCovarianceName = MigrationMatrixPath+"Overlay9NuWroWienerSVD_Total_CovarianceMatrices_Overlay9_Combined_"+UBCodeVersion+".root";
	cout << "File Alt CV Covariances = " << AltCVFileCovarianceName << endl;			
	TFile* AltCVFileCovarianceMatrices = new TFile(AltCVFileCovarianceName,"readonly");	

	//----------------------------------------//
		
	TString PathToFilesUBCodeExtension = PathToFiles+CutExtension;

	TString NuWroFileName = "STVStudies_Overlay9NuWro_Combined"+CutExtension+".root";
	TFile* NuWroFile = new TFile(PathToFilesUBCodeExtension+"/"+NuWroFileName,"readonly");

	TString TrueGENIEFileName = "TruthSTVAnalysis_Overlay9_Combined_"+UBCodeVersion+".root";
	TFile* TrueGENIEFile = new TFile(PathToFiles+"/"+TrueGENIEFileName,"readonly");

	TString TrueNuWroFileName = "TruthSTVAnalysis_Overlay9NuWro_Combined_"+UBCodeVersion+".root";
	TFile* TrueNuWroFile = new TFile(PathToFiles+"/"+TrueNuWroFileName,"readonly");			

	//----------------------------------------//

	// Loop over the plots

	for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

		//----------------------------------------//			

		UnfUnc.push_back((TH1D*)FileUnfUnc->Get("UnfUnc_"+PlotNames[WhichPlot]));			

		//----------------------------------------//

		ResponseMatrices.push_back((TH2D*)FileResponseMatrices->Get("POTScaledCC1pReco"+PlotNames[WhichPlot]+"2D"));

		// Already flux-averaged rates
		FullCovarianceMatrices.push_back((TH2D*)CVFileCovarianceMatrices->Get("TotalCovariance_"+PlotNames[WhichPlot]));
		MCStatCovarianceMatrices.push_back((TH2D*)CVFileCovarianceMatrices->Get("MCStatCovariance_"+PlotNames[WhichPlot]));
		StatCovarianceMatrices.push_back((TH2D*)CVFileCovarianceMatrices->Get("StatCovariance_"+PlotNames[WhichPlot]));	
		XSecCovarianceMatrices.push_back((TH2D*)CVFileCovarianceMatrices->Get("XSecCovariance_"+PlotNames[WhichPlot]));						
		PartialCovarianceMatrices.push_back((TH2D*)AltCVFileCovarianceMatrices->Get("TotalCovariance_"+PlotNames[WhichPlot]));					

		//----------------------------------------//

		TH1D* DataPlot = (TH1D*)(NuWroFile->Get("CC1pReco"+PlotNames[WhichPlot]));
		TH1D* NuWroTrue = (TH1D*)(TrueNuWroFile->Get("True"+PlotNames[WhichPlot]));
		TH1D* GENIETrue = (TH1D*)(TrueGENIEFile->Get("True"+PlotNames[WhichPlot]));
	
		//----------------------------------------//

		int n = GENIETrue->GetNbinsX();
		double Nuedges[n+1];
		double NuedgesOffset[n+1];		
			    
		for (int i = 0; i < n+1; i++) { 
			
			Nuedges[i] = GENIETrue->GetBinLowEdge(i+1);
			NuedgesOffset[i] = GENIETrue->GetBinLowEdge(i+1) + 0.2*GENIETrue->GetBinWidth(i+1);			 
			
		}

		//----------------------------------------//

		int m = DataPlot->GetNbinsX();			
		TString XTitle = DataPlot->GetXaxis()->GetTitle();
		TString YTitle = DataPlot->GetYaxis()->GetTitle();	

		// Flux-averaged event rates
		double SF = Units/(IntegratedFlux*NTargets);
		DataPlot->Scale(SF);
		GENIETrue->Scale(SF); // CV
		NuWroTrue->Scale(SF); // Alt CV		
			
		//----------------------------------------//

		// Construct vectors (for 1D histogram) and matrices (for 2D histogram) for input
		// Convert input into mathematical formats, easy and clean to be processed. 
		// Converted defined/implemented in source files, see include/Util.h

		TVectorD signal(n);			
		TVectorD altsignal(n);					
		TVectorD measure(m);
		TMatrixD response(m, n);
		TMatrixD fullcovariance(m, m);		
		TMatrixD partialcovariance(m, m);
		TMatrixD mcstatcovariance(m, m);			
		TMatrixD statcovariance(m, m);
		TMatrixD xseccovariance(m, m);

		H2V(GENIETrue, signal);
		H2V(NuWroTrue, altsignal);
		H2V(DataPlot, measure);
		H2M(ResponseMatrices[WhichPlot], response, kFALSE); // X axis: Reco, Y axis: True
		H2M(FullCovarianceMatrices[WhichPlot], fullcovariance, kTRUE); // X axis: True, Y axis: Reco
		H2M(PartialCovarianceMatrices[WhichPlot], partialcovariance, kTRUE);		
		H2M(MCStatCovarianceMatrices[WhichPlot], mcstatcovariance, kTRUE);
		H2M(StatCovarianceMatrices[WhichPlot], statcovariance, kTRUE);	
		H2M(XSecCovarianceMatrices[WhichPlot], xseccovariance, kTRUE);

		//----------------------------------------//

		// Construct to record additinal smearing matrix and wiener filter (diagomal matrix) elements. 
    
		TMatrixD FullAddSmear(n,n);
		TVectorD FullWF(n);
		TMatrixD FullUnfoldCov(n,n);
		TMatrixD FullCovRotation(n,n);

		TMatrixD PartialAddSmear(n,n);
		TVectorD PartialWF(n);
		TMatrixD PartialUnfoldCov(n,n);
		TMatrixD PartialCovRotation(n,n);

		//----------------------------------------//

		TH2D* fullunfcov = new TH2D("fullunfcov_"+PlotNames[WhichPlot]+"_Combined","Unfolded spectrum covariance", n, Nuedges, n, Nuedges);

		TH2D* partialunfcov = new TH2D("partialunfcov_"+PlotNames[WhichPlot]+"_Combined","Unfolded spectrum covariance", n, Nuedges, n, Nuedges);			

		//----------------------------------------//

		// Core implementation of Wiener-SVD
		// AddSmear and WF to record the core information in the unfolding.

		TVectorD fullunfold = WienerSVD(response, signal, measure, fullcovariance, 2, 0.5, FullAddSmear, FullWF, FullUnfoldCov, FullCovRotation);			
		TVectorD partialunfold = WienerSVD(response, signal, measure, partialcovariance, 2, 0.5, PartialAddSmear, PartialWF, PartialUnfoldCov, PartialCovRotation);

		//----------------------------------------//

		TString CanvasName = PlotNames[WhichPlot]+"_Combined";
		TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
		PlotCanvas->cd();
		PlotCanvas->SetBottomMargin(0.16);
		PlotCanvas->SetTopMargin(0.13);			
		PlotCanvas->SetLeftMargin(0.21);			
		
		TH1D* fullunf = new TH1D("fullunf_"+PlotNames[WhichPlot]+"_Combined",";"+XTitle+";"+YTitle,n,Nuedges);
		TH1D* partialunf = new TH1D("partialunf_"+PlotNames[WhichPlot]+"_Combined",";"+XTitle+";"+YTitle,n,NuedgesOffset);	

		//------------------------------//

		// Legend & POT Normalization

		double tor860_wcut = PeLEE_ReturnBeamOnRunPOT("Combined");
		TString Label = ToStringPOT(tor860_wcut)+" POT";

		TLegend* leg = new TLegend(0.23,0.89,0.95,0.98);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.05);
		leg->SetTextFont(FontStyle);
		leg->SetNColumns(2);
		leg->SetMargin(0.1);				

		//----------------------------------------//

		V2H(fullunfold, fullunf); 	
		M2H(FullUnfoldCov, fullunfcov);	

		V2H(partialunfold, partialunf); 	
		M2H(PartialUnfoldCov, partialunfcov);	

		TMatrixD FullCovRotation_T (TMatrixD::kTransposed, FullCovRotation); 
		TMatrixD FullMCStatCov = FullCovRotation*mcstatcovariance*FullCovRotation_T;
		TMatrixD FullStatCov = FullCovRotation*statcovariance*FullCovRotation_T;	
		TMatrixD FullXSecCov = FullCovRotation*xseccovariance*FullCovRotation_T;					 						

		//----------------------------------------//	

		ReweightXSec(fullunf);
		ReweightXSec(partialunf);							

		for (int i = 1; i <= n;i++ ) { 

			double BinUnfUnc = UnfUnc[WhichPlot]->GetBinContent(i);	
			double BinWidth = fullunf->GetBinWidth(i);				

			// Using the full covariance 
			// and separating out the Stat/MC Stat/XSec unc
			double FullUnc = TMath::Sqrt( FullStatCov(i-1,i-1) + FullMCStatCov(i-1,i-1) ) / BinWidth;
			double FullCovUnc = TMath::Sqrt( TMath::Power(FullUnc,2.) + TMath::Power(BinUnfUnc,2.) );
			fullunf->SetBinError(i, FullCovUnc );

			// Using the full covariance 
			// and separating out the Stat/MC Stat/XSec unc
			double PartialUnc = TMath::Sqrt( PartialUnfoldCov(i-1,i-1) ) / BinWidth;			
			double PartialCovUnc = TMath::Sqrt( TMath::Power(PartialUnc,2.)  + TMath::Power(BinUnfUnc,2.) );
			partialunf->SetBinError(i, PartialCovUnc );											

		}					

		//----------------------------------------//			

		fullunf->GetXaxis()->CenterTitle();
		fullunf->GetXaxis()->SetLabelSize(TextSize);
		fullunf->GetXaxis()->SetLabelFont(FontStyle);
		fullunf->GetXaxis()->SetTitleSize(TextSize);
		fullunf->GetXaxis()->SetTitleFont(FontStyle);			
		fullunf->GetXaxis()->SetNdivisions(6);	

		fullunf->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);
		fullunf->GetYaxis()->CenterTitle();
		fullunf->GetYaxis()->SetLabelSize(TextSize);
		fullunf->GetYaxis()->SetLabelFont(FontStyle);
		fullunf->GetYaxis()->SetTitleSize(TextSize);
		fullunf->GetYaxis()->SetTitleFont(FontStyle);			
		fullunf->GetYaxis()->SetNdivisions(6);

		fullunf->SetLineColor(OverlayColor);
		fullunf->SetLineWidth(LineWidth);		
		fullunf->SetMarkerColor(OverlayColor);
		fullunf->SetMarkerStyle(20);
		fullunf->SetMarkerSize(1.);

		partialunf->SetLineColor(kOrange+7);
		partialunf->SetLineWidth(LineWidth);		
		partialunf->SetMarkerColor(kOrange+7);
		partialunf->SetMarkerStyle(20);
		partialunf->SetMarkerSize(1.);			

		// ----------------------------------------//		

		// CV GENIE has to be multiplied by the additional smearing matrix Ac

		TH1D* FullTrueGENIE = new TH1D("FullTrueGENIE_"+PlotNames[WhichPlot]+"_Combined",";"+XTitle+";"+YTitle,n,Nuedges);
		TVectorD FullAcTrueGENIEUnfold = FullAddSmear * signal;
		V2H(FullAcTrueGENIEUnfold, FullTrueGENIE);		
		ReweightXSec(FullTrueGENIE);
		FullTrueGENIE->SetLineColor(OverlayColor);
		FullTrueGENIE->SetLineWidth(LineWidth);					
		FullTrueGENIE->SetLineStyle(kDashed);

		TH1D* PartialTrueGENIE = new TH1D("PartialTrueGENIE_"+PlotNames[WhichPlot]+"_Combined",";"+XTitle+";"+YTitle,n,Nuedges);
		TVectorD PartialAcTrueGENIEUnfold = PartialAddSmear * signal;
		V2H(PartialAcTrueGENIEUnfold, PartialTrueGENIE);		
		ReweightXSec(PartialTrueGENIE);
		PartialTrueGENIE->SetLineColor(kOrange+7);	
		PartialTrueGENIE->SetLineStyle(kDashed);	
		PartialTrueGENIE->SetLineWidth(LineWidth);						

		// ----------------------------------------//

		// NuWro has to be multiplied by the additional smearing matrix Ac

		TH1D* FullTrueNuWro = new TH1D("FullTrueNuWro_"+PlotNames[WhichPlot]+"_Combined",";"+XTitle+";"+YTitle,n,Nuedges);
		TVectorD FullAcTrueNuWroUnfold = FullAddSmear * altsignal;
		V2H(FullAcTrueNuWroUnfold, FullTrueNuWro);														
		ReweightXSec(FullTrueNuWro);															
		FullTrueNuWro->SetLineColor(OverlayColor);
		FullTrueNuWro->SetMarkerColor(OverlayColor);
		FullTrueNuWro->SetLineWidth(3);		

		TH1D* PartialTrueNuWro = new TH1D("PartialTrueNuWro_"+PlotNames[WhichPlot]+"_Combined",";"+XTitle+";"+YTitle,n,Nuedges);
		TVectorD PartialAcTrueNuWroUnfold = PartialAddSmear * altsignal;
		V2H(PartialAcTrueNuWroUnfold, PartialTrueNuWro);														
		ReweightXSec(PartialTrueNuWro);															
		PartialTrueNuWro->SetLineColor(kOrange+7);
		PartialTrueNuWro->SetMarkerColor(kOrange+7);	
		PartialTrueNuWro->SetLineWidth(3);				

		// ----------------------------------------//	

		double max = 1.2 * TMath::Max(fullunf->GetMaximum(), partialunf->GetMaximum());		
		fullunf->GetYaxis()->SetRangeUser(0.,max);			

		//------------------------------//					

		PlotCanvas->cd();
		fullunf->Draw("e1x0"); // using full covariance and separating out the stat/mc stat/xsec unc
		leg->AddEntry(fullunf,"Unf (Full Cov)","lep");
		FullTrueNuWro->Draw("same hist");
		leg->AddEntry(FullTrueNuWro,"True (Full Cov)","l");		


		partialunf->Draw("e1x0 same"); // using full covariance and separating out the stat/mc stat/xsec unc
		leg->AddEntry(partialunf,"Unf (Part Cov)","lep");		
		PartialTrueNuWro->Draw("same hist");	
		leg->AddEntry(PartialTrueNuWro,"True (Part Cov)","l");	

		fullunf->Draw("e1x0 same");	
		partialunf->Draw("e1x0 same");	

		leg->Draw();

		//------------------------------//

		TLatex *text = new TLatex();
		text->SetTextFont(FontStyle);
		text->SetTextSize(TextSize);
		text->DrawLatexNDC(0.24, 0.8, "NuWro");

		//------------------------------//			

/*
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
						double UnfUncBin = UnfUnc[WhichPlot]->GetBinContent(ix);
//						double UnfUncBin = 0.;

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
*/
		//------------------------------//

	} // End of the loop over the plots

} // End of the program 