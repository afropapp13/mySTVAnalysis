#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TSpline.h>
#include <TProfile.h>
#include <TAttFill.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

#include "../myClasses/Tools.h"
#include "../myClasses/Constants.h"

using namespace Constants;
using namespace std;

// ---------------------------------------- //

int LocateBinWithValue(TH1D* h, double Value) {

	int NBins = h->GetXaxis()->GetNbins();

	for (int i = 1; i <= NBins; i++) {

		double CurrentEntry = h->GetBinContent(i);
		if (CurrentEntry == Value) { return i; } 

	}

	return -99;

}

// ---------------------------------------- //

void PrettyPlot(TH1D* h,int LineWidth = 2, int FontStyle = 132, int Ndivisions = 6, double TextSize = 0.06) {

	// ---------------------------------------- //

	h->SetLineWidth(LineWidth);

	// ---------------------------------------- //

	// X-axis

	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetLabelSize(TextSize);
	h->GetXaxis()->SetTitleSize(TextSize);
	h->GetXaxis()->SetTitleOffset(1.05);
	h->GetXaxis()->SetNdivisions(Ndivisions);

	// ---------------------------------------- //

	// Y-axis

	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleSize(TextSize); 
	h->GetYaxis()->SetTickSize(0.02);
	h->GetYaxis()->SetLabelSize(TextSize);
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleOffset(1.05);
	h->GetYaxis()->SetNdivisions(Ndivisions);
	//h->GetYaxis()->SetTitle("POT Normalized Events");

	return;	

}

// ---------------------------------------- //

void MCC8_OverlayXSecMethods() {

	// ---------------------------------------- //

	int LineWidth = 3;
	int FontStyle = 132;
	int Ndivisions = 6; 
	double TextSize = 0.06; 

	gStyle->SetOptStat(0);	
	gStyle->SetEndErrorSize(4);	

	TString PathToFiles = "myXSec/";

	// ---------------------------------------- //

	TH1D::SetDefaultSumw2();
	vector<TString> FileNames; FileNames.clear();
	vector<TString> Label; Label.clear();

	// ---------------------------------------- //

	// Watch out, we did it only for the DeltaPT plot

	vector<TString> PlotNames; 
	PlotNames.push_back("CCQEMuonCosThetaPlot");

	const int NPlots = PlotNames.size();	

	// ---------------------------------------- //

//	FileNames.push_back("Overlay9_Run1"); Label.push_back("Run1"); 
	//FileNames.push_back("Overlay9_Run1_CV"); Label.push_back("Overlay9 Run1 CV");
//	FileNames.push_back("Overlay9_Run3"); Label.push_back("Run3"); 
	//FileNames.push_back("Overlay9_Run3_CV"); Label.push_back("Overlay9 Run3 CV");

	FileNames.push_back("Overlay9_Combined"); Label.push_back("Combined"); 

	const int NFiles = FileNames.size();

	// ---------------------------------------- //

	for (int WhichFile = 0; WhichFile < NFiles; WhichFile++) {

		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot++) {

			// ---------------------------------------- //

			TString CanvasName = "MCC8_Comparison_"+PlotNames[WhichPlot]+"_"+FileNames[WhichFile];
			TCanvas* can = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			can->SetBottomMargin(0.17);
			can->SetLeftMargin(0.17);

			// ---------------------------------------- //

			TLegend* leg = new TLegend(0.05,0.92,0.95,0.98);
			leg->SetNColumns(3);		

			// ---------------------------------------- //

			// CCQE Wiener SVD

			TFile* f = TFile::Open(PathToFiles+UBCodeVersion+"/CCQEWienerSVD_ExtractedXSec_"+FileNames[WhichFile]+"_"+UBCodeVersion+".root","readonly");
			TH1D* Plots = (TH1D*)(f->Get("RecoFullUnc"+PlotNames[WhichPlot])); // Total total unc
			
			Plots->SetLineColor(kBlue+2);
			Plots->SetLineStyle(kSolid);			
			Plots->SetMarkerColor(kBlue+2);
			PrettyPlot(Plots);
			Plots->SetMarkerStyle(20);
			Plots->SetMarkerSize(2.);			
			Plots->SetTitle("");			
			Plots->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);
			Plots->GetXaxis()->SetNdivisions(8);			

			// ---------------------------------------- //

			// STLV Wiener SVD

			TFile* Wf = TFile::Open(PathToFiles+UBCodeVersion+"/WienerSVD_ExtractedXSec_"+FileNames[WhichFile]+"_"+UBCodeVersion+".root","readonly");
			TH1D* Wh = (TH1D*)(Wf->Get("RecoFullUnc"+PlotNames[WhichPlot])); // Plots total unc

			int n = Wh->GetXaxis()->GetNbins();
			double Nuedges[n+1];   
			for (int i = 0; i < n+1; i++) { Nuedges[i] = Wh->GetBinLowEdge(i+1) + 0.2 * Wh->GetBinWidth(i+1); }

			TString Xtitle = Wh->GetXaxis()->GetTitle();
			TString Ytitle = Wh->GetYaxis()->GetTitle();		
			TH1D* Offset = new TH1D("Offset_"+PlotNames[WhichPlot],";"+Xtitle+";"+Ytitle,n,Nuedges);

			for (int WhichBin = 1; WhichBin <= n; WhichBin++ ) {

				Offset->SetBinContent(WhichBin, Wh->GetBinContent(WhichBin));
				Offset->SetBinError(WhichBin, Wh->GetBinError(WhichBin));

			}

			Offset->SetLineColor(kOrange+7);
			Offset->SetLineStyle(kSolid);			
			Offset->SetMarkerColor(kOrange+7);
			Offset->SetMarkerStyle(22);
			Offset->SetMarkerSize(2.);				
			PrettyPlot(Offset);			

			// ---------------------------------------- //

			// MCC8 data points

			TFile* MCC8f = TFile::Open(PathToFiles+UBCodeVersion+"/CCQE_MCC8Data_All.root","readonly");
			TH1D* MCC8h = (TH1D*)(MCC8f->Get("TrueMuonCosThetaPlot")); // Plots total unc	

			int MCC8n = Wh->GetXaxis()->GetNbins();
			double MCC8Nuedges[n+1];   
			for (int i = 0; i < n+1; i++) { MCC8Nuedges[i] = MCC8h->GetBinLowEdge(i+1) - 0.2 * MCC8h->GetBinWidth(i+1); }			

			TString MCC8fXtitle = MCC8h->GetXaxis()->GetTitle();
			TString MCC8fYtitle = MCC8h->GetYaxis()->GetTitle();		
			TH1D* MCC8Offset = new TH1D("MCC8Offset_"+PlotNames[WhichPlot],";"+MCC8fXtitle+";"+MCC8fYtitle,MCC8n,MCC8Nuedges);

			for (int WhichBin = 1; WhichBin <= MCC8n; WhichBin++ ) {

				MCC8Offset->SetBinContent(WhichBin, MCC8h->GetBinContent(WhichBin));
				MCC8Offset->SetBinError(WhichBin, MCC8h->GetBinError(WhichBin));

			}

			MCC8Offset->SetLineColor(kGreen-3);
			MCC8Offset->SetLineStyle(kSolid);			
			MCC8Offset->SetMarkerColor(kGreen-3);
			MCC8Offset->SetMarkerStyle(21);
			MCC8Offset->SetMarkerSize(2.);				
			PrettyPlot(MCC8Offset);									

			// ---------------------------------------- //				

			double min = 2. * MCC8Offset->GetMinimum();
			double max = 1.2 * Plots->GetMaximum();
//			double max = 1.2 * Offset->GetMaximum();				

			Plots->GetYaxis()->SetRangeUser(min,max);
			Plots->GetYaxis()->SetTitleOffset(1.2);

			Plots->Draw("e1 same");
			Offset->Draw("e1 same");
			MCC8Offset->Draw("e1 same");			

			leg->AddEntry(MCC8Offset,"CCQE MCC8","p");
			leg->AddEntry(Plots,"CCQE MCC9","p");
			leg->AddEntry(Offset,"STLV MCC9","p");

			// ---------------------------------------- //

			leg->SetBorderSize(0);
			leg->SetTextSize(TextSize);
			leg->SetTextFont(FontStyle);
			leg->Draw();

			can->SaveAs("myPlots/pdf/"+UBCodeVersion+"/BeamOn9/"+CanvasName+".pdf");
			delete can;

		} // End of the loop over the plots

	} // End of the loop over the runs

	// -----------------------------------------------------------------------------------------------------------------------------


} // End of the program
