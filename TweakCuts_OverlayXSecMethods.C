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

// -------------------------------------------------------------------------------------------------------------------------------------

int LocateBinWithValue(TH1D* h, double Value) {

	int NBins = h->GetXaxis()->GetNbins();

	for (int i = 1; i <= NBins; i++) {

		double CurrentEntry = h->GetBinContent(i);
		if (CurrentEntry == Value) { return i; } 

	}

	return -99;

}

// -------------------------------------------------------------------------------------------------------------------------------------

void PrettyPlot(TH1D* h,int LineWidth = 2, int FontStyle = 132, int Ndivisions = 6, double TextSize = 0.06) {

	// ----------------------------------------------------------------------------------------------------------------

	h->SetLineWidth(LineWidth);

	// ----------------------------------------------------------------------------------------------------------------

	// X-axis

	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetLabelSize(TextSize);
	h->GetXaxis()->SetTitleSize(TextSize);
	h->GetXaxis()->SetTitleOffset(1.05);
	h->GetXaxis()->SetNdivisions(Ndivisions);

	// ----------------------------------------------------------------------------------------------------------------

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

void TweakCuts_OverlayXSecMethods() {

	// ----------------------------------------------------------------------------------------------------------------------------------------

	int LineWidth = 3;
	int FontStyle = 132;
	int Ndivisions = 6; 
	double TextSize = 0.06; 

	gStyle->SetOptStat(0);	
	gStyle->SetEndErrorSize(4);	

	TString PathToFiles = "myXSec/";

	// ----------------------------------------------------------------------------------------------------------------------------------------

	const std::vector<int> Colors{kBlack,610,410,kRed+1,kGreen+3,kBlue};

	// -----------------------------------------------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	vector<TString> FileNames; FileNames.clear();
	vector<TString> Label; Label.clear();

	// -----------------------------------------------------------------------------------------------------------------------------

	// Watch out, we did it only for the DeltaPT plot

	vector<TString> PlotNames; 
	PlotNames.push_back("DeltaPTPlot");
	PlotNames.push_back("DeltaAlphaTPlot");
	PlotNames.push_back("DeltaPhiTPlot");		

	const int NPlots = PlotNames.size();	

	// -----------------------------------------------------------------------------------------------------------------------------

//	FileNames.push_back("Overlay9_Run1"); Label.push_back("Run1"); 
	//FileNames.push_back("Overlay9_Run1_CV"); Label.push_back("Overlay9 Run1 CV");
//	FileNames.push_back("Overlay9_Run3"); Label.push_back("Run3"); 
	//FileNames.push_back("Overlay9_Run3_CV"); Label.push_back("Overlay9 Run3 CV");

	FileNames.push_back("Overlay9_Combined"); Label.push_back("Combined"); 

	const int NFiles = FileNames.size();

	// -----------------------------------------------------------------------------------------------------------------------------

	for (int WhichFile = 0; WhichFile < NFiles; WhichFile++) {

		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot++) {

			// -----------------------------------------------------------------------------------------------------------------------------

			TString CanvasName = "TweakCuts_Comparison_"+PlotNames[WhichPlot]+"_"+FileNames[WhichFile];
			TCanvas* can = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			can->SetBottomMargin(0.17);
			can->SetLeftMargin(0.17);

			// -----------------------------------------------------------------------------------------------------------------------------

			TLegend* leg = new TLegend(0.15,0.92,0.95,0.98);
			leg->SetNColumns(2);		

			// -----------------------------------------------------------------------------------------------------------------------------		

			// Old proton pid score cut at 0

			TFile* f = TFile::Open(PathToFiles+UBCodeVersion+"/WienerSVD_ExtractedXSec_"+FileNames[WhichFile]+"_"+UBCodeVersion+"_ProtonPIDCutAtZero.root","readonly");

			TH1D* Plots = (TH1D*)(f->Get("RecoFullUnc"+PlotNames[WhichPlot]));
			Plots->SetLineColor(kBlack);
			Plots->SetLineStyle(kSolid);			
			Plots->SetMarkerColor(kBlack);
			PrettyPlot(Plots);
			Plots->SetMarkerStyle(20);
			Plots->SetTitle("");			
			Plots->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);

			// Old proton pid score cut at 0.05

			TFile* Wf = TFile::Open(PathToFiles+UBCodeVersion+"/WienerSVD_ExtractedXSec_"+FileNames[WhichFile]+"_"+UBCodeVersion+".root","readonly");
			
			TH1D* Wh = (TH1D*)(Wf->Get("RecoFullUnc"+PlotNames[WhichPlot])); // Plots with total uncertainty

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

			Offset->SetLineColor(kRed);
			Offset->SetLineStyle(kSolid);			
			Offset->SetMarkerColor(kRed);
			Offset->SetMarkerStyle(22);
			Offset->SetMarkerSize(2.);				
			PrettyPlot(Offset);

			double MaxValue = TMath::Max(Plots->GetMaximum(),Offset->GetMaximum());
			int MaxValueBin = LocateBinWithValue(Plots,MaxValue);
			double MaxValueError = TMath::Max(Plots->GetBinError(MaxValueBin),Offset->GetBinError(MaxValueBin));

			double MinValue = TMath::Min(Plots->GetMinimum(),Offset->GetMinimum());
			int MinValueBin = LocateBinWithValue(Offset,MinValue);
			double MinValueError = TMath::Max(Plots->GetBinError(MinValueBin),Offset->GetBinError(MinValueBin));		

			double min = TMath::Min(0., 0.8*(MinValue-MinValueError));
			double max = TMath::Max(0., 1.22*(MaxValue+MaxValueError));	

			Plots->GetYaxis()->SetRangeUser(min,max);
			Plots->GetYaxis()->SetTitleOffset(1.2);

			Plots->Draw("e1 same");
			Offset->Draw("e1 same");

			leg->AddEntry(Plots,"Score < 0","p");
			leg->AddEntry(Offset,"Score < 0.05","p");

			// -----------------------------------------------------------------------------------------------------------------------------

			leg->SetBorderSize(0);
			leg->SetTextSize(TextSize);
			leg->SetTextFont(FontStyle);
			leg->Draw();

			//can->SaveAs("myPlots/pdf/"+UBCodeVersion+"/BeamOn9/"+CanvasName+".pdf");
			//delete can;

		} // End of the loop over the plots

	} // End of the loop over the runs

	// -----------------------------------------------------------------------------------------------------------------------------


} // End of the program
