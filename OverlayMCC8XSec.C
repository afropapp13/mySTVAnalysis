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

void OverlayMCC8XSec() {

	// ----------------------------------------------------------------------------------------------------------------------------------------

	int LineWidth = 3;
	int FontStyle = 132;
	int Ndivisions = 6; 
	double TextSize = 0.06; 

	gStyle->SetOptStat(0);	

	TString PathToFiles = "myXSec/";

	// ----------------------------------------------------------------------------------------------------------------------------------------

	const std::vector<int> Colors{kBlack,kBlue-7,410,610,kRed+1,kGreen+3,kBlue};

	// -----------------------------------------------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	vector<TString> FileNames; FileNames.clear();
	vector<TString> Label; Label.clear();
	
	vector<TH1D*> Plots; Plots.clear();
	vector<TH1D*> PlotsCC1p; PlotsCC1p.clear();
	
	vector<TH1D*> WPlots; WPlots.clear();
	vector<TH1D*> WPlotsCC1p; WPlotsCC1p.clear();	

	// -----------------------------------------------------------------------------------------------------------------------------

//	TString PlotName = "DeltaPTPlot";
//	TString PlotName = "DeltaAlphaTPlot";
//	TString PlotName = "DeltaPhiTPlot";
//	TString PlotName = "MuonMomentumPlot";
	TString PlotName = "MuonCosThetaPlot";
//	TString PlotName = "MuonPhiPlot";
//	TString PlotName = "ProtonMomentumPlot";
//	TString PlotName = "ProtonCosThetaPlot";
//	TString PlotName = "ProtonPhiPlot";

	// -----------------------------------------------------------------------------------------------------------------------------

	FileNames.push_back("Overlay9_Run1"); Label.push_back("Run1"); 
	//FileNames.push_back("Overlay9_Run1_CV"); Label.push_back("Overlay9 Run1 CV");
//	FileNames.push_back("Overlay9_Run3"); Label.push_back("Run3"); 
	//FileNames.push_back("Overlay9_Run3_CV"); Label.push_back("Overlay9 Run3 CV");

	const int NFiles = FileNames.size();

	// -----------------------------------------------------------------------------------------------------------------------------

	TString CanvasName = "OverlayCanvas";
	TCanvas* can = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
	can->SetBottomMargin(0.17);
	can->SetLeftMargin(0.15);

	// -----------------------------------------------------------------------------------------------------------------------------

	TLegend* leg = new TLegend(0.5,0.6,0.7,0.8);
	leg->SetNColumns(1);

	// -----------------------------------------------------------------------------------------------------------------------------

	for (int WhichFile = 0; WhichFile < NFiles; WhichFile++) {

//		TFile* f = TFile::Open(PathToFiles+UBCodeVersion+"/ExtractedXSec_"+FileNames[WhichFile]+"_"+UBCodeVersion+".root","readonly");

//		TH1D* h = (TH1D*)(f->Get("TotalReco"+PlotName));
//		Plots.push_back(h);
//		Plots[WhichFile]->SetLineColor(Colors[WhichFile]);
//		Plots[WhichFile]->SetMarkerColor(Colors[WhichFile]);
//		PrettyPlot(Plots[WhichFile]);
//		Plots[WhichFile]->SetMarkerStyle(20);
//		Plots[WhichFile]->GetYaxis()->SetTitle(VarLabel[PlotName]);
//		Plots[WhichFile]->Draw("e same");
//		leg->AddEntry(Plots[WhichFile],"EE BeamOn " + Label[WhichFile],"lep");

//		TH1D* hCC1p = (TH1D*)(f->Get("CC1pReco"+PlotName));
//		PlotsCC1p.push_back(hCC1p);
//		PlotsCC1p[WhichFile]->SetLineColor(Colors[WhichFile]);
//		PlotsCC1p[WhichFile]->SetMarkerColor(Colors[WhichFile]);
//		PlotsCC1p[WhichFile]->SetFillColor(Colors[WhichFile]);
//		PrettyPlot(PlotsCC1p[WhichFile]);
//		//PlotsCC1p[WhichFile]->SetFillColorAlpha(Colors[WhichFile],0.35);
//		PlotsCC1p[WhichFile]->SetFillStyle(3001);
//		
//		PlotsCC1p[WhichFile]->Draw("e2 same");
//		leg->AddEntry(PlotsCC1p[WhichFile],"MC " + Label[WhichFile],"f");

		// Wiener SVD

		TFile* Wf = TFile::Open(PathToFiles+UBCodeVersion+"/WienerSVD_ExtractedXSec_"+FileNames[WhichFile]+"_"+UBCodeVersion+".root","readonly");
		
		TH1D* Wh = (TH1D*)(Wf->Get("Reco"+PlotName));
		WPlots.push_back(Wh);
		WPlots[WhichFile]->SetLineColor(Colors[WhichFile]);
		WPlots[WhichFile]->SetMarkerColor(Colors[WhichFile]);
		WPlots[WhichFile]->SetMarkerStyle(20);		
		PrettyPlot(WPlots[WhichFile]);
		WPlots[WhichFile]->Draw("e same");
		leg->AddEntry(WPlots[WhichFile],"SVD BeamOn " + Label[WhichFile],"lep");

		TH1D* WhCC1p = (TH1D*)(Wf->Get("True"+PlotName));
		WPlotsCC1p.push_back(WhCC1p);
		WPlotsCC1p[WhichFile]->SetLineColor(Colors[WhichFile+1]);
		WPlotsCC1p[WhichFile]->SetMarkerColor(Colors[WhichFile+1]);
//		WPlotsCC1p[WhichFile]->SetFillColor(Colors[WhichFile]);
		PrettyPlot(WPlotsCC1p[WhichFile]);
//		WPlotsCC1p[WhichFile]->SetFillColorAlpha(Colors[WhichFile],0.35);
//		WPlotsCC1p[WhichFile]->SetFillStyle(3001);
		
//		WPlotsCC1p[WhichFile]->Draw("e2 same");
		WPlotsCC1p[WhichFile]->Draw("hist same");
		leg->AddEntry(WPlotsCC1p[WhichFile],"MC " + Label[WhichFile],"l");			

	}

	// -----------------------------------------------------------------------------------------------------------------------------

	leg->SetBorderSize(0);
	leg->SetTextSize(TextSize);
	leg->SetTextFont(FontStyle);
	leg->Draw();

	// -----------------------------------------------------------------------------------------------------------------------------


} // End of the program
