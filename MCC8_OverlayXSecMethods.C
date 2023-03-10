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
	h->GetYaxis()->SetTitle("#frac{d#sigma}{cos#theta_{#mu}} [10^{-38} #frac{cm^{2}}{Ar}]");
	h->GetYaxis()->SetRangeUser(-0.5,21.);
	h->GetYaxis()->SetTitleOffset(1.2);	

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

			TLegend* leg = new TLegend(0.15,0.9,0.95,0.98);
			leg->SetNColumns(2);		

			// ---------------------------------------- //

			// CCQE (MCC9 with MCC8 signal def)

			// Wiener SVD
			//TFile* f = TFile::Open(PathToFiles+UBCodeVersion+"/CCQEWienerSVD_ExtractedXSec_"+FileNames[WhichFile]+"_"+UBCodeVersion+".root","readonly");
			//TH1D* Plots = (TH1D*)(f->Get("RecoFullUnc"+PlotNames[WhichPlot])); // Total total unc

			// Effective efficiency
			TFile* f = TFile::Open(PathToFiles+UBCodeVersion+"/CCQEExtractedXSec_"+FileNames[WhichFile]+"_"+UBCodeVersion+".root","readonly");			
			TH1D* Plots = (TH1D*)(f->Get("TotalReco"+PlotNames[WhichPlot])); // Total total unc
			
			Plots->SetLineColor(kBlue+2);
			Plots->SetLineStyle(kSolid);			
			Plots->SetMarkerColor(kBlue+2);
			PrettyPlot(Plots);
			Plots->SetMarkerStyle(20);
			Plots->SetMarkerSize(2.);			
			Plots->SetTitle("");			
			Plots->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);
			Plots->GetXaxis()->SetNdivisions(8);	

			
			Plots->Draw("e1 same");		
			leg->AddEntry(Plots,"\"PRL\" MCC9","p");	

			TH2D* Cov = (TH2D*)f->Get("UnfCovCCQEMuonCosThetaPlot");								
			

			// ---------------------------------------- //

			// STLV Wiener SVD (MCC9 with CCQE signal defintion)

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

			
			Offset->Draw("e1 same");
			leg->AddEntry(Offset,"Arxiv \"MCC8\"","p");								
			

			// ---------------------------------------- //

			// STLV Wiener SVD (MCC9 w/o CCQE signal defintion)

			TH1D* WhNoCCQE = (TH1D*)(Wf->Get("RecoFullUncMuonCosThetaPlot")); // Plots total unc

			int nNoCCQE = WhNoCCQE->GetXaxis()->GetNbins();
			double NuedgesNoCCQE[nNoCCQE+1];   
			for (int i = 0; i < nNoCCQE+1; i++) { NuedgesNoCCQE[i] = WhNoCCQE->GetBinLowEdge(i+1) + 0.2 * WhNoCCQE->GetBinWidth(i+1); }
		
			TH1D* OffsetNoCCQE = new TH1D("OffsetNoCCQE_"+PlotNames[WhichPlot],";"+Xtitle+";"+Ytitle,nNoCCQE,NuedgesNoCCQE);

			for (int WhichBin = 1; WhichBin <= nNoCCQE; WhichBin++ ) {

				OffsetNoCCQE->SetBinContent(WhichBin, WhNoCCQE->GetBinContent(WhichBin));
				OffsetNoCCQE->SetBinError(WhichBin, WhNoCCQE->GetBinError(WhichBin));

			}

			OffsetNoCCQE->SetLineColor(kMagenta-2);
			OffsetNoCCQE->SetLineStyle(kSolid);			
			OffsetNoCCQE->SetMarkerColor(kMagenta-2);
			OffsetNoCCQE->SetMarkerStyle(23);
			OffsetNoCCQE->SetMarkerSize(2.);				
			PrettyPlot(OffsetNoCCQE);


			OffsetNoCCQE->Draw("e1 same");	
			leg->AddEntry(OffsetNoCCQE,"Arxiv MCC9","p");		

			// ---------------------------------------- //

			// Effective efficiency MCC9

			TFile* WfEE = TFile::Open(PathToFiles+UBCodeVersion+"/ExtractedXSec_"+FileNames[WhichFile]+"_"+UBCodeVersion+".root","readonly");
			TH1D* WhEE = (TH1D*)(WfEE->Get("TotalRecoMuonCosThetaPlot")); // Plots total unc
		
			TH1D* OffsetEE = new TH1D("OffsetEE_"+PlotNames[WhichPlot],";"+Xtitle+";"+Ytitle,nNoCCQE,NuedgesNoCCQE);

			for (int WhichBin = 1; WhichBin <= nNoCCQE; WhichBin++ ) {

				OffsetEE->SetBinContent(WhichBin, WhEE->GetBinContent(WhichBin));
				OffsetEE->SetBinError(WhichBin, WhEE->GetBinError(WhichBin));

			}

			OffsetEE->SetLineColor(kRed+1);
			OffsetEE->SetLineStyle(kSolid);			
			OffsetEE->SetMarkerColor(kRed+1);
			OffsetEE->SetMarkerStyle(24);
			OffsetEE->SetMarkerSize(2.);				
			PrettyPlot(OffsetEE);


			OffsetEE->Draw("e1 same");	
			leg->AddEntry(OffsetEE,"Arxiv MCC9 EE","p");											

			// ---------------------------------------- //

			// MCC8 published data points (PRL 2020)

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


			MCC8Offset->Draw("e1 same");	
			leg->AddEntry(MCC8Offset,"PRL MCC8","p");															

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
