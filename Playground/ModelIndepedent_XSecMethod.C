#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

#include "../../myClasses/Constants.h"

using namespace Constants;
using namespace std;

//----------------------------------------//

double GetSpread(std::vector<double> BinContent) {

	double min = *std::min_element(BinContent.begin(), BinContent.end());
	double max = *std::max_element(BinContent.begin(), BinContent.end());

	// After discussion with S. Gardiner, dividing by sqrt(2) 
	// following the recipe outlined in BaBar's not
	// https://babar.heprc.uvic.ca/BFROOT/www/Statistics/Report/report.pdf

//	double spread = TMath::Abs(max-min) / TMath::Sqrt(12);
	double spread = TMath::Abs(max-min) / TMath::Sqrt(2);

	return spread;

}

//----------------------------------------//

int LocateBinWithValue(TH1D* h, double Value) {

	int NBins = h->GetXaxis()->GetNbins();

	for (int i = 1; i <= NBins; i++) {

		double CurrentEntry = h->GetBinContent(i);
		if (CurrentEntry == Value) { return i; } 

	}

	return -99;

}

//----------------------------------------//

void PrettyPlot(TH1D* h,int LineWidth = 2, int FontStyle = 132, int Ndivisions = 6, double TextSize = 0.06) {

	//----------------------------------------//

	h->SetLineWidth(LineWidth);

	//----------------------------------------//

	// X-axis

	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetLabelSize(TextSize);
	h->GetXaxis()->SetTitleSize(TextSize);
	h->GetXaxis()->SetTitleOffset(1.05);
	h->GetXaxis()->SetNdivisions(Ndivisions);

	//----------------------------------------//

	// Y-axis

	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleSize(TextSize); 
	h->GetYaxis()->SetTickSize(0.02);
	h->GetYaxis()->SetLabelSize(TextSize);
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleOffset(1.05);
	h->GetYaxis()->SetNdivisions(Ndivisions);

//	h->Scale(1./h->GetMaximum());

	h->GetYaxis()->SetRangeUser(0.,1.15*h->GetMaximum());	

	return;	

}

//----------------------------------------//

void ModelIndepedent_XSecMethod() {

	//----------------------------------------//

	TH1D::SetDefaultSumw2();

	int LineWidth = 3;
	int FontStyle = 132;
	int Ndivisions = 6; 
	double TextSize = 0.06; 

	gStyle->SetOptStat(0);	
	gStyle->SetEndErrorSize(4);	

	TString PathToFiles = "../myXSec/";

	//----------------------------------------//

	const std::vector<int> Colors{kBlue+1,kOrange+7,kGreen+1,610,410,kRed+1,kBlue};

	//----------------------------------------//

	// Watch out, we did it only for the DeltaPT plot

//	vector<TString> PlotNames; 
//	PlotNames.push_back("DeltaPTPlot");

	const int NPlots = PlotNames.size();	

	//----------------------------------------//

	vector<TString> MC; MC.clear();
	vector<TString> FileNames; FileNames.clear();
	vector<TString> Label; Label.clear();	

	//----------------------------------------//

	// Alternative MC Studies
	// Beam On as always

	FileNames.push_back(""); Label.push_back("Nominal");  MC.push_back("Overlay9");
	FileNames.push_back("AltMCNoTune"); Label.push_back("No Tune");  MC.push_back("Overlay9");
	FileNames.push_back("AltMCTwiceMEC"); Label.push_back("2xMEC");  MC.push_back("Overlay9");	

	// Fake Data Studies

	//FileNames.push_back("Overlay9NuWro"); Label.push_back("Nominal");  MC.push_back("Overlay9");
	//FileNames.push_back("Overlay9NuWro"); Label.push_back("No Tune");  MC.push_back("NoTuneOverlay9");
	//FileNames.push_back("Overlay9NuWro"); Label.push_back("2xMEC");  MC.push_back("TwiceMECOverlay9");				

	//----------------------------------------//	 

	const int NFiles = FileNames.size();

	TFile* f[NFiles];

	//----------------------------------------//

	// 1st index = WhichPlot, 2nd index = WhichFile
	std::vector<std::vector<TH1D*> > Plots; Plots.resize(NPlots);
	std::vector<std::vector<TH1D*> > TruePlots; TruePlots.resize(NPlots);

	// 1st index = WhichPlot, 2nd index = WhichBin, 3rd index = WhichFile 
	std::vector<std::vector<std::vector<double> > > BinEntries; BinEntries.resize(NPlots);		

	//----------------------------------------//

	// File to store the unfolding uncertainties

	TFile* fUnc = TFile::Open(PathToFiles+UBCodeVersion+"/WienerSVD_UnfoldingUnc_Combined_"+UBCodeVersion+".root","recreate");

	//----------------------------------------//	

	for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot++) {

		//----------------------------------------//

		Plots[WhichPlot].resize(NFiles);
		TruePlots[WhichPlot].resize(NFiles);		

		//----------------------------------------//		

		TString CanvasName = "ModelIndependent_Comparison_"+PlotNames[WhichPlot];
		TCanvas* can = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
		can->SetBottomMargin(0.17);
		can->SetLeftMargin(0.17);

		TString TrueCanvasName = "MC_ModelIndependent_Comparison_"+PlotNames[WhichPlot];
		TCanvas* Truecan = new TCanvas(TrueCanvasName,TrueCanvasName,205,34,1024,768);
		Truecan->SetBottomMargin(0.17);
		Truecan->SetLeftMargin(0.17);		

		//----------------------------------------//

		TLegend* leg = new TLegend(0.15,0.92,0.95,0.98);
		leg->SetNColumns(3);		

		//----------------------------------------//

		int n = -1.; // will be the number of bins for a given plot
		TH1D* Uncertainty;		

		//----------------------------------------//

		for (int WhichFile = 0; WhichFile < NFiles; WhichFile++) {	

			f[WhichFile] = TFile::Open(PathToFiles+UBCodeVersion+"/"+FileNames[WhichFile]+"WienerSVD_ExtractedXSec_"+MC[WhichFile]+"_Combined_"+UBCodeVersion+".root","readonly");				
			Plots[WhichPlot][WhichFile] = (TH1D*)(f[WhichFile]->Get("MCStatReco"+PlotNames[WhichPlot])); // Only MC Stat unc
			TruePlots[WhichPlot][WhichFile] = (TH1D*)(f[WhichFile]->Get("True"+PlotNames[WhichPlot])); // Ac * True	

			//------------------------------------//

			n = Plots[WhichPlot][WhichFile]->GetXaxis()->GetNbins();
			double Nuedges[n+1];   
			for (int i = 0; i < n+1; i++) { 
				
				// Introduce offset for plotting purposes
				Nuedges[i] = Plots[WhichPlot][WhichFile]->GetBinLowEdge(i+1) + 0.1 * WhichFile * Plots[WhichPlot][WhichFile]->GetBinWidth(i+1); 

			}			

			//------------------------------------//

			if (WhichFile == 0) { BinEntries[WhichPlot].resize(n);}

			//------------------------------------//

			TString Xtitle = Plots[WhichPlot][WhichFile]->GetXaxis()->GetTitle();
			TString Ytitle = Plots[WhichPlot][WhichFile]->GetYaxis()->GetTitle();		
			TH1D* Offset = new TH1D("Offset_"+PlotNames[WhichPlot],";"+Xtitle+";"+Ytitle,n,Nuedges);
			if (WhichFile == 0) { 
				
				Uncertainty = new TH1D("Uncertainty_"+PlotNames[WhichPlot],";"+Xtitle+";"+Ytitle,n,Nuedges); 
				
			}		

			for (int WhichBin = 1; WhichBin <= n; WhichBin++ ) {

				double BinContent = Plots[WhichPlot][WhichFile]->GetBinContent(WhichBin);
				double BinError = Plots[WhichPlot][WhichFile]->GetBinError(WhichBin);				 

				Offset->SetBinContent(WhichBin,BinContent);
				Offset->SetBinError(WhichBin,BinError);

				BinEntries[WhichPlot][WhichBin-1].push_back(BinContent);

			}			

			Plots[WhichPlot][WhichFile] = (TH1D*)(Offset->Clone());

			//------------------------------------//

			PrettyPlot(Plots[WhichPlot][WhichFile]);
			PrettyPlot(TruePlots[WhichPlot][WhichFile]);			

			Plots[WhichPlot][WhichFile]->SetLineColor(Colors[WhichFile]);
			Plots[WhichPlot][WhichFile]->SetMarkerColor(Colors[WhichFile]);
			Plots[WhichPlot][WhichFile]->SetMarkerSize(2.);
			Plots[WhichPlot][WhichFile]->SetMarkerStyle(20);	
			Plots[WhichPlot][WhichFile]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);			

			TruePlots[WhichPlot][WhichFile]->SetLineColor(Colors[WhichFile]);
			TruePlots[WhichPlot][WhichFile]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);	

			double max = TMath::Max(Plots[WhichPlot][WhichFile]->GetMaximum(),Plots[WhichPlot][0]->GetMaximum());
			double truemax = TMath::Max(TruePlots[WhichPlot][WhichFile]->GetMaximum(),TruePlots[WhichPlot][0]->GetMaximum());
			if (truemax > max) { max = truemax; }

			can->cd();
//			Plots[WhichPlot][WhichFile]->Draw("e1x0 same");
//			if (WhichFile == 0) { Plots[WhichPlot][WhichFile]->Draw("e1x0 same"); }
//			else { Plots[WhichPlot][WhichFile]->Draw("p0 hist same"); }
			Plots[WhichPlot][WhichFile]->Draw("p0 hist same");
			//TruePlots[WhichPlot][WhichFile]->Draw("hist same");			
			leg->AddEntry(Plots[WhichPlot][WhichFile],Label[WhichFile],"p");

			Plots[WhichPlot][0]->GetYaxis()->SetRangeUser(0.,max);
			Plots[WhichPlot][0]->Draw("p0 hist same");			

			Truecan->cd();
			TruePlots[WhichPlot][WhichFile]->Draw("hist same");	

			TruePlots[WhichPlot][0]->GetYaxis()->SetRangeUser(0.,max);
			TruePlots[WhichPlot][0]->Draw("hist same");					

			//------------------------------------//

		} // End of the loop over the files

		leg->SetBorderSize(0);
		leg->SetTextSize(TextSize);
		leg->SetTextFont(FontStyle);
		can->cd();
		leg->Draw();
		Truecan->cd();
		leg->Draw();	

		//------------------------------------//	

		// Assigning the uncertainty as the spread / sqrt(12)

		for (int WhichBin = 0; WhichBin < n; WhichBin++) {

			double spread = GetSpread(BinEntries[WhichPlot][WhichBin]);
			Uncertainty->SetBinContent(WhichBin+1,spread);

		}

		Uncertainty->SetFillColorAlpha(kRed+1,1.);	
		Uncertainty->SetLineColor(kRed+1);
		Uncertainty->SetMarkerColor(kRed+1);
		can->cd();		
		Uncertainty->Draw("e2 hist same");	

		fUnc->cd();		
		Uncertainty->Write("UnfUnc_"+PlotNames[WhichPlot]);

		//------------------------------------//				

		can->SaveAs("../myPlots/pdf/"+UBCodeVersion+"/BeamOn9/"+CanvasName+".pdf");
		Truecan->SaveAs("../myPlots/pdf/"+UBCodeVersion+"/BeamOn9/"+TrueCanvasName+".pdf");		
		delete can;

	} // End of the loop over the plots

	fUnc->Close();

	// -----------------------------------------------------------------------------------------------------------------------------


} // End of the program
