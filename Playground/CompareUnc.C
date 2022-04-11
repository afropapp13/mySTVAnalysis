#include <TCanvas.h>
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

void PrettyPlot(TH1D* h, int LineWidth = 2, int FontStyle = 132, int Ndivisions = 6, double TextSize = 0.06) {

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
	h->GetYaxis()->SetTitleOffset(0.95);
	h->GetYaxis()->SetNdivisions(Ndivisions);

//	h->Scale(1./h->GetMaximum());

	h->GetYaxis()->SetRangeUser(0.,1.15*h->GetMaximum());	

	return;	

}

//----------------------------------------//

double LocateHistoMaxValue(TH1D* h) {

	int NBins = h->GetXaxis()->GetNbins();
	double max = -1.;

	for (int i = 1; i <= NBins; i++) {

		double CurrentEntry = h->GetBinContent(i);
		if (CurrentEntry > max) { max = CurrentEntry; } 

	}

	return max;

}

//----------------------------------------//

void CompareUnc() {

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

	std::vector<int> Colors;

	//----------------------------------------//

//	vector<TString> PlotNames; 
//	PlotNames.push_back("DeltaPTPlot");
//	PlotNames.push_back("DeltaAlphaTPlot");

	const int NPlots = PlotNames.size();	

	//----------------------------------------//

	vector<TString> FileNames; FileNames.clear();
	vector<TString> Label; Label.clear();	
	vector<TString> PlotLabel; PlotLabel.clear();	

	//----------------------------------------//

	// Alternative MC Studies
	// Beam On as always

	FileNames.push_back("NuWro"); Label.push_back("NuWro Data Spread/sqrt(12)"); PlotLabel.push_back("NuWroDataUnfUnc_"); Colors.push_back(kOrange+7);
	FileNames.push_back("NuWro"); Label.push_back("NuWro (Data-MC)/sqrt(12)"); PlotLabel.push_back("NuWroDataMCUnfUnc_"); Colors.push_back(kGreen+1);
	FileNames.push_back(""); Label.push_back("Beam On Data Spread/sqrt(12)"); PlotLabel.push_back("UnfUnc_"); Colors.push_back(kBlue+1);

	//----------------------------------------//	 

	const int NFiles = FileNames.size();

	TFile* f[NFiles];

	//----------------------------------------//

	// 1st index = WhichPlot
	// 2nd index = WhichFile (0 = NuWro Data, 1 = NuWro Data-MC difference, 2 = BeamOn Data)
	std::vector<std::vector<TH1D*> > Plots; Plots.resize(NPlots);
	std::vector<std::vector<TH1D*> > TruePlots; TruePlots.resize(NPlots);	

	//----------------------------------------//	

	for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot++) {

		//----------------------------------------//

		Plots[WhichPlot].resize(NFiles);
		TruePlots[WhichPlot].resize(NFiles);				

		//----------------------------------------//		

		TString CanvasName = "CompareUnc_Comparison_"+PlotNames[WhichPlot];
		TCanvas* can = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
		can->SetTopMargin(0.18);
		can->SetBottomMargin(0.17);
		can->SetLeftMargin(0.17);	

		//----------------------------------------//

		TLegend* leg = new TLegend(0.15,0.83,0.65,0.98);
		leg->SetNColumns(1);	

		//----------------------------------------//

		double max = -1.;

		for (int WhichFile = 0; WhichFile < NFiles; WhichFile++) {	

			f[WhichFile] = TFile::Open(PathToFiles+UBCodeVersion+"/"+FileNames[WhichFile]+"WienerSVD_UnfoldingUnc_Combined_"+UBCodeVersion+".root","readonly");				
			Plots[WhichPlot][WhichFile] = (TH1D*)(f[WhichFile]->Get(PlotLabel[WhichFile]+PlotNames[WhichPlot]));
			TruePlots[WhichPlot][WhichFile] = (TH1D*)(f[WhichFile]->Get("AltTrue"+PlotNames[WhichPlot])); // Ac * True

			//------------------------------------//

			// Division by NuWro Truth & Percentage
			Plots[WhichPlot][WhichFile]->Divide(TruePlots[WhichPlot][0]);
			Plots[WhichPlot][WhichFile]->Scale(100.);

			//------------------------------------//

			PrettyPlot(Plots[WhichPlot][WhichFile]);
			Plots[WhichPlot][WhichFile]->SetLineColor(Colors[WhichFile]);
			Plots[WhichPlot][WhichFile]->SetMarkerColor(Colors[WhichFile]);
			Plots[WhichPlot][WhichFile]->SetMarkerSize(2.);
			Plots[WhichPlot][WhichFile]->SetMarkerStyle(20);	
//			Plots[WhichPlot][WhichFile]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);
			Plots[WhichPlot][WhichFile]->GetYaxis()->SetTitle("Fractional Uncertainty [%]");

			max = TMath::Max(LocateHistoMaxValue(Plots[WhichPlot][WhichFile]),max);

			can->cd();
			Plots[WhichPlot][WhichFile]->Draw("p0 hist same");
			leg->AddEntry(Plots[WhichPlot][WhichFile],Label[WhichFile],"p");

			Plots[WhichPlot][0]->GetYaxis()->SetRangeUser(0.,1.05*max);
			Plots[WhichPlot][0]->Draw("p0 hist same");							

			//------------------------------------//

		} // End of the loop over the files

		leg->SetBorderSize(0);
		leg->SetTextSize(TextSize);
		leg->SetTextFont(FontStyle);
		can->cd();
		leg->Draw();	

		//------------------------------------//				

		can->SaveAs("../myPlots/pdf/"+UBCodeVersion+"/BeamOn9/"+CanvasName+".pdf");
		delete can;

	} // End of the loop over the plots

	// -----------------------------------------------------------------------------------------------------------------------------


} // End of the program