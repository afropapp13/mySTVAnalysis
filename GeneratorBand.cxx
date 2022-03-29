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

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

//----------------------------------------//

void PrettyPlot(TH1D* h) {

	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetTitleSize(0.06);
	h->GetXaxis()->SetLabelSize(0.06);
	h->GetXaxis()->SetTitleOffset(1.05);
	h->GetXaxis()->SetNdivisions(8);
	h->GetXaxis()->CenterTitle();

	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetNdivisions(8);
	h->GetYaxis()->SetTitleOffset(1.35);
	h->GetYaxis()->SetTitleSize(0.06);
	h->GetYaxis()->SetLabelSize(0.06);
	h->GetYaxis()->CenterTitle();

}

//----------------------------------------//

void GetMinAndSpread(TH1D* MinValuePlot, TH1D* SpreadPlot, std::vector< std::vector<double> > BinEntries){

	int NBins = MinValuePlot->GetXaxis()->GetNbins();

	for (int ibin = 0; ibin < NBins; ibin++) {

		// For a given bin, get the min / max / spread values

		double min = *std::min_element(BinEntries[ibin].begin(), BinEntries[ibin].end());
		double max = *std::max_element(BinEntries[ibin].begin(), BinEntries[ibin].end());
		double spread = TMath::Abs(max-min);

		MinValuePlot->SetBinContent(ibin+1,min);
		SpreadPlot->SetBinContent(ibin+1,spread);		

	}

	return;

}

//----------------------------------------//

void GeneratorBand() {

	//----------------------------------------//

	int DecimalAccuracy = 2;
	int FontStyle = 132;
	double TextSize = 0.06;	

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(4);			

	//----------------------------------------//

	vector<TString> PlotNames;
	PlotNames = MultiDimXSec;

//	PlotNames.push_back("DeltaPTPlot"); 


	const int NPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << NPlots << endl;

	//----------------------------------------//

	// MC Samples to compare

	vector<TString> MCSampleBand;

	MCSampleBand.push_back("OverlayGENIE");
	MCSampleBand.push_back("GiBUU");
	MCSampleBand.push_back("NEUT");
	MCSampleBand.push_back("Overlay9NuWro");	

	int NMC = MCSampleBand.size();

	//----------------------------------------//

	vector<TString> Runs;
//	Runs.push_back("Run1");
//	Runs.push_back("Run2");	
//	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	//----------------------------------------//

	for (int irun = 0; irun < NRuns; irun++) {

		//----------------------------------------//

		// Open the file that contains all the xsecs

		TString XSecFileName = "myXSec/" + UBCodeVersion + "/GenXSec/All_XSecs_Combined_" + UBCodeVersion + ".root";
		TFile* fXSec = new TFile(XSecFileName,"readonly");

		//----------------------------------------//		

		// Loop over the plots

		for (int iplot = 0; iplot < NPlots; iplot ++) {

			//----------------------------------------//			

			// Canvas & legend

			TString CanvasName = PlotNames[iplot]+"_"+Runs[irun];
			TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.14);
			PlotCanvas->SetTopMargin(0.12);
			PlotCanvas->SetLeftMargin(0.19);
			PlotCanvas->SetRightMargin(0.03);			
			PlotCanvas->Draw();

			TLegend* leg = new TLegend(0.6,0.7,0.7,0.85);

			if (
				PlotNames[iplot] == "MuonCosThetaPlot" ||
				PlotNames[iplot] == "ProtonCosThetaPlot" ||				
				PlotNames[iplot] == "DeltaPtyPlot" ||
				PlotNames[iplot] == "SerialDeltaPty_DeltaPtxPlot_0" ||
				PlotNames[iplot] == "SerialDeltaPty_DeltaPtxPlot_1" ||
				PlotNames[iplot] == "SerialDeltaPty_DeltaPtxPlot_2" ||								
				PlotNames[iplot] == "SerialProtonCosTheta_MuonCosThetaPlot_0" ||
				PlotNames[iplot] == "SerialProtonCosTheta_MuonCosThetaPlot_1" ||
				PlotNames[iplot] == "SerialProtonCosTheta_MuonCosThetaPlot_2" ||
				PlotNames[iplot] == "SerialProtonCosTheta_MuonCosThetaPlot_3" ||												
				PlotNames[iplot] == "SerialDeltaPn_DeltaPTPlot_2" ||
				PlotNames[iplot] == "SerialDeltaAlphaT_DeltaPTPlot_1" ||
				PlotNames[iplot] == "SerialMuonMomentum_MuonCosThetaPlot_3" ||
				PlotNames[iplot] == "SerialProtonMomentum_ProtonCosThetaPlot_3"								
			) { 
					
					leg = new TLegend(0.22,0.7,0.32,0.85); 				
					
			}

			leg->SetBorderSize(0);
			leg->SetTextSize(0.03);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(1);
			leg->SetMargin(0.15);			

			//----------------------------------------//

			// Get the data only & plot

//			TH1D* BeamOnXSec = (TH1D*)( fXSec->Get("XSecOnly_" + PlotNames[iplot]) ); // Only xsec uncertainties
			TH1D* BeamOnXSec = (TH1D*)( fXSec->Get("FullUnc_" + PlotNames[iplot]) ); // Only xsec uncertainties			
	
			PrettyPlot(BeamOnXSec);
			BeamOnXSec->SetLineColor(BeamOnColor);
			BeamOnXSec->SetMarkerColor(BeamOnColor);
			BeamOnXSec->SetLineWidth(1);		
			BeamOnXSec->SetMarkerSize(1.);
			BeamOnXSec->SetMarkerStyle(20);	
			BeamOnXSec->GetYaxis()->SetTitle(VarLabel[PlotNames[iplot]]);		
			BeamOnXSec->GetYaxis()->SetRangeUser(XSecRange[ PlotNames[iplot] ].first,XSecRange[ PlotNames[iplot] ].second);																

			BeamOnXSec->Draw("e1x0 same");

			//----------------------------------------//

			// Loop over the MC gen predictions

			TH1D* MCPlot[NMC];

			for (int imc = 0; imc < NMC; imc++) {	

				MCPlot[imc] = (TH1D*)( fXSec->Get(MCSampleBand[imc] + "_" + PlotNames[iplot]) );

			}			

			//----------------------------------------//
		
			// Store the bin entries for each one of the MC samples

			// Index 1 = Bin number, Index 2 = mc sample
			std::vector< std::vector<double> > BinEntries;
			int NBins = BeamOnXSec->GetXaxis()->GetNbins();
			BinEntries.resize(NBins);

			for (int ibin = 0; ibin < NBins; ibin++) {

				BinEntries[ibin].resize(NMC);

				for (int imc = 0; imc < NMC; imc++) {

					BinEntries[ibin][imc] = MCPlot[imc]->GetBinContent(ibin+1);

				}

			}

			//----------------------------------------//
		
			// Store the min values & the spread = max - min 
			// On a bin-by-bin basis

			TH1D* MinValuePlot = (TH1D*)( BeamOnXSec->Clone() );
			TH1D* SpreadPlot = (TH1D*)( BeamOnXSec->Clone() );			

			GetMinAndSpread(MinValuePlot,SpreadPlot,BinEntries);

			TString THStackName = "THStack_" + PlotNames[iplot]+"_"+Runs[irun];
			THStack* thstack = new THStack(THStackName,"");

			MinValuePlot->SetLineColor(kWhite);
			MinValuePlot->SetFillColor(kWhite);			
			thstack->Add(MinValuePlot,"hist");
			thstack->Draw("same");

			SpreadPlot->SetLineColor(kWhite);
			SpreadPlot->SetFillColorAlpha(OverlayColor, 0.45);			
			thstack->Add(SpreadPlot,"hist");
			thstack->Draw("same");			

			//----------------------------------------//	

			// Plot the xsec beam on plot again 

			BeamOnXSec->Draw("e1x0 same");					

			//----------------------------------------//	

			leg->AddEntry(BeamOnXSec,"MicroBooNE Data","ep");	
			leg->AddEntry(BeamOnXSec,"(Total Unc)","");
			leg->AddEntry(SpreadPlot,"Generator Spread","f");				

			leg->Draw();

			TLatex *textSlice = new TLatex();
			textSlice->SetTextFont(FontStyle);
			textSlice->SetTextSize(0.06);
			textSlice->DrawLatexNDC(0.24, 0.92, LatexLabel[ MapUncorCor[PlotNames[iplot]] ]);	

			gPad->RedrawAxis();					

			//----------------------------------------//

			PlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/BeamOn9/GeneratorBand_"+PlotNames[iplot]+"_"+Runs[irun]+"_"+UBCodeVersion+".pdf");
			delete PlotCanvas;						

			//----------------------------------------//												

		} // End of the loop over the plots

		//----------------------------------------//					

	} // End of the loop over the runs	

} // End of the program 
