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
			PlotCanvas->Draw();

			TLegend* leg = new TLegend(0.6,0.65,0.7,0.85);

			if (
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
					
					leg = new TLegend(0.22,0.58,0.32,0.85); 				
					
			}

			leg->SetBorderSize(0);
			leg->SetTextSize(0.03);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(1);
			leg->SetMargin(0.15);			

			//----------------------------------------//

			// Get the shape + stat data plot & plot it

			TH1D* BeamOnShapeStat = (TH1D*)( fXSec->Get("StatShape_" + PlotNames[iplot]) );
			BeamOnShapeStat->Draw("e1x0 same");

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
			int NBins = BeamOnShapeStat->GetXaxis()->GetNbins();
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

			TH1D* MinValuePlot = (TH1D*)( BeamOnShapeStat->Clone() );
			TH1D* SpreadPlot = (TH1D*)( BeamOnShapeStat->Clone() );			

			GetMinAndSpread(MinValuePlot,SpreadPlot,BinEntries);

			TString THStackName = "THStack_" + PlotNames[iplot]+"_"+Runs[irun];
			THStack* thstack = new THStack(THStackName,"");

			MinValuePlot->SetLineColor(kWhite);
			MinValuePlot->SetFillColor(kWhite);			
			thstack->Add(MinValuePlot,"hist");
			thstack->Draw("same");

			SpreadPlot->SetLineColor(OverlayColor);
			SpreadPlot->SetFillColorAlpha(OverlayColor, 0.45);			
			thstack->Add(SpreadPlot,"hist");
			thstack->Draw("same");			

			//----------------------------------------//	

			// Plot the stat+shape again 
			// And then the stat only & norm only on top of that

			BeamOnShapeStat->Draw("e1x0 same");

			TH1D* BeamOnStatOnly = (TH1D*)( fXSec->Get("StatOnly_" + PlotNames[iplot]) );
			BeamOnStatOnly->Draw("e1x0 same");

			TH1D* BeamOnNormOnly = (TH1D*)( fXSec->Get("NormOnly_" + PlotNames[iplot]) );
			BeamOnNormOnly->SetFillColorAlpha(kGray+1, 0.45);	
			BeamOnNormOnly->SetLineColor(kGray+1);
			BeamOnNormOnly->SetMarkerColor(kGray+1);			
			BeamOnNormOnly->Draw("e2 same");						

			//----------------------------------------//	

			leg->AddEntry(BeamOnShapeStat,"MicroBooNE Data","ep");	
			leg->AddEntry(BeamOnShapeStat,"(Stat #oplus Shape Unc)","");
			leg->AddEntry(BeamOnNormOnly,"Norm Unc","f");
			leg->AddEntry(SpreadPlot,"Generator Spread","f");				

			leg->Draw();

			TLatex *textSlice = new TLatex();
			textSlice->SetTextFont(FontStyle);
			textSlice->SetTextSize(0.06);
			textSlice->DrawLatexNDC(0.24, 0.92, LatexLabel[ MapUncorCor[PlotNames[iplot]] ]);			

			//----------------------------------------//

			PlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/BeamOn9/GeneratorBand_"+PlotNames[iplot]+"_"+Runs[irun]+"_"+UBCodeVersion+".pdf");
			delete PlotCanvas;						

			//----------------------------------------//												

		} // End of the loop over the plots

		//----------------------------------------//					

	} // End of the loop over the runs	

} // End of the program 
