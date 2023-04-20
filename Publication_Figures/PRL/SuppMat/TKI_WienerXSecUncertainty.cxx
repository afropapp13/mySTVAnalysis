#include <TFile.h>
#include <TH1D.h>
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

#include "../../../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

//----------------------------------------//

void TKI_WienerXSecUncertainty() {

	//----------------------------------------//

	int DecimalAccuracy = 2;
	int TextFont = 132;
	double TextSize = 0.06;	

	TH1D::SetDefaultSumw2();	

	TString PathToFiles = "../../../myXSec/";

	gStyle->SetPaintTextFormat("4.2f");		

	//----------------------------------------//

	vector<TString> PlotNames;  vector<TString> SaveFig; vector<double> MaxValue;
	PlotNames.push_back("DeltaPTPlot");  SaveFig.push_back("Fig133"); MaxValue.push_back(48.);
	PlotNames.push_back("DeltaAlphaTPlot");  SaveFig.push_back("Fig134"); MaxValue.push_back(17.);
	PlotNames.push_back("DeltaPtxPlot");  SaveFig.push_back("Fig135"); MaxValue.push_back(39.9);		

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	//----------------------------------------//

	vector<TString> Runs;
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	//----------------------------------------//

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		//----------------------------------------//

		vector<TString> UncSource; UncSource.clear();
		vector<int> Colors; Colors.clear();	

		UncSource.push_back("NuWro"); Colors.push_back(kCyan);
		UncSource.push_back("Stat"); Colors.push_back(kRed+1);
		UncSource.push_back("LY"); Colors.push_back(kGreen+2);
		UncSource.push_back("TPC"); Colors.push_back(kOrange+1);
		UncSource.push_back("SCERecomb2"); Colors.push_back(kBlue);	
		UncSource.push_back("XSec"); Colors.push_back(kMagenta);
		UncSource.push_back("G4"); Colors.push_back(kViolet+1);		
		UncSource.push_back("Flux"); Colors.push_back(kYellow+2);
		UncSource.push_back("Dirt"); Colors.push_back(kCyan+2);
		UncSource.push_back("POT"); Colors.push_back(kGreen);
		UncSource.push_back("NTarget"); Colors.push_back(kMagenta-10);
		UncSource.push_back("MCStat"); Colors.push_back(kGray);			

		const int NSources = UncSource.size();																		

		//----------------------------------------//

		// XSec file
		TFile* fXSec = TFile::Open(PathToFiles+UBCodeVersion+"/WienerSVD_ExtractedXSec_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root","readonly");

		// Unfolding uncertainty file
//		TFile* fUnc = TFile::Open(PathToFiles+UBCodeVersion+"/WienerSVD_UnfoldingUnc_Combined_"+UBCodeVersion+".root","readonly");

		//----------------------------------------//

		// Loop over the plots

		for (int iplot = 0; iplot < N1DPlots; iplot ++) {		

			//----------------------------------------//

			// CV XSec plot
			TH1D* CV = (TH1D*)(fXSec->Get("RecoFullUnc"+PlotNames[iplot]));	
			int NBins = CV->GetXaxis()->GetNbins();				

			//----------------------------------------//		

			// Canvas & legend declaration

			TCanvas* PlotCanvas = new TCanvas(PlotNames[iplot]+"_"+Runs[WhichRun],PlotNames[iplot]+"_"+Runs[WhichRun],205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.14);
			//PlotCanvas->SetTopMargin(0.12);
			//PlotCanvas->SetLeftMargin(0.19);			

			TLegend* leg = new TLegend(0.1,0.9,0.9,0.99);
			leg->SetBorderSize(0);
			leg->SetTextSize(TextSize-0.01);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(7);
			leg->SetMargin(0.15);	

			//----------------------------------------//			

			// First plot the total uncertainty		

			TH1D* TotalUnc = (TH1D*)(CV->Clone());

			for (int ibin = 1; ibin <= NBins; ibin++) {

				// CV entry in a given bin
				double CVEntry = CV->GetBinContent(ibin);
				// Uncertainty in a given bin
				double TotalError = CV->GetBinError(ibin);	
				double FracTotalUnc = 100. * TotalError / CVEntry;	
				TotalUnc->SetBinContent(ibin,FracTotalUnc);

			}

			TotalUnc->SetLineColor(kBlack);
			TotalUnc->SetLineWidth(2);

			TotalUnc->GetXaxis()->SetNdivisions(8);
			TotalUnc->GetXaxis()->SetTitleOffset(1.);
			TotalUnc->GetXaxis()->SetTitleSize(TextSize);
			TotalUnc->GetXaxis()->SetLabelSize(TextSize);
			TotalUnc->GetXaxis()->SetTitleFont(TextFont);
			TotalUnc->GetXaxis()->SetLabelFont(TextFont);
			TotalUnc->GetXaxis()->CenterTitle();

			TotalUnc->GetYaxis()->SetNdivisions(8);
			TotalUnc->GetYaxis()->SetTitleOffset(0.8);
			TotalUnc->GetYaxis()->SetTitleSize(TextSize);
			TotalUnc->GetYaxis()->SetLabelSize(TextSize);
			TotalUnc->GetYaxis()->SetTitleFont(TextFont);
			TotalUnc->GetYaxis()->SetLabelFont(TextFont);			
			TotalUnc->GetYaxis()->SetTitle("Uncertainty [%]");
			TotalUnc->GetYaxis()->CenterTitle();			
			TotalUnc->GetYaxis()->SetRangeUser(0.,MaxValue[iplot]);																		

			TotalUnc->Draw("hist text0 same");	

			leg->AddEntry(TotalUnc,"Total","l");

			//----------------------------------------//

			// Loop over the sources of uncertainty

			TH1D* UncPlot[NSources];	

			for (int isource = 0; isource < NSources; isource++) {

				UncPlot[isource] = (TH1D*)(fXSec->Get(UncSource[isource] + "Reco" + PlotNames[iplot]));

				for (int ibin = 1; ibin <= NBins; ibin++) {

					// CV entry in a given bin
					double CVEntry = CV->GetBinContent(ibin);
					// Uncertainty in a given bin
					double UncError = UncPlot[isource]->GetBinError(ibin);	
					double FracUnc = 100. * UncError / CVEntry;	
					UncPlot[isource]->SetBinContent(ibin,FracUnc);

				}

				UncPlot[isource]->SetLineColor(Colors[isource]);
				UncPlot[isource]->SetMarkerColor(Colors[isource]);
				UncPlot[isource]->SetFillStyle(0);								
				UncPlot[isource]->SetLineWidth(2);

				UncPlot[isource]->Draw("hist text0 same");

				leg->AddEntry(UncPlot[isource],UncSource[isource],"l");
				
			}												

			gPad->RedrawAxis();
			leg->Draw();	

			//----------------------------------------//

			// Saving the canvas with the uncertainties on the final xsec
			PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance_PRD_Rename/"+SaveFig[iplot]+".pdf");
//			delete PlotCanvas;					

			//----------------------------------------//

		} // End of the loop over the plots

		//----------------------------------------//					

	} // End of the loop over the runs	

} // End of the program 
