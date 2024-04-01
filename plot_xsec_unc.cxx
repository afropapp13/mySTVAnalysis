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

#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

//----------------------------------------//

void plot_xsec_unc() {

	//----------------------------------------//

	int DecimalAccuracy = 2;
	int TextFont = 132;
	double TextSize = 0.06;	
	gStyle->SetPaintTextFormat("4.2f");		

	TH1D::SetDefaultSumw2();	

	//----------------------------------------//

	//vector<TString> PlotNames;
	//PlotNames.push_back("MuonCosThetaPlot");
	//PlotNames.push_back("MuonCosThetaSingleBinPlot");

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
		TFile* fXSec = TFile::Open(PathToExtractedXSec + "/WienerSVD_ExtractedXSec_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root","readonly");

		//----------------------------------------//

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {		

			//----------------------------------------//

			// CV XSec plot
			TH1D* CV = (TH1D*)(fXSec->Get("RecoFullUnc"+PlotNames[WhichPlot]));	
			int NBins = CV->GetXaxis()->GetNbins();				

			//----------------------------------------//		

			// Canvas & legend declaration

			TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+"_"+Runs[WhichRun],PlotNames[WhichPlot]+"_"+Runs[WhichRun],205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.15);

			TLegend* leg = new TLegend(0.15,0.91,0.85,0.99);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.03);
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

			if ( PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot")  { TotalUnc->GetXaxis()->SetNdivisions(0); }		
			else { TotalUnc->GetXaxis()->SetNdivisions(8); }	
			
			TotalUnc->GetXaxis()->SetTitleSize(TextSize);
			TotalUnc->GetXaxis()->SetLabelSize(TextSize);
			TotalUnc->GetXaxis()->SetTitleFont(TextFont);
			TotalUnc->GetXaxis()->SetLabelFont(TextFont);			
	
			TotalUnc->GetYaxis()->SetNdivisions(8);
			TotalUnc->GetYaxis()->SetTitleOffset(0.8);
			TotalUnc->GetYaxis()->SetTitleSize(TextSize);
			TotalUnc->GetYaxis()->SetLabelSize(TextSize);
			TotalUnc->GetYaxis()->SetTitleFont(TextFont);
			TotalUnc->GetYaxis()->SetLabelFont(TextFont);			
			TotalUnc->GetYaxis()->SetTitle("Uncertainty [%]");
			TotalUnc->GetYaxis()->CenterTitle();			
			if ( PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot")  { TotalUnc->GetYaxis()->SetRangeUser(0.,14.9); }		
			else {  TotalUnc->GetYaxis()->SetRangeUser(0.,49.9); }	

			TotalUnc->Draw("hist text0 same");	
			leg->AddEntry(TotalUnc,"Total","l");
			
			//----------------------------------------//

			// Loop over the sources of uncertainty

			TH1D* UncPlot[NSources];	

			for (int isource = 0; isource < NSources; isource++) {

				UncPlot[isource] = (TH1D*)(fXSec->Get(UncSource[isource] + "Reco" + PlotNames[WhichPlot]));

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

			//----------------------------------------//

			// For the single bin study, print out the corresponding uncertainties

			if (PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot") {

				double NuWroUnc = UncPlot[0]->GetBinContent(1);
				double StatUnc = TMath::Sqrt( TMath::Power(UncPlot[1]->GetBinContent(1),2.) + TMath::Power(UncPlot[11]->GetBinContent(1),2.) );		
				double DetUnc = TMath::Sqrt( TMath::Power(UncPlot[2]->GetBinContent(1),2.) + TMath::Power(UncPlot[3]->GetBinContent(1),2.) + TMath::Power(UncPlot[4]->GetBinContent(1),2.) );				
				double XSecUnc = TMath::Sqrt( TMath::Power(UncPlot[5]->GetBinContent(1),2.) + TMath::Power(NuWroUnc,2.) );
				double FluxUnc = UncPlot[7]->GetBinContent(1);
				double G4Unc = UncPlot[6]->GetBinContent(1);
				double DirtUnc = UncPlot[8]->GetBinContent(1);
				double POTUnc = UncPlot[9]->GetBinContent(1);
				double NTargetUnc = UncPlot[10]->GetBinContent(1);

				cout << " Flux = " <<  FluxUnc << " %" << endl;
				cout << " XSec = " <<  XSecUnc << " %" << endl;
				cout << " Det = " <<  DetUnc << " %" << endl;
				cout << " NuWro = " <<  NuWroUnc << " %" << endl;
				cout << " POT = " <<  POTUnc << " %" << endl;				
				cout << " Stat = " <<  StatUnc << " %" << endl;	
				cout << " NTarget = " <<  NTargetUnc << " %" << endl;						
				cout << " G4 = " <<  G4Unc << " %" << endl;		
				cout << " Dirt = " <<  DirtUnc << " %" << endl;			
	

				double ManualTotalUnc = 0;

				for (int isource = 0; isource < NSources; isource++){

					ManualTotalUnc += TMath::Power(UncPlot[isource]->GetBinContent(1),2.);

				}		

				ManualTotalUnc = TMath::Sqrt(ManualTotalUnc);

				cout << endl << " Manual Total = " << ManualTotalUnc << " %" << endl;	

				TH1D* TotalClone = (TH1D*)(TotalUnc->Clone());	
				TotalClone->SetLineColor(kBlack);
				TotalClone->SetLineWidth(2);
				TotalClone->SetBinContent(1,ManualTotalUnc);																		
			}

			gPad->RedrawAxis();
			leg->Draw();	

			//----------------------------------------//

			// Saving the canvas with the uncertainties on the final xsec
			PlotCanvas->SaveAs(PlotPath + "/Data9/XSecUnc_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
			delete PlotCanvas;					

			//----------------------------------------//

		} // End of the loop over the plots

		//----------------------------------------//					

	} // End of the loop over the runs	

} // End of the program 
