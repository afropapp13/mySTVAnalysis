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
#include <THStack.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "/home/afroditi/Dropbox/PhD/Secondary_Code/mySimFunctions.cpp"

//----------------------------------------//

void PRD_DeltaPTInCosThetaMu_InteBreakdown() {

	//----------------------------------------//

	int DecimalAccuracy = 2;
	int FontStyle = 132;
	double TextSize = 0.06;	

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(4);			

	//----------------------------------------//

	vector<TString> PlotNames; vector<double> Min; vector<double> Max; vector<TString> Label;
	//PlotNames.push_back("DeltaPTPlot");
	PlotNames.push_back("SerialDeltaPT_MuonCosThetaPlot_0"); Min.push_back(0.); Max.push_back(7.99); Label.push_back("(a)");
	PlotNames.push_back("SerialDeltaPT_MuonCosThetaPlot_1"); Min.push_back(0.); Max.push_back(22.); Label.push_back("(b)");
	PlotNames.push_back("SerialDeltaPT_MuonCosThetaPlot_2"); Min.push_back(0.); Max.push_back(44.); Label.push_back("(c)");
	PlotNames.push_back("SerialDeltaPT_MuonCosThetaPlot_3"); Min.push_back(0.); Max.push_back(63.);	Label.push_back("(d)");

	const int NPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << NPlots << endl;

	//----------------------------------------//

	// MC Samples to compare

	vector<TString> MCSampleBand;

	MCSampleBand.push_back("OverlayGENIE");
//	MCSampleBand.push_back("GiBUU");
//	MCSampleBand.push_back("GiBUUNoFSI");	
//	MCSampleBand.push_back("GiBUUTscaling");	
//	MCSampleBand.push_back("NEUT");
//	MCSampleBand.push_back("NEUTv5401_RFG");	
//	MCSampleBand.push_back("Overlay9NuWro");
//	MCSampleBand.push_back("GENIEv2");
	// MCSampleBand.push_back("GENIEv2LFG");
	// MCSampleBand.push_back("GENIEv2EffSF");		
	// MCSampleBand.push_back("Genie_v3_0_6_Out_Of_The_Box");					
	// MCSampleBand.push_back("SuSav2");
	// MCSampleBand.push_back("G21hA");
	// MCSampleBand.push_back("G21G4");
	// MCSampleBand.push_back("G21NoFSI");		
	// MCSampleBand.push_back("Genie_v3_0_6_hN2018");
	// MCSampleBand.push_back("Genie_v3_0_6_NoFSI");
	// MCSampleBand.push_back("Genie_v3_0_6_NoRPA");
	// MCSampleBand.push_back("Genie_v3_0_6_RFG");
	// MCSampleBand.push_back("Genie_v3_0_6_EffSF");			

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

		TString XSecFileName = "../../myXSec/" + UBCodeVersion + "/GenXSec/All_XSecs_Combined_" + UBCodeVersion + ".root";
		TFile* fXSec = new TFile(XSecFileName,"readonly");

		//----------------------------------------//		

		// Loop over the plots

		for (int iplot = 0; iplot < NPlots; iplot ++) {

			//----------------------------------------//

			// Loop over the generators

			for (int igen = 0; igen < NMC; igen++) {

				//----------------------------------------//			

				// Canvas & legend

				TString CanvasName = MCSampleBand[igen] + "_" + PlotNames[iplot]+"_"+Runs[irun];
				TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
				PlotCanvas->cd();
				PlotCanvas->SetBottomMargin(0.14);
				PlotCanvas->SetTopMargin(0.12);
				PlotCanvas->SetLeftMargin(0.19);
				PlotCanvas->SetRightMargin(0.03);				
				PlotCanvas->Draw();

				TLegend* leg = new TLegend(0.62,0.76,0.85,0.85);
				leg->SetBorderSize(0);
				leg->SetTextSize(0.05);
				leg->SetTextFont(FontStyle);
				leg->SetNColumns(2);
				leg->SetMargin(0.13);

				TLegend* legData = new TLegend(0.61,0.59,0.94,0.76);
				legData->SetBorderSize(0);
				legData->SetTextSize(0.05);
				legData->SetTextFont(FontStyle);
				legData->SetNColumns(1);
				legData->SetMargin(0.08);		

				//----------------------------------------//

				// Get the shape + stat data plot & plot it

				TH1D* BeamOnShapeStat = (TH1D*)( fXSec->Get("StatShape_" + PlotNames[iplot]) );
				BeamOnShapeStat->GetYaxis()->SetRangeUser(Min[iplot],Max[iplot]);
				BeamOnShapeStat->Draw("e1x0 same");

				//----------------------------------------//

				// Loop over the MC interaction channel predictions

				TH1D* MCPlot[NInte];

				TString THStackName = MCSampleBand[igen] + "THStack_" + PlotNames[iplot]+"_"+Runs[irun];
				THStack* thstack = new THStack(THStackName,"");				

				// Starting from 1, 0 corresponds to all the events 
				// 1 = QE, 2 = MEC, 3 = RES, 4 = DIS, 5 = COH 
				// Ignore COH

				for (int imc = 1; imc < NInte-1; imc++) {	

					TString InteName = InteractionLabels[imc] + MCSampleBand[igen] + "_" + PlotNames[iplot];
					MCPlot[imc] = (TH1D*)( fXSec->Get(InteName) );
					MCPlot[imc]->SetLineColor(InteBreakColors[imc]);
					MCPlot[imc]->SetFillColor(InteBreakColors[imc]);					
					thstack->Add(MCPlot[imc],"hist");
					thstack->Draw("same");

					TLegendEntry* lGenie = leg->AddEntry(MCPlot[imc],InteractionLabels[imc],"l");
					lGenie->SetTextColor(InteBreakColors[imc]);					

				}			

				//----------------------------------------//	

				// Plot the stat+shape again 
				// And then the stat only & norm only on top of that

				BeamOnShapeStat->Draw("e1x0 same");

				TH1D* BeamOnStatOnly = (TH1D*)( fXSec->Get("StatOnly_" + PlotNames[iplot]) );
				BeamOnStatOnly->Draw("e1x0 same");

				TH1D* BeamOnNormOnly = (TH1D*)( fXSec->Get("NormOnly_" + PlotNames[iplot]) );
				BeamOnNormOnly->SetFillColor(kGray+1);
				BeamOnNormOnly->SetLineColor(kGray+1);
				BeamOnNormOnly->SetMarkerColor(kGray+1);		
				BeamOnNormOnly->Draw("e2 same");						

				//----------------------------------------//	

				legData->AddEntry(BeamOnShapeStat,"MicroBooNE Data","ep");
				legData->AddEntry(BeamOnShapeStat,"(Stat #oplus Shape)","");					
				legData->AddEntry(BeamOnShapeStat,ToString(Fulltor860_wcut_Combined)+" POT","");
				legData->AddEntry(BeamOnNormOnly,"Norm","f");				

				leg->Draw();
				legData->Draw();

				TLatex *textSlice = new TLatex();
				textSlice->SetTextFont(FontStyle);
				textSlice->SetTextSize(TextSize);
				textSlice->DrawLatexNDC(0.2, 0.92, Label[iplot] + " " +LatexLabel[ MapUncorCor[PlotNames[iplot]] ]);			

				gPad->RedrawAxis();

				//----------------------------------------//

				PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance/Figures/PRD_InteBreakDown_" + MCSampleBand[igen] + "_"+PlotNames[iplot]+"_"+Runs[irun]+"_"+UBCodeVersion+".pdf");
				delete PlotCanvas;						

			} // End of the loop over the generators

			//----------------------------------------//												

		} // End of the loop over the plots

		//----------------------------------------//					

	} // End of the loop over the runs	

} // End of the program 
