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

#include "../../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

//----------------------------------------//

void PRL_Supp_XSec_DeltaPtxInDeltaPty_IntBreakdown() {

	//----------------------------------------//

	int DecimalAccuracy = 2;
	int FontStyle = 132;
	double TextSize = 0.06;	

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(4);			

	//----------------------------------------//

	vector<TString> PlotNames;
	PlotNames.push_back("DeltaPtxPlot");
	PlotNames.push_back("SerialDeltaPtx_DeltaPtyPlot_0");
	PlotNames.push_back("SerialDeltaPtx_DeltaPtyPlot_1");
//	PlotNames.push_back("SerialDeltaPtx_DeltaPtyPlot_2");			

	const int NPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << NPlots << endl;

	//----------------------------------------//

	// MC Samples to compare

	vector<TString> MCSampleBand; vector<TString> Label;

	MCSampleBand.push_back("OverlayGENIE"); Label.push_back("G18 LFG");	
	MCSampleBand.push_back("Genie_v3_0_6_RFG"); Label.push_back("G18 RFG");
	MCSampleBand.push_back("Genie_v3_0_6_EffSF"); Label.push_back("G18 EffSF");	
	MCSampleBand.push_back("Genie_v3_0_6_NoRPA"); Label.push_back("G18 No RPA");	

//	MCSampleBand.push_back("G21hA"); Label.push_back("G21 hA");	
//	MCSampleBand.push_back("SuSav2"); Label.push_back("G21 hN");	
//	MCSampleBand.push_back("G21G4"); Label.push_back("G21 G4");	
//	MCSampleBand.push_back("G21NoFSI"); Label.push_back("G21 No FSI");
//	MCSampleBand.push_back("GiBUU"); Label.push_back("GiBUU");
//	MCSampleBand.push_back("GiBUUNoFSI"); Label.push_back("GiBUU No FSI");	
//	MCSampleBand.push_back("Genie_v3_0_6_NoFSI"); Label.push_back("G18 No FSI");	
	// MCSampleBand.push_back("GiBUUTscaling"); Label.push_back("GiBUUTscaling");	
	// MCSampleBand.push_back("NEUT");  Label.push_back("NEUT");
	// MCSampleBand.push_back("NEUTv5401_RFG");  Label.push_back("NEUTv5401_RFG");	
	// MCSampleBand.push_back("Overlay9NuWro"); Label.push_back("NuWro");
	// MCSampleBand.push_back("GENIEv2"); Label.push_back("Gv2");
	// MCSampleBand.push_back("GENIEv2LFG"); Label.push_back("Gv2 LFG");
	// MCSampleBand.push_back("GENIEv2EffSF"); Label.push_back("Gv2 EffSF");		
	// MCSampleBand.push_back("Genie_v3_0_6_Out_Of_The_Box"); Label.push_back("G18 No Tune");						
	// MCSampleBand.push_back("Genie_v3_0_6_hN2018"); Label.push_back("G18 hN Tune");			

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

				TLegend* leg = new TLegend(0.25,0.75,0.55,0.85);
				leg->SetBorderSize(0);
				leg->SetTextSize(TextSize);
				leg->SetTextFont(FontStyle);
				leg->SetNColumns(2);
				leg->SetMargin(0.15);

				TLegend* legData = new TLegend(0.25,0.6,0.55,0.75);
				legData->SetBorderSize(0);
				legData->SetTextSize(TextSize);
				legData->SetTextFont(FontStyle);
				legData->SetNColumns(1);
				legData->SetMargin(0.15);							

				//----------------------------------------//

				// Get the shape + stat data plot & plot it

				TH1D* BeamOnShapeStat = (TH1D*)( fXSec->Get("StatShape_" + PlotNames[iplot]) );
				BeamOnShapeStat->Draw("e1x0 same");

				//----------------------------------------//

				// Loop over the MC interaction channel predictions

				TH1D* MCPlot[NInte];

				TString THStackName = MCSampleBand[igen] + "THStack_" + PlotNames[iplot]+"_"+Runs[irun];
				THStack* thstack = new THStack(THStackName,"");				

				// Starting from 1, 0 corresponds to all the events 
				// 1 = QE, 2 = MEC, 3 = RES, 4 = DIS, 5 = COH 
				// Skip COH

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
				BeamOnNormOnly->SetFillColorAlpha(kMagenta-9, 0.45);
				BeamOnNormOnly->SetFillStyle(3008);					
				BeamOnNormOnly->SetLineColor(kMagenta-9);
				BeamOnNormOnly->SetMarkerColor(kMagenta-9);			
				BeamOnNormOnly->Draw("e2 same");						

				//----------------------------------------//	

				legData->AddEntry(BeamOnShapeStat,"MicroBooNE Data","ep");	
				legData->AddEntry(BeamOnShapeStat,"(Stat #oplus Shape Unc)","");
				legData->AddEntry(BeamOnNormOnly,"Norm Unc","f");				

				if (iplot == 0 && igen == 0) { 

					leg->Draw();
					legData->Draw();					 
					
				}

				TLatex *textSlice = new TLatex();
				textSlice->SetTextFont(FontStyle);
				textSlice->SetTextSize(TextSize);
				textSlice->DrawLatexNDC(0.2, 0.92, Label[igen] + "    " +LatexLabel[ MapUncorCor[PlotNames[iplot]] ]);			

				//----------------------------------------//

				PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance/Figures/PRL_Supp_"+MCSampleBand[igen]+"_"+PlotNames[iplot]+"_"+Runs[irun]+"_"+UBCodeVersion+".pdf");
				PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/Papadopoulou_MITThesis/templates/Figures/PRL_Supp_"+MCSampleBand[igen]+"_"+PlotNames[iplot]+"_"+Runs[irun]+"_"+UBCodeVersion+".pdf");	
				//delete PlotCanvas;						

			} // End of the loop over the generators

			//----------------------------------------//												

		} // End of the loop over the plots

		//----------------------------------------//					

	} // End of the loop over the runs	

} // End of the program 
