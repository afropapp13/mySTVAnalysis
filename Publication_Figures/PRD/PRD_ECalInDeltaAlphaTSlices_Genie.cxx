#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TVectorD.h>
#include <TMatrixD.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../../../myClasses/Constants.h"
#include "../../../myClasses/Util.h"

using namespace std;
using namespace Constants;

#include "/home/afroditi/Dropbox/PhD/Secondary_Code/mySimFunctions.cpp"

static std::map<TString,TString> Mapping = {

		{ "ECalPlot", "ECalPlot" },
		{ "SerialECal_DeltaAlphaTPlot_0", "ECal_DeltaAlphaT_0_00To45_00Plot" },
		{ "SerialECal_DeltaAlphaTPlot_1", "ECal_DeltaAlphaT_45_00To90_00Plot" },
		{ "SerialECal_DeltaAlphaTPlot_2", "ECal_DeltaAlphaT_90_00To135_00Plot" },
		{ "SerialECal_DeltaAlphaTPlot_3", "ECal_DeltaAlphaT_135_00To180_00Plot" }

};

//----------------------------------------//

void PRD_ECalInDeltaAlphaTSlices_Genie() {

	//----------------------------------------//

	int DecimalAccuracy = 2;
	int FontStyle = 132;
	double TextSize = 0.06;	

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(4);		

	gStyle->SetPalette(55); 
	const Int_t NCont = 999; 
	gStyle->SetNumberContours(NCont); 
	gStyle->SetTitleSize(TextSize,"t"); 
	gStyle->SetTitleFont(FontStyle,"t");
	gStyle->SetOptStat(0);	
	gStyle->SetPaintTextFormat("4.2f");		

	//----------------------------------------//

	vector<TString> PlotNames; vector<TString> PanelLabels; vector<double> Min; vector<double> Max; vector<TString> SaveFig;
//	PlotNames.push_back("ECalPlot"); PanelLabels.push_back("(a)"); Min.push_back(0.); Max.push_back(24.);
	PlotNames.push_back("SerialECal_DeltaAlphaTPlot_0"); PanelLabels.push_back("(a)"); Min.push_back(0.); Max.push_back(0.099); SaveFig.push_back("Fig112");
	PlotNames.push_back("SerialECal_DeltaAlphaTPlot_1");PanelLabels.push_back("(b)"); Min.push_back(0.); Max.push_back(0.11); SaveFig.push_back("Fig113");	
	PlotNames.push_back("SerialECal_DeltaAlphaTPlot_2"); PanelLabels.push_back("(c)"); Min.push_back(0.); Max.push_back(0.159); SaveFig.push_back("Fig114");
	PlotNames.push_back("SerialECal_DeltaAlphaTPlot_3");	PanelLabels.push_back("(d)"); Min.push_back(0.); Max.push_back(0.18); SaveFig.push_back("Fig115");

	const int NPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << NPlots << endl;

	//----------------------------------------//

	// MC Samples to compare

	vector<TString> MCSampleBand; vector<TString> Label; vector<int> MCColors;  vector<int> LineStyle;

//	MCSampleBand.push_back("GiBUUNoFSI"); Label.push_back("GiB No FSI"); MCColors.push_back(NuWroColor); LineStyle.push_back(kDashed);
//	MCSampleBand.push_back("Genie_v3_0_6_NoFSI"); Label.push_back("G18 No FSI"); MCColors.push_back(OverlayColor); LineStyle.push_back(kDashed);	
	MCSampleBand.push_back("OverlayGENIE"); Label.push_back("G18"); MCColors.push_back(OverlayColor); LineStyle.push_back(G18LineStyle);
//	MCSampleBand.push_back("GiBUU"); Label.push_back("GiB     "); MCColors.push_back(GiBUUColor); LineStyle.push_back(kSolid);	
//	MCSampleBand.push_back("GiBUUTscaling"); Label.push_back("GiBUUTscaling");	
//	MCSampleBand.push_back("NEUT");  Label.push_back("NEUT"); MCColors.push_back(NEUTColor); LineStyle.push_back(kSolid);
//	MCSampleBand.push_back("Overlay9NuWro");  Label.push_back("NuWro"); MCColors.push_back(NuWroColor); LineStyle.push_back(kSolid);	
//	MCSampleBand.push_back("NEUTv5401_RFG");  Label.push_back("NEUTv5401_RFG");	
//	MCSampleBand.push_back("Overlay9NuWro"); Label.push_back("NuWro");
//	MCSampleBand.push_back("GENIEv2LFG"); Label.push_back("Gv2 LFG");
//	MCSampleBand.push_back("GENIEv2EffSF"); Label.push_back("Gv2 EffSF");		
	MCSampleBand.push_back("Genie_v3_0_6_Out_Of_The_Box"); Label.push_back("Untuned"); MCColors.push_back(kMagenta); LineStyle.push_back(UntunedLineStyle);						
	MCSampleBand.push_back("SuSav2"); Label.push_back("G21"); MCColors.push_back(kOrange+6); LineStyle.push_back(G21LineStyle);
	MCSampleBand.push_back("GENIEv2"); Label.push_back("Gv2"); MCColors.push_back(kBlue); LineStyle.push_back(Gv2LineStyle);	
//	MCSampleBand.push_back("G21hA"); Label.push_back("G21hA");	
//	MCSampleBand.push_back("G21G4"); Label.push_back("G21G4");
//	MCSampleBand.push_back("G21NoFSI"); Label.push_back("G21NoFSI");		
//	MCSampleBand.push_back("Genie_v3_0_6_hN2018"); Label.push_back("G18 hN Tune");
//	MCSampleBand.push_back("Genie_v3_0_6_NoRPA"); Label.push_back("G18 No RPA Tune");
//	MCSampleBand.push_back("Genie_v3_0_6_RFG"); Label.push_back("G18 RFG Tune");
//	MCSampleBand.push_back("Genie_v3_0_6_EffSF"); Label.push_back("G18 EffSF Tune");				

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

			// Canvas & legend

			TString CanvasName = "PRD_" + PlotNames[iplot]+"_"+Runs[irun];
			TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.14);
			PlotCanvas->SetTopMargin(0.08);
			PlotCanvas->SetLeftMargin(0.19);
			PlotCanvas->SetRightMargin(0.05);				
			PlotCanvas->Draw();

//			TLegend* legData = new TLegend(0.51,0.68,0.97,0.8);
			TLegend* legData = new TLegend(0.24,0.76,0.66,0.88);			
			legData->SetBorderSize(0);
			legData->SetTextSize(TextSize);
			legData->SetTextFont(FontStyle);
			legData->SetNColumns(1);
			legData->SetMargin(0.08);
			legData->SetFillStyle(0);

			TLegend* legUnc = new TLegend(0.232,0.7,0.642,0.76);			
			legUnc->SetBorderSize(0);
			legUnc->SetTextSize(TextSize);
			legUnc->SetTextFont(FontStyle);
			legUnc->SetNColumns(2);
			legUnc->SetMargin(0.2);
			legUnc->SetFillStyle(0);							

//			TLegend* leg = new TLegend(0.52,0.81,0.98,0.9);
			TLegend* leg = new TLegend(0.64,0.7,0.9,0.88);			
			leg->SetBorderSize(0);
			leg->SetTextSize(TextSize);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(1);
			leg->SetMargin(0.08);	
			leg->SetFillStyle(0);										

			//----------------------------------------//

			// Get the shape + stat data plot & plot it

			TH1D* BeamOnShapeStat = (TH1D*)( fXSec->Get("StatShape_" + PlotNames[iplot]) );
			BeamOnShapeStat->GetYaxis()->SetRangeUser(XSecRange[ Mapping[PlotNames[iplot]] ].first, XSecRange[ Mapping[PlotNames[iplot]] ].second);			
			BeamOnShapeStat->GetYaxis()->SetLabelOffset(0.006);
			BeamOnShapeStat->GetXaxis()->SetLabelOffset(0.015);	
			BeamOnShapeStat->GetYaxis()->SetNdivisions(6);			
			BeamOnShapeStat->SetLineWidth(2);
			BeamOnShapeStat->Draw("e1x0 same");	

			TH2D* Cov = (TH2D*)fXSec->Get("UnfCov_"+PlotNames[iplot]);				

			//----------------------------------------//						

			TH1D* MCPlot[NMC];
			double Chi2[NMC];
			int Ndof[NMC];
			double pval[NMC];

			// Loop over the MC predictions

			for (int igen = 0; igen < NMC; igen++) {			

				TString InteName = MCSampleBand[igen] + "_" + PlotNames[iplot];
				MCPlot[igen] = (TH1D*)( fXSec->Get(InteName) );
				MCPlot[igen]->SetLineColor(MCColors[igen]);		
				MCPlot[igen]->SetLineStyle(LineStyle[igen]);
				MCPlot[igen]->SetLineWidth(3);					
				MCPlot[igen]->Draw("same hist");	

				CalcChiSquared(MCPlot[igen],BeamOnShapeStat,Cov,Chi2[igen],Ndof[igen],pval[igen]);
				TString Chi2NdofAlt = " (" + to_string_with_precision(Chi2[igen],1) + "/" + TString(std::to_string(Ndof[igen])) +")";									

				TLegendEntry* lGenie = leg->AddEntry(MCPlot[igen],Label[igen] + Chi2NdofAlt,"l");
				lGenie->SetTextColor(MCColors[igen]);

			} // End of the loop over the generators

			//----------------------------------------//	

			// Plot the stat+shape again 
			// And then the stat only & norm only on top of that

			BeamOnShapeStat->GetYaxis()->SetRangeUser(Min[iplot],Max[iplot]);
			BeamOnShapeStat->Draw("e1x0 same");

			TH1D* BeamOnStatOnly = (TH1D*)( fXSec->Get("StatOnly_" + PlotNames[iplot]) );
			BeamOnStatOnly->SetLineWidth(2);
			BeamOnStatOnly->Draw("e1x0 same");

			TH1D* BeamOnNormOnly = (TH1D*)( fXSec->Get("NormOnly_" + PlotNames[iplot]) );
			//BeamOnNormOnly->SetFillColorAlpha(kGray+1, 0.45);
			BeamOnNormOnly->SetFillColor(kGray+1);
			BeamOnNormOnly->SetLineColor(kGray+1);
			BeamOnNormOnly->SetMarkerColor(kGray+1);			
			BeamOnNormOnly->Draw("e2 same");

			//----------------------------------------//

			// Legend & Run / POT

			double tor860_wcut = Fulltor860_wcut_Combined;
			TString Label = ToString(tor860_wcut).ReplaceAll("e"," #times 10").ReplaceAll("+","^{")+"} POT";	

			TLatex *textPOT = new TLatex();
			textPOT->SetTextFont(FontStyle);
			textPOT->SetTextSize(0.06);
			//if (iplot == 0) { textPOT->DrawLatexNDC(0.7, 0.94,Label);	}							

			gPad->RedrawAxis();									

			//----------------------------------------//	

			legData->AddEntry(BeamOnShapeStat,"MicroBooNE Data","");
			legData->AddEntry(BeamOnShapeStat,Label,"");			

			legUnc->AddEntry(BeamOnShapeStat,"Stat#oplusShape","ep");
			legUnc->AddEntry(BeamOnNormOnly,"Norm","f");				

			legData->Draw();
			legUnc->Draw();	
			leg->Draw();			

			TLatex *text = new TLatex();
			text->SetTextFont(FontStyle);
			text->SetTextSize(0.06);
			text->DrawLatexNDC(0.2, 0.94, PanelLabels[iplot] + " " + LatexLabel[ Mapping[PlotNames[iplot]] ]);	

			//----------------------------------------//

			PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance_PRD_Rename/"+SaveFig[iplot]+".pdf");
			//PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance/Figures/PRD_Genie_"+PlotNames[iplot]+"_"+Runs[irun]+"_"+UBCodeVersion+".pdf");
			//delete PlotCanvas;	

			//----------------------------------------//

		} // End of the loop over the plots

		//----------------------------------------//					

	} // End of the loop over the runs	

} // End of the program 
