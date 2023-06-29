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

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../../../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "/home/afroditi/Dropbox/PhD/Secondary_Code/mySimFunctions.cpp"

static std::map<TString,TString> Mapping = {

		{ "DeltaAlphaTPlot", "DeltaAlphaTPlot" },
		{ "SerialDeltaAlphaT_DeltaPTPlot_0", "DeltaAlphaT_DeltaPT_0_00To0_20Plot" },
		{ "SerialDeltaAlphaT_DeltaPTPlot_1", "DeltaAlphaT_DeltaPT_0_20To0_40Plot" },
		{ "SerialDeltaAlphaT_DeltaPTPlot_2", "DeltaAlphaT_DeltaPT_0_40To1_00Plot" },

};

//----------------------------------------//

void PRL_SuppMat_GiBUU_Fig2_DeltaAlphaTInDeltaPTSlices() {

	//----------------------------------------//

	int DecimalAccuracy = 2;
	int FontStyle = 132;
	double TextSize = 0.05;	
	double LegendTextSize = 0.06;		

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(4);			

	//----------------------------------------//

	vector<TString> PlotNames; vector<TString> PanelLabels; vector<double> Min; vector<double> Max; vector<TString> Units; vector<TString> SaveFig;
	PlotNames.push_back("DeltaAlphaTPlot"); PanelLabels.push_back("(a)"); Min.push_back(0.); Max.push_back(0.104); Units.push_back("[$10^{-38}\\frac{cm^{2}}{deg\\,^{40}Ar}$]"); SaveFig.push_back("Fig84");
	PlotNames.push_back("SerialDeltaAlphaT_DeltaPTPlot_0"); PanelLabels.push_back("(b)"); Min.push_back(0.); Max.push_back(0.39); Units.push_back("[$10^{-38}\\frac{cm^{2}}{deg\\,(GeV/\\textit{c})\\,^{40}Ar}$]"); SaveFig.push_back("Fig85");
	PlotNames.push_back("SerialDeltaAlphaT_DeltaPTPlot_2");	PanelLabels.push_back("(c)"); Min.push_back(0.); Max.push_back(0.054); Units.push_back("[$10^{-38}\\frac{cm^{2}}{deg\\,(GeV/\\textit{c})\\,^{40}Ar}$]"); SaveFig.push_back("Fig86");

	const int NPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << NPlots << endl;

	//----------------------------------------//

	// MC Samples to compare

	vector<TString> MCSampleBand; vector<TString> Label; vector<int> MCColors;  vector<int> LineStyle;

	MCSampleBand.push_back("GiBUUNoFSI"); Label.push_back("GiBUU no-FSI"); MCColors.push_back(NuWroColor); LineStyle.push_back(kDashed);
	MCSampleBand.push_back("GiBUU"); Label.push_back("GiBUU FSI"); MCColors.push_back(NuWroColor); LineStyle.push_back(kSolid);	
	MCSampleBand.push_back("Genie_v3_0_6_NoFSI"); Label.push_back("G18 no-FSI"); MCColors.push_back(OverlayColor); LineStyle.push_back(kDashed);		
	MCSampleBand.push_back("OverlayGENIE"); Label.push_back("G18 FSI"); MCColors.push_back(OverlayColor); LineStyle.push_back(kSolid);			

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

		TString XSecFileName = "../../../myXSec/" + UBCodeVersion + "/GenXSec/All_XSecs_Combined_" + UBCodeVersion + ".root";
		TFile* fXSec = new TFile(XSecFileName,"readonly");	

		//----------------------------------------//		

		// Loop over the plots

		for (int iplot = 0; iplot < NPlots; iplot ++) {

			//----------------------------------------//			

			// Canvas & legend

			TString CanvasName = "Fig2_" + PlotNames[iplot]+"_"+Runs[irun];
			TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.14);
			PlotCanvas->SetTopMargin(0.1);
			PlotCanvas->SetLeftMargin(0.18);
			PlotCanvas->SetRightMargin(0.03);				
			PlotCanvas->Draw();

//			TLegend* legData = new TLegend(0.22,0.68,0.77,0.78);
			TLegend* legData = new TLegend(0.27,0.3,0.82,0.42);	
			if (iplot == 1) { legData = new TLegend(0.28,0.7,0.82,0.82); }
			if (iplot == 2) { legData = new TLegend(0.22,0.76,0.77,0.88); }					
			legData->SetBorderSize(0);
			legData->SetTextSize(LegendTextSize);
			legData->SetTextFont(FontStyle);
			legData->SetNColumns(1);
			legData->SetMargin(0.07);
			legData->SetFillStyle(0);

			TLegend* legUnc = new TLegend(0.28,0.24,0.68,0.3);	
			if (iplot == 1) { legUnc = new TLegend(0.28,0.64,0.68,0.7); }
			if (iplot == 2) { legUnc = new TLegend(0.22,0.7,0.63,0.76); }						
			legUnc->SetBorderSize(0);
			legUnc->SetTextSize(LegendTextSize);
			legUnc->SetTextFont(FontStyle);
			legUnc->SetNColumns(2);
			legUnc->SetMargin(0.18);
			legUnc->SetFillStyle(0);						

//			TLegend* leg = new TLegend(0.23,0.78,0.72,0.88);
			TLegend* leg = new TLegend(0.19,0.22,0.96,0.34);
			if (iplot == 1) { leg = new TLegend(0.21,0.72,0.96,0.84); }
			if (iplot == 2) { leg = new TLegend(0.225,0.46,0.8,0.7); }			
			leg->SetBorderSize(0);
			leg->SetTextSize(LegendTextSize);
			leg->SetTextFont(FontStyle);			
			leg->SetNColumns(2);			
			leg->SetMargin(0.13);
			if (iplot != 2) { leg->SetMargin(0.07); }
			leg->SetFillStyle(0);	

			if (iplot == 2) { 
				
				leg->SetNColumns(1); 
				leg->SetMargin(0.06);

			}

			//----------------------------------------//

			// Get the shape + stat data plot & plot it

			TH1D* BeamOnShapeStat = (TH1D*)( fXSec->Get("StatShape_" + PlotNames[iplot]) );
			BeamOnShapeStat->GetYaxis()->SetRangeUser(Min[iplot],Max[iplot]);

			BeamOnShapeStat->Draw("e1x0 same");	

			TH2D* Cov = (TH2D*)fXSec->Get("UnfCov_"+PlotNames[iplot]);	

			//----------------------------------------//

			TH2D* Ac = (TH2D*)fXSec->Get("Ac_"+PlotNames[iplot]);

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

			TH1D* BeamOnNormOnly = (TH1D*)( fXSec->Get("NormOnly_" + PlotNames[iplot]) );
			BeamOnNormOnly->SetFillColorAlpha(NormUncBandColor, 0.45);
			BeamOnNormOnly->SetLineColor(kWhite);
			BeamOnNormOnly->SetMarkerColor(NormUncBandColor);
			BeamOnNormOnly->Draw("e2 same");	

			// Plot the stat+shape again 
			// And then the stat only & norm only on top of that

			BeamOnShapeStat->Draw("e1x0 same");

			TH1D* BeamOnStatOnly = (TH1D*)( fXSec->Get("StatOnly_" + PlotNames[iplot]) );
			BeamOnStatOnly->Draw("e1x0 same");

			//----------------------------------------//

			// Legend & Run / POT

			double tor860_wcut = Fulltor860_wcut_Combined;
			TString Label = ToString(tor860_wcut).ReplaceAll("e"," #times 10").ReplaceAll("+","^{")+"} POT";	

			TLatex *textPOT = new TLatex();
			textPOT->SetTextFont(FontStyle);
			textPOT->SetTextSize(0.06);
			//if (iplot == 0) { textPOT->DrawLatexNDC(0.7, 0.94,Label);	}												

			//----------------------------------------//	

			legData->AddEntry(BeamOnShapeStat,"MicroBooNE Data","");
			legData->AddEntry(BeamOnShapeStat,Label,"");			
//			legData->AddEntry(BeamOnShapeStat,Label,"");

			legUnc->AddEntry(BeamOnShapeStat,"Stat#oplusShape","ep");
			legUnc->AddEntry(BeamOnNormOnly,"Norm","f");				

			if (iplot == 2) { legData->Draw(); }
			if (iplot == 2) { legUnc->Draw();}					
			leg->Draw();

			TLatex *text = new TLatex();
			text->SetTextFont(FontStyle);
			text->SetTextSize(0.06);
			text->DrawLatexNDC(0.2, 0.94,  PanelLabels[iplot] + " " + LatexLabel[ Mapping[PlotNames[iplot]] ]);	

			//----------------------------------------//
/*
			// Transparency pad

			TPad* pad = new TPad("pad","pad",0.3,0.2,0.85,0.38, 21);
			if (iplot == 2) { pad = new TPad("pad","pad",0.2,0.47,0.75,0.65, 21); }

			pad->SetFillColor(kWhite);
			pad->SetFillStyle(0);			 
			PlotCanvas->cd();
			pad->Draw();
			pad->cd();
			pad->SetBottomMargin(0.17);
			pad->SetTopMargin(0.0);
			pad->SetRightMargin(0.08);
			pad->SetLeftMargin(0.1);		

			TH1D* hA = (TH1D*)(MCPlot[1]->Clone());
			hA->Divide(MCPlot[0]);
			hA->SetLineStyle(kSolid);
			//hA->GetXaxis()->SetRangeUser(0.,0.7);			

			hA->GetXaxis()->SetNdivisions(6);
			hA->GetXaxis()->SetTitle("");
			hA->GetXaxis()->SetLabelSize(0.2);
			hA->GetXaxis()->SetLabelFont(FontStyle);

			hA->GetYaxis()->SetNdivisions(5);		
			hA->GetYaxis()->SetLabelSize(0.2);
			hA->GetYaxis()->SetLabelFont(FontStyle);	
			hA->GetYaxis()->SetRangeUser(0.51,1.39);
			if (iplot == 1) { hA->GetYaxis()->SetRangeUser(0.41,0.99); }
			if (iplot == 2) { hA->GetYaxis()->SetRangeUser(0.61,2.79); }																		

			hA->Draw("hist same");

			TH1D* hN = (TH1D*)(MCPlot[2]->Clone());
			hN->Divide(MCPlot[0]);
			hN->SetLineStyle(kSolid);	
			hN->Draw("hist same");

			TH1D* G4 = (TH1D*)(MCPlot[3]->Clone());
			G4->Divide(MCPlot[0]);
			G4->SetLineStyle(kSolid);	
			G4->Draw("hist same");			

			TLatex *textFSI = new TLatex();
			textFSI->SetTextFont(FontStyle);
			textFSI->SetTextSize(0.2);							
			textFSI->DrawLatexNDC(0.2, 0.84, "FSI/No FSI");			

			gPad->RedrawAxis();										
*/
			//----------------------------------------//

			gPad->RedrawAxis();			

			PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance_PRL_Rename/"+SaveFig[iplot]+".pdf");
			//PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_Neutrino2022_PublicNote/Figures/PRL_Fig2_"+PlotNames[iplot]+"_"+Runs[irun]+"_"+UBCodeVersion+".pdf");			
			//PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/Papadopoulou_MITThesis/templates/Figures/PRL_Fig2_"+PlotNames[iplot]+"_"+Runs[irun]+"_"+UBCodeVersion+".pdf");	
			delete PlotCanvas;															

		} // End of the loop over the plots

		//----------------------------------------//					

	} // End of the loop over the runs	

} // End of the program 