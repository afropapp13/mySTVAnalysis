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

		{ "DeltaAlphaTPlot", "DeltaAlphaTPlot" },
		{ "SerialDeltaAlphaT_ProtonCosThetaPlot_0", "DeltaAlphaT_ProtonCosTheta_Minus1_00To0_00Plot" },
		{ "SerialDeltaAlphaT_ProtonCosThetaPlot_1", "DeltaAlphaT_ProtonCosTheta_0_00To0_50Plot" },
		{ "SerialDeltaAlphaT_ProtonCosThetaPlot_2", "DeltaAlphaT_ProtonCosTheta_0_50To0_75Plot" },
		{ "SerialDeltaAlphaT_ProtonCosThetaPlot_3", "DeltaAlphaT_ProtonCosTheta_0_75To1_00Plot" }

};

//----------------------------------------//

void PRD_DeltaAlphaTInCosThetaPSlices_Genie() {

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

	vector<TString> PlotNames; vector<TString> PanelLabels; vector<double> Min; vector<double> Max;
//	PlotNames.push_back("DeltaAlphaTPlot"); PanelLabels.push_back("(a)"); Min.push_back(0.); Max.push_back(44.);
	PlotNames.push_back("SerialDeltaAlphaT_ProtonCosThetaPlot_0"); PanelLabels.push_back("(a)"); Min.push_back(0.); Max.push_back(0.013);
	PlotNames.push_back("SerialDeltaAlphaT_ProtonCosThetaPlot_1");PanelLabels.push_back("(b)"); Min.push_back(0.); Max.push_back(0.059);	
	PlotNames.push_back("SerialDeltaAlphaT_ProtonCosThetaPlot_2"); PanelLabels.push_back("(c)"); Min.push_back(0.); Max.push_back(0.13);
	PlotNames.push_back("SerialDeltaAlphaT_ProtonCosThetaPlot_3");	PanelLabels.push_back("(d)"); Min.push_back(0.); Max.push_back(0.21);

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
			PlotCanvas->SetTopMargin(0.09);
			PlotCanvas->SetLeftMargin(0.19);
			PlotCanvas->SetRightMargin(0.05);				
			PlotCanvas->Draw();

//			TLegend* legData = new TLegend(0.23,0.67,0.67,0.78);
			TLegend* legData = new TLegend(0.23,0.76,0.67,0.88);			
			legData->SetBorderSize(0);
			legData->SetTextSize(TextSize);
			legData->SetTextFont(FontStyle);
			legData->SetNColumns(1);
			legData->SetMargin(0.08);	
			legData->SetFillStyle(0);	

//			TLegend* legUnc = new TLegend(0.23,0.67,0.67,0.78);
			TLegend* legUnc = new TLegend(0.235,0.7,0.63,0.76);			
			legUnc->SetBorderSize(0);
			legUnc->SetTextSize(TextSize);
			legUnc->SetTextFont(FontStyle);
			legUnc->SetNColumns(2);
			legUnc->SetMargin(0.15);	
			legUnc->SetFillStyle(0);					

//			TLegend* leg = new TLegend(0.24,0.79,0.68,0.88);
			TLegend* leg = new TLegend(0.23,0.5,0.68,0.7);
			if (iplot == 1 || iplot == 2) { leg = new TLegend(0.63,0.7,0.83,0.88); }			
			leg->SetBorderSize(0);
			leg->SetTextSize(TextSize);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(1);
			leg->SetMargin(0.075);
			if (iplot == 1 || iplot == 2) { leg->SetMargin(0.15); }
			leg->SetFillStyle(0);									

			//----------------------------------------//

			// Get the shape + stat data plot & plot it

			TH1D* BeamOnShapeStat = (TH1D*)( fXSec->Get("StatShape_" + PlotNames[iplot]) );
			BeamOnShapeStat->GetYaxis()->SetRangeUser(XSecRange[ Mapping[PlotNames[iplot]] ].first, XSecRange[ Mapping[PlotNames[iplot]] ].second);
			BeamOnShapeStat->GetYaxis()->SetNdivisions(6);
			BeamOnShapeStat->GetYaxis()->SetTitleOffset(1.47);									
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

			// Plot norm only

			TH1D* BeamOnNormOnly = (TH1D*)( fXSec->Get("NormOnly_" + PlotNames[iplot]) );
			BeamOnNormOnly->SetFillColor(kGray+1);
			BeamOnNormOnly->SetLineColor(kGray+1);
			BeamOnNormOnly->SetMarkerColor(kGray+1);			
			BeamOnNormOnly->Draw("e2 same");			

			// Plot the stat+shape again 
			// And then the stat only on top of that

			BeamOnShapeStat->GetYaxis()->SetRangeUser(Min[iplot],Max[iplot]);
			BeamOnShapeStat->Draw("e1x0 same");

			TH1D* BeamOnStatOnly = (TH1D*)( fXSec->Get("StatOnly_" + PlotNames[iplot]) );
			BeamOnStatOnly->SetLineWidth(2);			
			BeamOnStatOnly->Draw("e1x0 same");

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

			PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance/Figures/PRD_Genie_"+PlotNames[iplot]+"_"+Runs[irun]+"_"+UBCodeVersion+".pdf");
			//delete PlotCanvas;	

			//----------------------------------------//

		} // End of the loop over the plots

		//----------------------------------------//					

	} // End of the loop over the runs	

} // End of the program 
