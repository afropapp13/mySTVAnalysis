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
#include <TF1.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../../../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "/home/afroditi/Dropbox/PhD/Secondary_Code/mySimFunctions.cpp"

static std::map<TString,TString> Mapping = {

		{ "DeltaPtxPlot", "DeltaPtxPlot" },
		{ "SerialDeltaPtx_DeltaPtyPlot_0", "DeltaPtx_DeltaPty_Minus0_75ToMinus0_15Plot" },
		{ "SerialDeltaPtx_DeltaPtyPlot_1", "DeltaPtx_DeltaPty_Minus0_15To0_15Plot" },
		{ "SerialDeltaPtx_DeltaPtyPlot_2", "DeltaPtx_DeltaPty_0_15To0_45Plot" },

};

//----------------------------------------//

void PRL_SuppMat_Fig3_DeltaPtxInDeltaPtySlices() {

	//----------------------------------------//

	int DecimalAccuracy = 2;
	int FontStyle = 132;
	double TextSize = 0.05;
	double LegendTextSize = 0.05;		

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(4);			

	//----------------------------------------//

	vector<TString> PlotNames; vector<TString> PanelLabels; vector<double> Min; vector<double> Max; vector<TString> Units; vector<TString> SaveFig;
	PlotNames.push_back("DeltaPtxPlot"); PanelLabels.push_back("(a)"); Min.push_back(0.); Max.push_back(42.); Units.push_back("[$10^{-38}\\frac{cm^{2}}{(GeV/\\textit{c})\\,^{40}Ar}$]"); SaveFig.push_back("Fig66");
	PlotNames.push_back("SerialDeltaPtx_DeltaPtyPlot_0"); PanelLabels.push_back("(b)"); Min.push_back(0.); Max.push_back(14.9); Units.push_back("[$10^{-38}\\frac{cm^{2}}{(GeV/\\textit{c})^2\\,^{40}Ar}$]"); SaveFig.push_back("Fig67");
	PlotNames.push_back("SerialDeltaPtx_DeltaPtyPlot_1"); PanelLabels.push_back("(c)"); Min.push_back(0.); Max.push_back(109.); Units.push_back("[$10^{-38}\\frac{cm^{2}}{(GeV/\\textit{c})^2\\,^{40}Ar}$]"); SaveFig.push_back("Fig68");

	const int NPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << NPlots << endl;

	//----------------------------------------//

	// MC Samples to compare

	vector<TString> MCSampleBand; vector<TString> Label; vector<int> MCColors;  vector<int> LineStyle;

//	MCSampleBand.push_back("OverlayGENIE"); Label.push_back("G18 hA"); MCColors.push_back(OverlayColor); LineStyle.push_back(kSolid);	
//	MCSampleBand.push_back("Genie_v3_0_6_NoFSI"); Label.push_back("G18 no-FSI");  MCColors.push_back(OverlayColor); LineStyle.push_back(kDashed);
//	MCSampleBand.push_back("Genie_v3_0_6_NoRPA"); Label.push_back("no-RPA");  MCColors.push_back(GiBUUColor); LineStyle.push_back(kSolid);
//	MCSampleBand.push_back("Genie_v3_0_6_RFG"); Label.push_back("RFG");  MCColors.push_back(kMagenta-9); LineStyle.push_back(kSolid);
//	MCSampleBand.push_back("Genie_v3_0_6_EffSF"); Label.push_back("EffSF"); MCColors.push_back(NEUTColor); LineStyle.push_back(kSolid);

//	MCSampleBand.push_back("GiBUUNoFSI"); Label.push_back("GiB No FSI"); MCColors.push_back(NuWroColor); LineStyle.push_back(kDashed);
//	MCSampleBand.push_back("GiBUU"); Label.push_back("GiB FSI"); MCColors.push_back(NuWroColor); LineStyle.push_back(kSolid);	
//	MCSampleBand.push_back("GiBUUTscaling"); Label.push_back("GiBUUTscaling");	
//	MCSampleBand.push_back("NEUT");  Label.push_back("NEUT");
//	MCSampleBand.push_back("NEUTv5401_RFG");  Label.push_back("NEUTv5401_RFG");	
//	MCSampleBand.push_back("Overlay9NuWro"); Label.push_back("NuWro");
//	MCSampleBand.push_back("GENIEv2"); Label.push_back("Gv2");
//	MCSampleBand.push_back("GENIEv2LFG"); Label.push_back("Gv2 LFG");
//	MCSampleBand.push_back("GENIEv2EffSF"); Label.push_back("Gv2 EffSF");		
//	MCSampleBand.push_back("Genie_v3_0_6_Out_Of_The_Box"); Label.push_back("G18 No Tune");							
//	MCSampleBand.push_back("Genie_v3_0_6_hN2018"); Label.push_back("G18 hN");  MCColors.push_back(GiBUUColor); LineStyle.push_back(kSolid);	

	// MCSampleBand.push_back("G21NoFSI"); Label.push_back("G21 no-FSI");  MCColors.push_back(kBlue+1); LineStyle.push_back(kDashed);
	// MCSampleBand.push_back("G21hA"); Label.push_back("G21 hA");  MCColors.push_back(kViolet-4); LineStyle.push_back(kSolid);
	// MCSampleBand.push_back("SuSav2"); Label.push_back("G21 hN"); MCColors.push_back(kGreen-2); LineStyle.push_back(kSolid);
	// MCSampleBand.push_back("G21G4"); Label.push_back("G21 G4");  MCColors.push_back(kOrange+6); LineStyle.push_back(kSolid);

	MCSampleBand.push_back("Genie_v3_0_6_NoFSI"); Label.push_back("G18 no-FSI"); MCColors.push_back(OverlayColor); LineStyle.push_back(kDashed);		
	MCSampleBand.push_back("OverlayGENIE"); Label.push_back("G18 hA"); MCColors.push_back(OverlayColor); LineStyle.push_back(kSolid);
	MCSampleBand.push_back("Genie_v3_0_6_hN2018"); Label.push_back("G18 hN"); MCColors.push_back(kOrange+6); LineStyle.push_back(kSolid);
	MCSampleBand.push_back("v3_2_0_G18_10d_02_11a"); Label.push_back("G18 G4"); MCColors.push_back(kGreen-2); LineStyle.push_back(kSolid);	

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

			TString CanvasName = "Fig3_" + PlotNames[iplot]+"_"+Runs[irun];
			TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.14);
			PlotCanvas->SetTopMargin(0.1);
			PlotCanvas->SetLeftMargin(0.185);
			PlotCanvas->SetRightMargin(0.005);				
			PlotCanvas->Draw();			

			TLegend* legData = new TLegend(0.22,0.64,0.77,0.88);
			legData->SetBorderSize(0);
			legData->SetTextSize(LegendTextSize);
			legData->SetTextFont(FontStyle);
			legData->SetNColumns(1);
			legData->SetMargin(0.05);
			legData->SetFillStyle(0);				

			TLegend* leg = new TLegend(0.645,0.64,0.945,0.88);
			leg->SetBorderSize(0);
			leg->SetTextSize(LegendTextSize);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(1);
			leg->SetMargin(0.12);
			leg->SetFillStyle(0);		

			TLegend* legSigma = new TLegend(0.45,0.185,0.61,0.285);
			if (iplot == 2) { legSigma = new TLegend(0.45,0.17,0.61,0.27); }
			legSigma->SetBorderSize(0);
			legSigma->SetTextSize(TextSize);
			legSigma->SetTextFont(FontStyle);
			legSigma->SetNColumns(1);
			legSigma->SetMargin(0.08);
			legSigma->SetFillStyle(0);															

			//----------------------------------------//

			// Get the shape + stat data plot & plot it

			TH1D* BeamOnShapeStat = (TH1D*)( fXSec->Get("StatShape_" + PlotNames[iplot]) );
			BeamOnShapeStat->GetYaxis()->SetRangeUser(Min[iplot],Max[iplot]);			
			BeamOnShapeStat->Draw("e1x0 same");	

			TH2D* Cov = (TH2D*)fXSec->Get("UnfCov_"+PlotNames[iplot]);			

			//----------------------------------------//

			TH2D* Ac = (TH2D*)fXSec->Get("Ac_"+PlotNames[iplot]);					

			//----------------------------------------//				
			//----------------------------------------//						

			TH1D* MCPlot[NMC];
			double Chi2[NMC];
			int Ndof[NMC];
			double pval[NMC];

			// + 1 for data
			TF1* Function = new TF1("f","gaus",-0.5,0.5);
			// Dummy parameters
			Function->SetParameter(0,30);
			Function->SetParameter(1,0.);
			Function->SetParameter(2,0.2);			
			double sigma[NMC+1];						

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

				//MCPlot[igen]->Fit(Function,"QR0");
				//sigma[igen] = Function->GetParameter(2); // gaussian sigma

				//TString SigmaTString = "#sigma = " + to_string_with_precision(sigma[igen],2); 
				//TLegendEntry* lSigma = legSigma->AddEntry(MCPlot[igen], SigmaTString,"");
				//lSigma->SetTextColor(MCColors[igen]);				

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

			gPad->RedrawAxis();												

			//----------------------------------------//	

			BeamOnStatOnly->Fit(Function,"QR0");
			sigma[NMC] = Function->GetParameter(2); // sigma

			legData->AddEntry(BeamOnShapeStat,"MicroBooNE Data","");
			legData->AddEntry(BeamOnShapeStat,Label,"");				
			legData->AddEntry(BeamOnShapeStat,"Stat#oplusShape","ep");
			legData->AddEntry(BeamOnNormOnly,"Norm","f");
			//legData->AddEntry(BeamOnShapeStat,"","");			
			//legData->AddEntry(BeamOnStatOnly,"#sigma_{Data} = " + to_string_with_precision(sigma[NMC],2),"");						

			TLegendEntry* lSigma = legSigma->AddEntry(BeamOnStatOnly,"#sigma_{Data} = " + to_string_with_precision(sigma[NMC],2) + " GeV/c","");
			lSigma->SetTextColor(kBlack);							

			leg->Draw();
			legData->Draw();
			legSigma->Draw();			

			TLatex *text = new TLatex();
			text->SetTextFont(FontStyle);
			text->SetTextSize(0.06);
			text->DrawLatexNDC(0.2, 0.94,  PanelLabels[iplot] + " " + LatexLabel[ Mapping[PlotNames[iplot]] ]);

			TLatex *textPrel = new TLatex();
			textPrel->SetTextFont(FontStyle);
			textPrel->SetTextSize(0.04);
			//textPrel->DrawLatexNDC(0.65, 0.68,"MicroBooNE Preliminary");										

			//----------------------------------------//

			PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance_PRL_Rename/"+ SaveFig[iplot] +".pdf");
			//PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_Neutrino2022_PublicNote/Figures/PRL_Fig3_"+PlotNames[iplot]+"_"+Runs[irun]+"_"+UBCodeVersion+".pdf");			
			//PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/Papadopoulou_MITThesis/templates/Figures/PRL_Fig3_"+PlotNames[iplot]+"_"+Runs[irun]+"_"+UBCodeVersion+".pdf");	
//			delete PlotCanvas;															

		} // End of the loop over the plots

		//----------------------------------------//					

	} // End of the loop over the runs	

} // End of the program 