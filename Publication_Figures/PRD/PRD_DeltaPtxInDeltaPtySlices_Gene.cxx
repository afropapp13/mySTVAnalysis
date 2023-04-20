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

		{ "DeltaPtxPlot", "DeltaPtxPlot" },
		{ "SerialDeltaPtx_DeltaPtyPlot_0", "DeltaPtx_DeltaPty_Minus0_75ToMinus0_15Plot" },
		{ "SerialDeltaPtx_DeltaPtyPlot_1", "DeltaPtx_DeltaPty_Minus0_15To0_15Plot" },
		{ "SerialDeltaPtx_DeltaPtyPlot_2", "DeltaPtx_DeltaPty_0_15To0_45Plot" }

};

//----------------------------------------//

void PRD_DeltaPtxInDeltaPtySlices_Gene() {

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

	vector<TString> PlotNames; vector<TString> PanelLabels; vector<double> Min; vector<double> Max; vector<TString> Units; vector<TString> SaveFig;
	PlotNames.push_back("DeltaPtxPlot"); PanelLabels.push_back("(a)"); Min.push_back(0.); Max.push_back(39.9); Units.push_back("[$10^{-38}\\frac{cm^{2}}{(GeV/\\textit{c})\\,^{40}Ar}$]"); SaveFig.push_back("Fig88");
	PlotNames.push_back("SerialDeltaPtx_DeltaPtyPlot_0"); PanelLabels.push_back("(b)"); Min.push_back(0.); Max.push_back(15.49); Units.push_back("[$10^{-38}\\frac{cm^{2}}{(GeV/\\textit{c})^{2}\\,^{40}Ar}$]"); SaveFig.push_back("Fig89");
	PlotNames.push_back("SerialDeltaPtx_DeltaPtyPlot_1");PanelLabels.push_back("(c)"); Min.push_back(0.); Max.push_back(104.); Units.push_back("[$10^{-38}\\frac{cm^{2}}{(GeV/\\textit{c})^{2}\\,^{40}Ar}$]"); SaveFig.push_back("Fig90");	
	PlotNames.push_back("SerialDeltaPtx_DeltaPtyPlot_2"); PanelLabels.push_back("(d)"); Min.push_back(0.); Max.push_back(13.); Units.push_back("[$10^{-38}\\frac{cm^{2}}{(GeV/\\textit{c})^{2}\\,^{40}Ar}$]"); SaveFig.push_back("Fig91");

	const int NPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << NPlots << endl;

	//----------------------------------------//

	// MC Samples to compare

	vector<TString> MCSampleBand; vector<TString> Label; vector<int> MCColors;  vector<int> LineStyle;

//	MCSampleBand.push_back("GiBUUNoFSI"); Label.push_back("GiB No FSI"); MCColors.push_back(NuWroColor); LineStyle.push_back(kDashed);
//	MCSampleBand.push_back("Genie_v3_0_6_NoFSI"); Label.push_back("G18 No FSI"); MCColors.push_back(OverlayColor); LineStyle.push_back(kDashed);	
	MCSampleBand.push_back("OverlayGENIE"); Label.push_back("G18    "); MCColors.push_back(OverlayColor); LineStyle.push_back(G18LineStyle);
	MCSampleBand.push_back("GiBUU"); Label.push_back("GiBUU"); MCColors.push_back(GiBUUColor); LineStyle.push_back(GiBUULineStyle);	
//	MCSampleBand.push_back("GiBUUTscaling"); Label.push_back("GiBUUTscaling");	
	MCSampleBand.push_back("NEUT");  Label.push_back("NEUT"); MCColors.push_back(kMagenta-9); LineStyle.push_back(NEUTLineStyle);
	MCSampleBand.push_back("Overlay9NuWro");  Label.push_back("NuWro"); MCColors.push_back(NEUTColor); LineStyle.push_back(NuWroLineStyle);	
//	MCSampleBand.push_back("NEUTv5401_RFG");  Label.push_back("NEUTv5401_RFG");	
//	MCSampleBand.push_back("Overlay9NuWro"); Label.push_back("NuWro");
//	MCSampleBand.push_back("GENIEv2"); Label.push_back("Gv2");
//	MCSampleBand.push_back("GENIEv2LFG"); Label.push_back("Gv2 LFG");
//	MCSampleBand.push_back("GENIEv2EffSF"); Label.push_back("Gv2 EffSF");		
//	MCSampleBand.push_back("Genie_v3_0_6_Out_Of_The_Box"); Label.push_back("G18 No Tune");					
//	MCSampleBand.push_back("SuSav2"); Label.push_back("G21hN");
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

		// Data release

		TString TxtName = "/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance_PRD_Rename/PRD_XSec_DeltaPtxInDeltaPty.tex";
		ofstream myTxtFile;
		myTxtFile.open(TxtName);		

		//----------------------------------------//		

		// Loop over the plots

		for (int iplot = 0; iplot < NPlots; iplot ++) {

			//----------------------------------------//			

			// Canvas & legend

			TString CanvasName = "PRD_" + PlotNames[iplot]+"_"+Runs[irun];
			TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.14);
			PlotCanvas->SetTopMargin(0.1);
			PlotCanvas->SetLeftMargin(0.19);
			PlotCanvas->SetRightMargin(0.01);				
			PlotCanvas->Draw();

//			TLegend* legData = new TLegend(0.51,0.68,0.97,0.8);
			TLegend* legData = new TLegend(0.23,0.75,0.65,0.87);		
			legData->SetBorderSize(0);
			legData->SetTextSize(TextSize);
			legData->SetTextFont(FontStyle);
			legData->SetNColumns(1);
			legData->SetMargin(0.08);
			legData->SetFillStyle(0);		

			TLegend* legUnc = new TLegend(0.222,0.69,0.622,0.75);
			legUnc->SetBorderSize(0);
			legUnc->SetTextSize(TextSize);
			legUnc->SetTextFont(FontStyle);
			legUnc->SetNColumns(2);
			legUnc->SetMargin(0.2);
			legUnc->SetFillStyle(0);							

//			TLegend* leg = new TLegend(0.52,0.81,0.98,0.9);
			TLegend* leg = new TLegend(0.65,0.69,0.91,0.86);					
			leg->SetBorderSize(0);
			leg->SetTextSize(TextSize);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(1);
			leg->SetMargin(0.13);	
			leg->SetFillStyle(0);						

			//----------------------------------------//

			// Get the shape + stat data plot & plot it

			TH1D* BeamOnShapeStat = (TH1D*)( fXSec->Get("StatShape_" + PlotNames[iplot]) );
			BeamOnShapeStat->GetYaxis()->SetRangeUser(XSecRange[ Mapping[PlotNames[iplot]] ].first, XSecRange[ Mapping[PlotNames[iplot]] ].second);			
			BeamOnShapeStat->SetLineWidth(2);
			BeamOnShapeStat->Draw("e1x0 same");	

			TH2D* Cov = (TH2D*)fXSec->Get("UnfCov_"+PlotNames[iplot]);				

			//----------------------------------------//
			//----------------------------------------//			

			// Data release

			TH1D* BeamOnFullUnc = (TH1D*)( fXSec->Get("FullUnc_" + PlotNames[iplot]) );	
			int NBins = BeamOnFullUnc->GetXaxis()->GetNbins();

			TString LatexLabelString = "$"+LatexLabel[ Mapping[ PlotNames[iplot] ] ]+"$";
			LatexLabelString.ReplaceAll("#","\\").ReplaceAll(" ","\\,").ReplaceAll("\\deltap","\\delta p").ReplaceAll("/c","/\\textit{c}");
			myTxtFile << "\\begin{table}[H]" << endl;
			myTxtFile << "\\raggedright" << endl;	
			myTxtFile << "\\begin{adjustbox}{width=\\textwidth}" << endl;						
			myTxtFile << "\\small" << endl;
			myTxtFile << "\\begin{tabular}{ |c|c|c|c|c| }" << endl;	
			myTxtFile << "\\hline" << endl;						
			myTxtFile << "\\multicolumn{5}{|c|}{Cross Section $\\delta p_{T,x}$, " << LatexLabelString << "} \\\\" << endl;
			myTxtFile << "\\hline" << endl;
			myTxtFile << "\\hline" << endl;			
			myTxtFile << "Bin \\# & Low edge [GeV/$\\textit{c}$] & High edge [GeV/$\\textit{c}$] & Cross Section " << Units[iplot] << " & Uncertainty " << Units[iplot] << " \\\\" << endl;			
			myTxtFile << "\\hline" << endl;
			myTxtFile << "\\hline" << endl;				

			for (int ibin = 1; ibin <= NBins; ibin++) {

				double BinLow = BeamOnFullUnc->GetBinLowEdge(ibin);
				double BinWidth = BeamOnFullUnc->GetBinWidth(ibin);		
				double BinHigh = BinLow + BinWidth;	
				double BinValue = BeamOnFullUnc->GetBinContent(ibin);
				double BinError = BeamOnFullUnc->GetBinError(ibin);				

				myTxtFile << ibin << std::setprecision(4) << " & " << BinLow << " & " << BinHigh << std::setprecision(8) << " & " << BinValue << " & " <<  BinError << "\\\\" << endl;

			}
			
			myTxtFile << "\\hline" << endl;			
			myTxtFile << "\\end{tabular}" << endl;
			myTxtFile << "\\end{adjustbox}" << endl;		
			myTxtFile << "\\end{table}" << endl;				
			myTxtFile << endl << endl;

			//----------------------------------------//			

			// myTxtFile << "\\begin{table}[H]" << endl;
			// myTxtFile << "\\centering" << endl;	
			// myTxtFile << "\\begin{adjustbox}{width=\\textwidth}" << endl;		
			// myTxtFile << "\\small" << endl;						
			// myTxtFile << "\\begin{tabular}{ " << PrintMultipleTimes(NBins+1,"|c") << "| }" << endl;
			// myTxtFile << "\\hline" << endl;						
			// myTxtFile << "\\multicolumn{" << NBins+1 << "}{|c|}{Unfolded Covariance Matrix $\\delta p_{T,x}$, " << LatexLabelString << "} \\\\" << endl;
			// myTxtFile << "\\hline" << endl;
			// myTxtFile << "\\hline" << endl;
			// myTxtFile << "Units in " << Units[iplot] << "$^{2}$" << endl;			
			// for (int ybin = 1; ybin <= NBins; ybin++) { myTxtFile << " & Bin " << ybin; }
			// myTxtFile << "\\\\"<< endl;
			// myTxtFile << "\\hline"<< endl;

			// for (int xbin = 1; xbin <= NBins; xbin++) {

			// 	myTxtFile << "Bin " << xbin << "  ";

			// 	for (int ybin = 1; ybin <= NBins; ybin++) {	

			// 		double CovBinValue = Cov->GetBinContent(xbin,ybin);
			// 		myTxtFile << std::setprecision(6) << " & " << CovBinValue;

			// 	}	
				
			// 	myTxtFile << "\\\\" << endl;

			// }	

			// myTxtFile << "\\hline" << endl;			
			// myTxtFile << "\\end{tabular}" << endl;
			// myTxtFile << "\\end{adjustbox}" << endl;		
			// myTxtFile << "\\end{table}" << endl;				
			// myTxtFile << endl << endl;	

			//----------------------------------------//

			TH2D* Ac = (TH2D*)fXSec->Get("Ac_"+PlotNames[iplot]);			

			// myTxtFile << "\\begin{table}[H]" << endl;
			// myTxtFile << "\\centering" << endl;	
			// myTxtFile << "\\begin{adjustbox}{width=\\textwidth}" << endl;		
			// myTxtFile << "\\small" << endl;						
			// myTxtFile << "\\begin{tabular}{ " << PrintMultipleTimes(NBins+1,"|c") << "| }" << endl;
			// myTxtFile << "\\hline" << endl;						
			// myTxtFile << "\\multicolumn{" << NBins+1 << "}{|c|}{Additional Smearing Matrix ($A_{C}$) $\\delta p_{T,x}$, " << LatexLabelString << "} \\\\" << endl;
			// myTxtFile << "\\hline" << endl;
			// myTxtFile << "\\hline" << endl;
			// for (int ybin = 1; ybin <= NBins; ybin++) { myTxtFile << " & Bin " << ybin; }
			// myTxtFile << "\\\\"<< endl;
			// myTxtFile << "\\hline"<< endl;

			// for (int xbin = 1; xbin <= NBins; xbin++) {

			// 	myTxtFile << "Bin " << xbin << "  ";				

			// 	for (int ybin = 1; ybin <= NBins; ybin++) {	

			// 		double AcBinValue = Ac->GetBinContent(xbin,ybin);
			// 		myTxtFile << std::setprecision(6) << " & " << AcBinValue;

			// 	}	
				
			// 	myTxtFile << "\\\\" << endl;

			// }	

			// myTxtFile << "\\hline" << endl;			
			// myTxtFile << "\\end{tabular}" << endl;
			// myTxtFile << "\\end{adjustbox}" << endl;		
			// myTxtFile << "\\end{table}" << endl;				
			// myTxtFile << endl << endl;			

			//----------------------------------------//
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

			// Legend & Run / POT

			double tor860_wcut = Fulltor860_wcut_Combined;
			TString Label = ToString(tor860_wcut).ReplaceAll("e"," #times 10").ReplaceAll("+","^{")+"} POT";	

			TLatex *textPOT = new TLatex();
			textPOT->SetTextFont(FontStyle);
			textPOT->SetTextSize(0.06);
			//if (iplot == 0) { textPOT->DrawLatexNDC(0.7, 0.94,Label);	}																

			//----------------------------------------//	

			// Norm only

			TH1D* BeamOnNormOnly = (TH1D*)( fXSec->Get("NormOnly_" + PlotNames[iplot]) );
			BeamOnNormOnly->SetFillColor(kGray+1);
			BeamOnNormOnly->SetLineColor(kGray+1);
			BeamOnNormOnly->SetMarkerColor(kGray+1);						

			TLatex *text = new TLatex();
			text->SetTextFont(FontStyle);
			text->SetTextSize(0.06);
			text->DrawLatexNDC(0.2, 0.94, PanelLabels[iplot] + " " + LatexLabel[ Mapping[PlotNames[iplot]] ]);	

			//----------------------------------------//	

			// Norm only

			BeamOnNormOnly->Draw("e2 same");			

			// Plot the stat+shape again 
			// And then the stat only on top of that

			if (PlotNames[iplot] == "DeltaPTPlot") { BeamOnShapeStat->GetYaxis()->SetTitleOffset(1.1); }
			BeamOnShapeStat->GetYaxis()->SetRangeUser(Min[iplot],Max[iplot]);
			BeamOnShapeStat->Draw("e1x0 same");

			TH1D* BeamOnStatOnly = (TH1D*)( fXSec->Get("StatOnly_" + PlotNames[iplot]) );
			BeamOnStatOnly->SetLineWidth(2);			
			BeamOnStatOnly->Draw("e1x0 same");			

			gPad->RedrawAxis();
			
			//----------------------------------------//
			
			legData->AddEntry(BeamOnShapeStat,"MicroBooNE Data","");
			legData->AddEntry(BeamOnShapeStat,Label,"");			

			legUnc->AddEntry(BeamOnShapeStat,"Stat#oplusShape","ep");
			legUnc->AddEntry(BeamOnNormOnly,"Norm","f");				

			legData->Draw();
			legUnc->Draw();	
			leg->Draw();						

			//----------------------------------------//

			PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance_PRD_Rename/"+SaveFig[iplot]+".pdf");
			//PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance/Figures/PRD_Generators_"+PlotNames[iplot]+"_"+Runs[irun]+"_"+UBCodeVersion+".pdf");
			//delete PlotCanvas;	

			//----------------------------------------//

		} // End of the loop over the plots

		//----------------------------------------//					

	} // End of the loop over the runs	

} // End of the program 
