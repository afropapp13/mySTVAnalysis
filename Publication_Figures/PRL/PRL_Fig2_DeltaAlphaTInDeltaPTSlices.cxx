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

#include "../../../myClasses/Constants.h"

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

void PRL_Fig2_DeltaAlphaTInDeltaPTSlices() {

	//----------------------------------------//

	int DecimalAccuracy = 2;
	int FontStyle = 132;
	double TextSize = 0.06;	

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(4);			

	//----------------------------------------//

	vector<TString> PlotNames;  vector<TString> PanelLabels;
	PlotNames.push_back("DeltaAlphaTPlot");  PanelLabels.push_back("(a)");
	PlotNames.push_back("SerialDeltaAlphaT_DeltaPTPlot_0"); PanelLabels.push_back("(b)");
	PlotNames.push_back("SerialDeltaAlphaT_DeltaPTPlot_2");		PanelLabels.push_back("(c)");

	const int NPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << NPlots << endl;

	//----------------------------------------//

	// MC Samples to compare

	vector<TString> MCSampleBand; vector<TString> Label; vector<int> MCColors;  vector<int> LineStyle;

	MCSampleBand.push_back("G21NoFSI"); Label.push_back("G21 No FSI"); MCColors.push_back(OverlayColor); LineStyle.push_back(kSolid);
	MCSampleBand.push_back("G21hA"); Label.push_back("G21 hA");  MCColors.push_back(NuWroColor); LineStyle.push_back(kSolid);	
	MCSampleBand.push_back("SuSav2"); Label.push_back("G21 hN");	 MCColors.push_back(NEUTColor); LineStyle.push_back(kSolid);
	MCSampleBand.push_back("G21G4"); Label.push_back("G21 G4");  MCColors.push_back(GiBUUColor); LineStyle.push_back(kSolid);	

//	MCSampleBand.push_back("GiBUUNoFSI"); Label.push_back("GiB No FSI"); MCColors.push_back(NuWroColor); LineStyle.push_back(kDashed);
//	MCSampleBand.push_back("GiBUU"); Label.push_back("GiB FSI"); MCColors.push_back(NuWroColor); LineStyle.push_back(kSolid);
//	MCSampleBand.push_back("Genie_v3_0_6_NoFSI"); Label.push_back("G18 No FSI"); MCColors.push_back(OverlayColor); LineStyle.push_back(kDashed);	
//	MCSampleBand.push_back("OverlayGENIE"); Label.push_back("G18 FSI"); MCColors.push_back(OverlayColor); LineStyle.push_back(kSolid);
//	MCSampleBand.push_back("GiBUUTscaling"); Label.push_back("GiBUUTscaling");	
//	MCSampleBand.push_back("NEUT");  Label.push_back("NEUT");
//	MCSampleBand.push_back("NEUTv5401_RFG");  Label.push_back("NEUTv5401_RFG");	
//	MCSampleBand.push_back("Overlay9NuWro"); Label.push_back("NuWro");
//	MCSampleBand.push_back("GENIEv2"); Label.push_back("Gv2");
//	MCSampleBand.push_back("GENIEv2LFG"); Label.push_back("Gv2 LFG");
//	MCSampleBand.push_back("GENIEv2EffSF"); Label.push_back("Gv2 EffSF");		
//	MCSampleBand.push_back("Genie_v3_0_6_Out_Of_The_Box"); Label.push_back("G18 No Tune");							
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

		TString TxtName = "/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance/XSec/XSec_DeltaAlphaTInDeltaPT.txt";
		ofstream myTxtFile;
		myTxtFile.open(TxtName);		

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
			PlotCanvas->SetLeftMargin(0.19);
			PlotCanvas->SetRightMargin(0.04);				
			PlotCanvas->Draw();

			TLegend* leg = new TLegend(0.23,0.78,0.72,0.88);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.04);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(2);
			leg->SetMargin(0.13);

			TLegend* legData = new TLegend(0.22,0.68,0.77,0.78);
			legData->SetBorderSize(0);
			legData->SetTextSize(0.04);
			legData->SetTextFont(FontStyle);
			legData->SetNColumns(1);
			legData->SetMargin(0.08);						

			//----------------------------------------//

			// Get the shape + stat data plot & plot it

			TH1D* BeamOnShapeStat = (TH1D*)( fXSec->Get("StatShape_" + PlotNames[iplot]) );
			BeamOnShapeStat->GetYaxis()->SetRangeUser(XSecRange[ Mapping[PlotNames[iplot]] ].first, XSecRange[ Mapping[PlotNames[iplot]] ].second);

			BeamOnShapeStat->Draw("e1x0 same");	

			TH2D* Cov = (TH2D*)fXSec->Get("UnfCov_"+PlotNames[iplot]);	

			//----------------------------------------//
			//----------------------------------------//			

			// Data release

			TH1D* BeamOnFullUnc = (TH1D*)( fXSec->Get("FullUnc_" + PlotNames[iplot]) );	
			int NBins = BeamOnFullUnc->GetXaxis()->GetNbins();

			myTxtFile << "Cross Section \\delta alpha_{T} " << LatexLabel[ Mapping[ PlotNames[iplot] ] ] << endl;
			myTxtFile << "Bin # & Low bin edge & High bin edge & Cross Section [10^{-38} \\frac{cm^{2}}{GeV/c ^{40}Ar}] & Uncertainty [10^{-38} \\frac{cm^{2}}{GeV/c ^{40}Ar}]" << endl << endl;			

			for (int ibin = 1; ibin <= NBins; ibin++) {

				double BinLow = BeamOnFullUnc->GetBinLowEdge(ibin);
				double BinWidth = BeamOnFullUnc->GetBinWidth(ibin);		
				double BinHigh = BinLow + BinWidth;	
				double BinValue = BeamOnFullUnc->GetBinContent(ibin);
				double BinError = BeamOnFullUnc->GetBinError(ibin);				

				myTxtFile << ibin << std::setprecision(4) << " & " << BinLow << " & " << BinHigh << std::setprecision(8) << " & " << BinValue << " & " <<  BinError << endl;

			}
			
			myTxtFile << endl << endl;

			//----------------------------------------//			

			myTxtFile << "Unfolded Covariance Matrix \\delta\\alpha_{T} " << LatexLabel[ Mapping[ PlotNames[iplot] ] ] << endl;
			myTxtFile << "       ";
			for (int ybin = 1; ybin <= NBins; ybin++) { myTxtFile << "Bin " << ybin << "  "; }
			myTxtFile << endl;

			for (int xbin = 1; xbin <= NBins; xbin++) {

				myTxtFile << "Bin " << xbin << "  ";

				for (int ybin = 1; ybin <= NBins; ybin++) {	

					double CovBinValue = Cov->GetBinContent(xbin,ybin);
					myTxtFile << std::setprecision(8) << CovBinValue << "   ";

				}	
				
				myTxtFile << endl;

			}	

			myTxtFile << endl << endl;	

			//----------------------------------------//

			TH2D* Ac = (TH2D*)fXSec->Get("Ac_"+PlotNames[iplot]);
			myTxtFile << "Additional Smearing Matrix \\delta\\alpha_{T} " << LatexLabel[ Mapping[ PlotNames[iplot] ] ] << endl;
			myTxtFile << "       ";
			for (int ybin = 1; ybin <= NBins; ybin++) { myTxtFile << "Bin " << ybin << "  "; }
			myTxtFile << endl;

			for (int xbin = 1; xbin <= NBins; xbin++) {

				myTxtFile << "Bin " << xbin << "  ";				

				for (int ybin = 1; ybin <= NBins; ybin++) {	

					double AcBinValue = Ac->GetBinContent(xbin,ybin);
					myTxtFile << std::setprecision(8) << AcBinValue << "   ";

				}	
				
				myTxtFile << endl;

			}	

			myTxtFile << endl << endl;			

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
				MCPlot[igen]->Draw("same hist");	

				CalcChiSquared(MCPlot[igen],BeamOnShapeStat,Cov,Chi2[igen],Ndof[igen],pval[igen]);
				TString Chi2NdofAlt = " (" + to_string_with_precision(Chi2[igen],1) + "/" + TString(std::to_string(Ndof[igen])) +")";									

				TLegendEntry* lGenie = leg->AddEntry(MCPlot[igen],Label[igen] + Chi2NdofAlt,"l");
				lGenie->SetTextColor(MCColors[igen]);

			} // End of the loop over the generators

			//----------------------------------------//	

			TH1D* BeamOnNormOnly = (TH1D*)( fXSec->Get("NormOnly_" + PlotNames[iplot]) );
			BeamOnNormOnly->SetFillColorAlpha(kGray+1, 0.45);
			BeamOnNormOnly->SetLineColor(kGray+1);
			BeamOnNormOnly->SetMarkerColor(kGray+1);			
			BeamOnNormOnly->Draw("e2 same");	

			// Plot the stat+shape again 
			// And then the stat only & norm only on top of that

			BeamOnShapeStat->Draw("e1x0 same");

			TH1D* BeamOnStatOnly = (TH1D*)( fXSec->Get("StatOnly_" + PlotNames[iplot]) );
			BeamOnStatOnly->Draw("e1x0 same");								

			//----------------------------------------//	

			legData->AddEntry(BeamOnShapeStat,"MicroBooNE Data (Stat #oplus Shape)","ep");	
			legData->AddEntry(BeamOnNormOnly,"Norm","f");				

			leg->Draw();
			if (iplot == 0) { legData->Draw(); }			

			TLatex *text = new TLatex();
			text->SetTextFont(FontStyle);
			text->SetTextSize(0.06);
			text->DrawLatexNDC(0.2, 0.94,  PanelLabels[iplot] + " " + LatexLabel[ Mapping[PlotNames[iplot]] ] + ", Preliminary");	

			//----------------------------------------//

			// Legend & Run / POT

			double tor860_wcut = Fulltor860_wcut_Combined;
			TString Label = ToString(tor860_wcut)+" POT";	

			TLatex *textPOT = new TLatex();
			textPOT->SetTextFont(FontStyle);
			textPOT->SetTextSize(0.06);
			if (iplot == 0) { textPOT->DrawLatexNDC(0.7, 0.94,Label);	}	

			//----------------------------------------//

			// Transparency pad

			TPad* pad = new TPad("pad","pad",0.3,0.2,0.85,0.38, 21);
			if (iplot == 2) { pad = new TPad("pad","pad",0.2,0.52,0.75,0.7, 21); }

			pad->SetFillColor(kWhite); 
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
			if (iplot == 2) { hA->GetYaxis()->SetRangeUser(0.61,2.29); }																		

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

			//----------------------------------------//

			gPad->RedrawAxis();			

			PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance/Figures/PRL_Fig2_"+PlotNames[iplot]+"_"+Runs[irun]+"_"+UBCodeVersion+".pdf");
			PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_Neutrino2022_PublicNote/Figures/PRL_Fig2_"+PlotNames[iplot]+"_"+Runs[irun]+"_"+UBCodeVersion+".pdf");			
			PlotCanvas->SaveAs("/home/afroditi/Dropbox/Apps/Overleaf/Papadopoulou_MITThesis/templates/Figures/PRL_Fig2_"+PlotNames[iplot]+"_"+Runs[irun]+"_"+UBCodeVersion+".pdf");	
			delete PlotCanvas;															

		} // End of the loop over the plots

		//----------------------------------------//					

	} // End of the loop over the runs	

} // End of the program 