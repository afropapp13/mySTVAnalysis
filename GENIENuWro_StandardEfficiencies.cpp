#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TGaxis.h>

#include <iostream>
#include <vector>

#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

void GENIENuWro_StandardEfficiencies() {

	//----------------------------------------//

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);	
	TGaxis::SetMaxDigits(4);
	TGaxis::SetExponentOffset(-0.05, 0., "y");	
	
	double text_size = 0.07;

	//----------------------------------------//

	TString cut = "_NoCuts_PID_NuScore";
	TString plot_name = "ProtonCosThetaPlot";
	TString run = "Combined";
	TString file_path = PathToFiles+"/";
	TString reco_file_path = PathToFiles+"/"+cut+"/";		

	//----------------------------------------//

	// GENIE file and plots

	TString genie_reco_file_name = "STVStudies_Overlay9_"+run+cut+".root";
	TFile*  genie_reco_file = new TFile(reco_file_path+genie_reco_file_name);
			
	TString genie_true_file_name = "TruthSTVAnalysis_Overlay9_"+run+"_"+UBCodeVersion+".root";
	TFile*  genie_true_file = new TFile(file_path+genie_true_file_name);

	TH1* genie_reco_plot = (TH1D*)(genie_reco_file->Get("CC1pTrue"+plot_name));
	TH1* genie_true_plot = (TH1D*)(genie_true_file->Get("True"+plot_name));	

	genie_reco_plot->Divide(genie_true_plot);
	genie_reco_plot->Scale(100.);		

	//----------------------------------------//

	// NuWro file and plots	

	TString nuwro_reco_file_name = "STVStudies_Overlay9NuWro_"+run+cut+".root";
	TFile*  nuwro_reco_file = new TFile(reco_file_path+nuwro_reco_file_name);
			
	TString nuwro_true_file_name = "TruthSTVAnalysis_Overlay9NuWro_"+run+"_"+UBCodeVersion+".root";
	TFile*  nuwro_true_file = new TFile(file_path+nuwro_true_file_name);

	TH1* nuwro_reco_plot = (TH1D*)(nuwro_reco_file->Get("CC1pTrue"+plot_name));
	TH1* nuwro_true_plot = (TH1D*)(nuwro_true_file->Get("True"+plot_name));	

	nuwro_reco_plot->Divide(nuwro_true_plot);	
	nuwro_reco_plot->Scale(100.);

	//----------------------------------------//

	TString canvas_name = "canvas_" + plot_name;
	TCanvas* tcanvas = new TCanvas(canvas_name,canvas_name,205,34,1024,768);
	tcanvas->cd();
	tcanvas->SetBottomMargin(0.15);
	tcanvas->SetLeftMargin(0.15);								

	TLegend* leg = new TLegend(0.25,0.92,0.9,1.);
	leg->SetBorderSize(0);
	leg->SetTextSize(text_size);
	leg->SetTextFont(FontStyle);
	leg->SetNColumns(2);
	leg->SetMargin(0.15);

	//----------------------------------------//

	// GENIE plot

	genie_reco_plot->SetLineColor(kBlack);
	genie_reco_plot->SetLineWidth(2);
	genie_reco_plot->SetMarkerStyle(20);
	genie_reco_plot->SetMarkerColor(kBlack);	
	genie_reco_plot->SetMarkerSize(2.);		

	genie_reco_plot->GetXaxis()->CenterTitle();
	genie_reco_plot->GetXaxis()->SetTitleFont(FontStyle);
	genie_reco_plot->GetXaxis()->SetLabelFont(FontStyle);
	genie_reco_plot->GetXaxis()->SetTitleSize(text_size);
	genie_reco_plot->GetXaxis()->SetLabelSize(text_size);
	genie_reco_plot->GetXaxis()->SetNdivisions(9);	

	genie_reco_plot->GetYaxis()->SetRangeUser(0.01,16);	
	genie_reco_plot->GetYaxis()->CenterTitle();
	genie_reco_plot->GetYaxis()->SetTitleFont(FontStyle);
	genie_reco_plot->GetYaxis()->SetLabelFont(FontStyle);
	genie_reco_plot->GetYaxis()->SetTitleSize(text_size);
	genie_reco_plot->GetYaxis()->SetLabelSize(text_size);
	genie_reco_plot->GetYaxis()->SetNdivisions(6);	
	genie_reco_plot->GetYaxis()->SetTitle("Efficiency [%]");		

	genie_reco_plot->Draw("ep");
	TLegendEntry* genie_leg = leg->AddEntry(genie_reco_plot,"G18","ep");	
	genie_leg->SetTextColor(kBlack);	

	//----------------------------------------//


	// NuWro plot

	nuwro_reco_plot->SetLineColor(kOrange+7);
	nuwro_reco_plot->SetLineWidth(2);
	nuwro_reco_plot->SetMarkerStyle(21);
	nuwro_reco_plot->SetMarkerSize(2.);
	nuwro_reco_plot->SetMarkerColor(kOrange+7);	

	nuwro_reco_plot->Draw("ep same");
	TLegendEntry* nuwro_leg = leg->AddEntry(nuwro_reco_plot,"NuWro","ep");	
	nuwro_leg->SetTextColor(kOrange+7);	

	//----------------------------------------//

	leg->Draw("same");

	//----------------------------------------//		

/*



					pEffPlot->SetMarkerColor(kBlack);
					pEffPlot->SetMarkerStyle(20);
					pEffPlot->SetMarkerSize(2.);
					pEffPlot->SetLineColor(kBlack);

					pEffPlot->Draw("e1");

					TLatex *textEff = new TLatex();
					textEff->SetTextFont(FontStyle);
					textEff->SetTextSize(TextSize);					
				
					TString CanvasEffPath = PlotPath+NameOfSamples[WhichSample]+"/";
					TString CanvasEffRatioName = "StandardEff"+PlotNamesClone[WhichPlot]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".pdf";
					PlotEffCanvas->SaveAs(CanvasEffPath + CanvasEffRatioName);

					delete PlotEffCanvas;
									
				}				

			} // End of the loop over the plots
*/

} // End of the program 
