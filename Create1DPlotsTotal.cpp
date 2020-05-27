#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TLatex.h>

#include <iostream>
#include <vector>

//#include  "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
//#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
//#include "/home/afroditi/Dropbox/PhD/Secondary_Code/ToString.cpp"
#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

void Create1DPlotsTotal(TString OverlaySample) {

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;
	gStyle->SetOptStat(0);	

	TString PathToFiles = "../myEvents/OutputFiles";

	// -------------------------------------------------------------------------------------

	int NEventsPassingSelectionCuts = 0;
	TString CutExtension = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
	VectorCuts.push_back("_NuScore");
	VectorCuts.push_back("_ThreePlaneLogChi2");
	VectorCuts.push_back("_Collinearity");

	int NCuts = (int)(VectorCuts.size());	
	for (int i = 0; i < NCuts; i++) { CutExtension = CutExtension + VectorCuts[i]; }

	// -------------------------------------------------------------------------------------

	PlotNames.push_back("DeltaPTPlot"); 
	PlotNames.push_back("DeltaAlphaTPlot"); 
	PlotNames.push_back("DeltaPhiTPlot");
	PlotNames.push_back("MuonMomentumPlot"); 
	PlotNames.push_back("MuonCosThetaPlot"); 
	PlotNames.push_back("MuonPhiPlot");
	PlotNames.push_back("ProtonMomentumPlot"); 
	PlotNames.push_back("ProtonCosThetaPlot"); 
	PlotNames.push_back("ProtonPhiPlot");
	PlotNames.push_back("ECalPlot");
	PlotNames.push_back("EQEPlot");
	PlotNames.push_back("Q2Plot");

	// -------------------------------------------------------------------------------------

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ------------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// ----------------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
		vector<vector<TH1D*> > PlotsTrueReco; PlotsTrueReco.clear();

		gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); //SetOffsetAndSize();

		vector<TString> LabelsOfSamples;
		vector<TString> NameOfSamples;

		NameOfSamples.push_back("Overlay9");

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();
		vector<TFile*> TruthFileSample; TruthFileSample.clear();

		TString Name = "";

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			FileSample.push_back(TFile::Open(PathToFiles+"/"+UBCodeVersion+"/"+CutExtension+"/STVStudies_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+CutExtension+".root"));
			TruthFileSample.push_back(TFile::Open(PathToFiles+"/"+UBCodeVersion+"/TruthSTVAnalysis_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root"));

			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();
			vector<TH1D*> CurrentPlotsTrueReco; CurrentPlotsTrueReco.clear();

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

				TH1D* histTrue = (TH1D*)(TruthFileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
				CurrentPlotsTrue.push_back(histTrue);

				TH1D* histTrueReco = (TH1D*)(FileSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]));
				CurrentPlotsTrueReco.push_back(histTrueReco);
		
			}

			PlotsTrue.push_back(CurrentPlotsTrue);
			PlotsTrueReco.push_back(CurrentPlotsTrueReco);

		}

		// Loop over the samples

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			Name = "myEfficiencies/"+UBCodeVersion+"/FileEfficiences_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
			TFile* FileEfficiences = new TFile(Name,"recreate");

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){
	
				TCanvas* PlotCanvas = new TCanvas(NameOfSamples[WhichSample]+"_"+PlotNames[WhichPlot],NameOfSamples[WhichSample]+"_"+PlotNames[WhichPlot],205,34,1024,768);
				PlotCanvas->cd();

				TLegend* leg = new TLegend(0.15,0.92,0.9,1.);
				leg->SetBorderSize(0);
				leg->SetTextSize(0.07);
				leg->SetTextFont(FontStyle);
				leg->SetNColumns(2);
				leg->SetMargin(0.15);

				PlotsTrue[WhichSample][WhichPlot]->SetLineColor(kRed);
				PlotsTrue[WhichSample][WhichPlot]->SetLineWidth(3);
				PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->CenterTitle();
				PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
				PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
				PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->SetTitleSize(0.06);
				PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(0.04);
				PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->SetTitleOffset(0.7);
				PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->SetNdivisions(5);

				PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->CenterTitle();
				PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
				PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(0.12);
				PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
				PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetRangeUser(0.,1.2*PlotsTrue[WhichSample][WhichPlot]->GetMaximum());
				PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(6);
				PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(0.8);
				PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(0.06);
				PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(0.04);
				PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetTitle("# events");

				PlotsTrueReco[WhichSample][WhichPlot]->SetLineColor(kBlue);
				PlotsTrueReco[WhichSample][WhichPlot]->SetLineWidth(3);

				// ----------------------------------------------------------------------------------------

				// We need the raw number of events, not the POT normalized one
				double tor860_wcut = -99.;
				if (Runs[WhichRun] == "Run1") { tor860_wcut = tor860_wcut_Run1; }
//				if (Runs[WhichRun] == "Run2") { tor860_wcut = tor860_wcut_Run2; }
//				if (Runs[WhichRun] == "Run3") { tor860_wcut = tor860_wcut_Run3; }
//				if (Runs[WhichRun] == "Run4") { tor860_wcut = tor860_wcut_Run4; }
//				if (Runs[WhichRun] == "Run5") { tor860_wcut = tor860_wcut_Run5; }

				double OverlayPOT = -99.;

				if (string(NameOfSamples[WhichSample]).find("Overlay") != std::string::npos) {

						TString PathToPOTFile = "../myEvents/mySamples/"+UBCodeVersion+"/PreSelection_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+"_POT.root";
						TFile* POTFile = TFile::Open(PathToPOTFile,"readonly");
						TH1D* POTCountHist = (TH1D*)(POTFile->Get("POTCountHist"));
						OverlayPOT = POTCountHist->GetBinContent(1);
						POTFile->Close();
				}

				PlotsTrueReco[WhichSample][WhichPlot]->Scale(OverlayPOT/tor860_wcut);

				// ----------------------------------------------------------------------------------------

				PlotsTrue[WhichSample][WhichPlot]->Draw();
				PlotsTrueReco[WhichSample][WhichPlot]->Draw("same");

				leg->AddEntry(PlotsTrue[WhichSample][WhichPlot],"True CC1p");
				leg->AddEntry(PlotsTrueReco[WhichSample][WhichPlot],"Reco CC1p");
				leg->Draw();

				TLatex *text = new TLatex();
				text->SetTextFont(FontStyle);
				text->SetTextSize(0.08);
				text->DrawTextNDC(0.14, 0.8, Runs[WhichRun]);		

				if (WhichSample == 0 && OverlaySample == "") 
					{ PlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/"+NameOfSamples[WhichSample]+"/"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".pdf"); }
				delete PlotCanvas;

				// ---------------------------------------------------------------------------------------------------------------------------

				TCanvas* PlotEffCanvas = new TCanvas(NameOfSamples[WhichSample]+"_"+"Eff"+PlotNames[WhichPlot],NameOfSamples[WhichSample]+"_"+"Eff"+PlotNames[WhichPlot],205,34,1024,768);
				PlotEffCanvas->cd();

				TH1D* pEffPlot = (TH1D*)PlotsTrueReco[WhichSample][WhichPlot]->Clone();
				pEffPlot->Divide(PlotsTrue[WhichSample][WhichPlot]);
				pEffPlot->SetLineWidth(3);
				pEffPlot->SetLineColor(kBlack);
				pEffPlot->SetMarkerStyle(20);

				pEffPlot->GetXaxis()->CenterTitle();
				pEffPlot->GetXaxis()->SetTitleFont(FontStyle);
				pEffPlot->GetXaxis()->SetLabelFont(FontStyle);
				pEffPlot->GetXaxis()->SetTitle(PlotsTrue[WhichSample][WhichPlot]->GetXaxis()->GetTitle());
				pEffPlot->GetXaxis()->SetTitleSize(0.06);
				pEffPlot->GetXaxis()->SetLabelSize(0.04);
				pEffPlot->GetXaxis()->SetTitleOffset(0.7);
				pEffPlot->GetXaxis()->SetNdivisions(5);

				pEffPlot->GetYaxis()->CenterTitle();
				pEffPlot->GetYaxis()->SetTitleFont(FontStyle);
				pEffPlot->GetYaxis()->SetTitleSize(0.12);
				pEffPlot->GetYaxis()->SetLabelFont(FontStyle);
				pEffPlot->GetYaxis()->SetNdivisions(6);
				pEffPlot->GetYaxis()->SetTitleOffset(0.45);
				pEffPlot->GetYaxis()->SetTitleSize(0.08);
				pEffPlot->GetYaxis()->SetLabelSize(0.04);
				pEffPlot->GetYaxis()->SetTitle("Efficiency");
				pEffPlot->GetYaxis()->SetRangeUser(0.,1.2*pEffPlot->GetMaximum());

				pEffPlot->Draw();

				TLatex *textEff = new TLatex();
				textEff->SetTextFont(FontStyle);
				textEff->SetTextSize(0.08);
				textEff->DrawTextNDC(0.14, 0.8, Runs[WhichRun]);

				FileEfficiences->cd();
				pEffPlot->Write();

				if (WhichSample == 0 && OverlaySample == "") 
					{ PlotEffCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/"+NameOfSamples[WhichSample]+"/Eff"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".pdf"); }
				delete PlotEffCanvas;

			} // End of the loop over the plots

			FileEfficiences->Close();

			std::cout << std::endl << "Efficiency file " << Name << " created" << std::endl << std::endl;
			std::cout << std::endl << "---------------------------------------------------------------------------------------------------------" << std::endl << std::endl;

		} // End of the loop over the samples

	} // End of the loop over the runs	

} // End of the program 
