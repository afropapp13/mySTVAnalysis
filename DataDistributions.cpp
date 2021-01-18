#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TEfficiency.h>

#include <iostream>
#include <vector>

#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

TString ToStringDouble(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

void DataDistributions(TString BeamOnSample) {

	// -------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;
	gStyle->SetOptStat(0);	

	// -------------------------------------------------------------------------------------

	int NEventsPassingSelectionCuts = 0;
	TString CutExtension = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();

	// v52
	VectorCuts.push_back("");
	VectorCuts.push_back("_PID");
	VectorCuts.push_back("_NuScore");

	/*
	// up to v43
	VectorCuts.push_back("");
	VectorCuts.push_back("_NuScore");
	VectorCuts.push_back("_ThreePlaneLogChi2");
	VectorCuts.push_back("_Collinearity");
	*/

	int NCuts = (int)(VectorCuts.size());	

	for (int i = 0; i < NCuts; i++) {

		CutExtension = CutExtension + VectorCuts[i];

	}

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
//	PlotNames.push_back("ECalPlot");
//	PlotNames.push_back("EQEPlot");
//	PlotNames.push_back("Q2Plot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ----------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");
//	Runs.push_back("Run2");
	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// ----------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {


		vector<TH1D*> PlotsReco; PlotsReco.clear();
		gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");// SetOffsetAndSize();

		vector<TString> LabelsOfSamples;
		vector<TString> NameOfSamples;
	
		NameOfSamples.push_back("BeamOn9");

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();

		vector<TEfficiency*> pEff;

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			FileSample.push_back(TFile::Open(PathToFiles+"/"+CutExtension+"/STVStudies_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+CutExtension+".root"));

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
				PlotsReco.push_back(histReco);
		
			}

		}

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {
	
			TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot],PlotNames[WhichPlot],205,34,1024,768);
			PlotCanvas->cd();

			TLegend* leg = new TLegend(0.12,0.92,0.85,1.);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.07);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(3);

//			MakeMyPlotPretty(PlotsReco[WhichPlot]);
			PlotsReco[WhichPlot]->SetLineColor(kBlack);
			PlotsReco[WhichPlot]->SetMarkerStyle(20);
			PlotsReco[WhichPlot]->SetMarkerSize(2.);
			PlotsReco[WhichPlot]->SetMarkerColor(kBlack);

			PlotsReco[WhichPlot]->GetXaxis()->CenterTitle();
			PlotsReco[WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
			PlotsReco[WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
			PlotsReco[WhichPlot]->GetXaxis()->SetNdivisions(6);									
			
			PlotsReco[WhichPlot]->GetYaxis()->SetRangeUser(0.,1.2*PlotsReco[WhichPlot]->GetMaximum());
			PlotsReco[WhichPlot]->GetYaxis()->SetTitle("# events");
			PlotsReco[WhichPlot]->GetYaxis()->CenterTitle();	
			PlotsReco[WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
			PlotsReco[WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);			
			PlotsReco[WhichPlot]->GetYaxis()->SetNdivisions(6);			
			
			PlotsReco[WhichPlot]->Draw("ex0 same");

			double tor860_wcut = -99.;
				
			if (Runs[WhichRun] == "Run1") { tor860_wcut = tor860_wcut_Run1; }
			if (Runs[WhichRun] == "Run2") { tor860_wcut = tor860_wcut_Run2; }
			if (Runs[WhichRun] == "Run3") { tor860_wcut = tor860_wcut_Run3; }
			if (Runs[WhichRun] == "Run4") { tor860_wcut = tor860_wcut_Run4; }
			if (Runs[WhichRun] == "Run5") { tor860_wcut = tor860_wcut_Run5; }

			TString Label = ToStringDouble(tor860_wcut)+" POT";

			TString DataLabel = "MicroBooNE " + Runs[WhichRun] + " " + Label;
			leg->AddEntry(PlotsReco[WhichPlot],DataLabel,"ep");
			leg->Draw();		

			PlotCanvas->SaveAs(PlotPath+NameOfSamples[0]+"/DataDistributions_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");

			delete PlotCanvas;

		} // End of the loop over the plots

	} // End of the loop over the runs

} // End of the program 
