#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TEfficiency.h>

#include <iostream>
#include <vector>

#include  "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/ToString.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/MakeMyPlotPretty.cpp"
#include "./Constants.h"

using namespace std;
using namespace Constants;

void DataDistributions() {

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;

//	TString PathToFiles = "./myFiles/";
//	TString PathToFiles = "../myAnalysis/OutputFiles";
	TString PathToFiles = "../myEvents/OutputFiles";

	// -------------------------------------------------------------------------------------

	int NEventsPassingSelectionCuts = 0;
	TString CutExtension = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
	VectorCuts.push_back("_NuScore");
	VectorCuts.push_back("_ThreePlaneLogChi2");
	VectorCuts.push_back("_MatchedFlash");
	VectorCuts.push_back("_Collinearity");

//	VectorCuts.push_back("_Chi2");
//	VectorCuts.push_back("_Distance");
//	VectorCuts.push_back("_Coplanarity");
//	VectorCuts.push_back("_TransImb");

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

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

//	vector<TCanvas*> PlotCanvas; PlotCanvas.clear();
	vector<TH1D*> PlotsReco; PlotsReco.clear();
	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();

	vector<TString> LabelsOfSamples;
	vector<TString> NameOfSamples;
	
	NameOfSamples.push_back("Run1Data9");

	const int NSamples = NameOfSamples.size();
	vector<TFile*> FileSample; FileSample.clear();

	vector<TEfficiency*> pEff;

	for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

		FileSample.push_back(TFile::Open(PathToFiles+"/"+UBCodeVersion+"/CCQEStudies_"+NameOfSamples[WhichSample]+CutExtension+".root"));

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

			TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
			PlotsReco.push_back(histReco);
		
		}

	}

	// Loop over the plots

	for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {
	
		TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot],PlotNames[WhichPlot],205,34,1024,768);
		PlotCanvas->cd();

		TLegend* leg = new TLegend(0.15,0.92,0.85,1.);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.07);
		leg->SetTextFont(FontStyle);
		leg->SetNColumns(3);

		MakeMyPlotPretty(PlotsReco[WhichPlot]);
		PlotsReco[WhichPlot]->SetLineColor(kBlack);
		PlotsReco[WhichPlot]->SetMarkerStyle(20);
		PlotsReco[WhichPlot]->SetMarkerSize(2.);
		PlotsReco[WhichPlot]->SetMarkerColor(kBlack);
		PlotsReco[WhichPlot]->GetYaxis()->SetRangeUser(0.,1.2*PlotsReco[WhichPlot]->GetMaximum());
		PlotsReco[WhichPlot]->GetYaxis()->SetTitle("# events");
		PlotsReco[WhichPlot]->Draw("same");

		leg->AddEntry(PlotsReco[WhichPlot],"BNB Run 1 Data","lep");
		leg->Draw();		

		PlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/"+NameOfSamples[0]+"/DataDistributions_"+PlotNames[WhichPlot]+"_"+WhichRun+"_"+UBCodeVersion+".pdf");
		PlotCanvas->SaveAs("./myPlots/eps/"+UBCodeVersion+"/"+NameOfSamples[0]+"/DataDistributions_"+PlotNames[WhichPlot]+"_"+WhichRun+"_"+UBCodeVersion+".eps");
		//delete PlotCanvas[WhichPlot];

	} // End of the loop over the plots

} // End of the program 
