#include <TFile.h>
#include <TF1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <iostream>
#include <vector>

#include  "./SecondaryCode/CenterAxisTitle.cpp"
#include "./SecondaryCode/SetOffsetAndSize.cpp"
#include "./SecondaryCode/ToString.cpp"
#include "./Constants.h"

using namespace std;
using namespace Constants;

void Create2DPlots() {

	TString PathToFiles = "myFiles/";

	TH2D::SetDefaultSumw2();

	vector<TString> PlotNames;

	PlotNames.push_back("RecoTrueDeltaPTPlot2D"); 
	PlotNames.push_back("RecoTrueDeltaAlphaTPlot2D"); 
	PlotNames.push_back("RecoTrueDeltaPhiTPlot2D"); 
	PlotNames.push_back("RecoTrueMuonMomentumPlot2D");
	PlotNames.push_back("RecoTrueMuonPhiPlot2D");
	PlotNames.push_back("RecoTrueMuonCosThetaPlot2D");
	PlotNames.push_back("RecoTrueProtonMomentumPlot2D");
	PlotNames.push_back("RecoTrueProtonPhiPlot2D");
	PlotNames.push_back("RecoTrueProtonCosThetaPlot2D");
	PlotNames.push_back("RecoTrueECalPlot2D"); 
	PlotNames.push_back("RecoTrueQ2Plot2D");

	const int N2DPlots = PlotNames.size();
	cout << "Number of 2D Plots = " << N2DPlots << endl;

	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); gStyle->SetTitleFont(FontStyle,"t"); SetOffsetAndSize();

	vector<TString> NameOfSamples;

	NameOfSamples.push_back("Overlay9");

	const int NSamples = NameOfSamples.size();
	TCanvas* PlotCanvas[NSamples][N2DPlots] = {}; TFile* FileSample[NSamples] = {};
	TH2D* Plots[NSamples][N2DPlots] = {};

	for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

		FileSample[WhichSample] = TFile::Open(PathToFiles+"/"+UBCodeVersion+"/CCQEAnalysis_"+NameOfSamples[WhichSample]+"_"+WhichRun+"_"+UBCodeVersion+".root");

	}

	for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

		for (int WhichPlot = 0; WhichPlot < N2DPlots; WhichPlot ++) {

			TString LogScale = "";

			PlotCanvas[WhichSample][WhichPlot] = new TCanvas(PlotNames[WhichPlot]+NameOfSamples[WhichSample],
					    PlotNames[WhichPlot]+NameOfSamples[WhichSample],205,34,1024,768);
			PlotCanvas[WhichSample][WhichPlot]->cd();
			Plots[WhichSample][WhichPlot] = (TH2D*)(FileSample[WhichSample]->Get("CC1p"+PlotNames[WhichPlot]));
			Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetZaxis()->SetLabelFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetZaxis()->SetLabelSize(0.03);
			double ScalingFactor = double(Plots[0][WhichPlot]->GetEntries()) / double(Plots[WhichSample][WhichPlot]->GetEntries());
			Plots[WhichSample][WhichPlot]->Scale(ScalingFactor);
			Plots[WhichSample][WhichPlot]->GetZaxis()->SetRangeUser(0,Plots[0][WhichPlot]->GetMaximum());

			CenterAxisTitle(Plots[WhichSample][WhichPlot]);
			PlotCanvas[WhichSample][WhichPlot]->SetLogz();
			Plots[WhichSample][WhichPlot]->Draw("colz");

			//PlotCanvas[WhichSample][WhichPlot]->SaveAs("/home/afroditi/Dropbox/Papers_Analyses/2019/TransverseVariables/Support/"
			//	+PlotNames[WhichPlot]+NameOfSamples[WhichSample]+".pdf");
			PlotCanvas[WhichSample][WhichPlot]->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/"+PlotNames[WhichPlot]+NameOfSamples[WhichSample]+LogScale+"_"+WhichRun+"_"+UBCodeVersion+".pdf");
			PlotCanvas[WhichSample][WhichPlot]->SaveAs("./myPlots/eps/"+UBCodeVersion+"/"+PlotNames[WhichPlot]+NameOfSamples[WhichSample]+LogScale+"_"+WhichRun+"_"+UBCodeVersion+".eps");
			//delete PlotCanvas[WhichSample][WhichPlot];

		} // End of the loop over the plots

	} // End of the loop over the samples

	//delete []PlotNames;

	TString FileName = "myMigrationMatrices/FileMigrationMatrices_"+WhichRun+"_"+UBCodeVersion+".root";
	TFile* FileMigrationMatrices = new TFile(FileName,"recreate");
	std::cout << std::endl << "File " << FileName << " has been created" << std::endl;
	FileMigrationMatrices->cd();
	FileMigrationMatrices->Write();

} // End of the program
