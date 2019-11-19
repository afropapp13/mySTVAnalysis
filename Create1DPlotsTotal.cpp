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
#include "./Constants.h"

using namespace std;
using namespace Constants;

void Create1DPlotsTotal() {

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;

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

	// -------------------------------------------------------------------------------------

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
	vector<vector<TH1D*> > PlotsTrueReco; PlotsTrueReco.clear();

	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();

	vector<TString> LabelsOfSamples;
	vector<TString> NameOfSamples;

	NameOfSamples.push_back("Overlay9"); // CV

	// Detector Systematics

	NameOfSamples.push_back("Overlay9_SCE");
	NameOfSamples.push_back("Overlay9_DLdown");

	const int NSamples = NameOfSamples.size();
	vector<TFile*> FileSample; FileSample.clear();
	vector<TFile*> TruthFileSample; TruthFileSample.clear();

	TString Name = "";
	TFile* FileEfficiences;

	for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

		FileSample.push_back(TFile::Open(PathToFiles+"/"+UBCodeVersion+"/CCQEStudies_"+NameOfSamples[WhichSample]+CutExtension+".root"));
		TruthFileSample.push_back(TFile::Open(PathToFiles+"/"+UBCodeVersion+"/TruthCCQEAnalysis_"+NameOfSamples[WhichSample]+"_"+WhichRun+"_"+UBCodeVersion+".root"));

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

		Name = "myEfficiencies/"+UBCodeVersion+"/FileEfficiences_"+NameOfSamples[WhichSample]+"_"+WhichRun+"_"+UBCodeVersion+".root";
		FileEfficiences = new TFile(Name,"recreate");

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){
	
			TCanvas* PlotCanvas = new TCanvas(NameOfSamples[WhichSample]+"_"+PlotNames[WhichPlot],NameOfSamples[WhichSample]+"_"+PlotNames[WhichPlot],205,34,1024,768);
			PlotCanvas->cd();

			TLegend* leg = new TLegend(0.1,0.92,0.9,1.);
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
			PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetRangeUser(0.,1.1*PlotsTrue[WhichSample][WhichPlot]->GetMaximum());
			PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(6);
			PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(0.8);
			PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(0.06);
			PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(0.04);
			PlotsTrue[WhichSample][WhichPlot]->GetYaxis()->SetTitle("# events");

			PlotsTrueReco[WhichSample][WhichPlot]->SetLineColor(kBlue);
			PlotsTrueReco[WhichSample][WhichPlot]->SetLineWidth(3);

			PlotsTrue[WhichSample][WhichPlot]->Draw();
			PlotsTrueReco[WhichSample][WhichPlot]->Draw("same");

			leg->AddEntry(PlotsTrue[WhichSample][WhichPlot],"True CC1p");
			leg->AddEntry(PlotsTrueReco[WhichSample][WhichPlot],"Candidate Reco CC1p");
			leg->Draw();		

			if (WhichSample == 0) 
				{ PlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/"+NameOfSamples[WhichSample]+"/"+PlotNames[WhichPlot]+"_"+WhichRun+"_"+UBCodeVersion+".pdf"); }
			delete PlotCanvas;

			// ---------------------------------------------------------------------------------------------------------------------------

			TCanvas* PlotEffCanvas = new TCanvas(NameOfSamples[WhichSample]+"_"+"Eff"+PlotNames[WhichPlot],NameOfSamples[WhichSample]+"_"+"Eff"+PlotNames[WhichPlot],205,34,1024,768);
			PlotEffCanvas->cd();

			PlotsTrueReco[WhichSample][WhichPlot]->Divide(PlotsTrue[WhichSample][WhichPlot]);
			TH1D* pEffPlot = PlotsTrueReco[WhichSample][WhichPlot];
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
			pEffPlot->GetYaxis()->SetRangeUser(0.,1.1*pEffPlot->GetMaximum());

			pEffPlot->Draw();
			FileEfficiences->cd();
			pEffPlot->Write();

			if (WhichSample == 0) 
				{ PlotEffCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/"+NameOfSamples[WhichSample]+"/Eff"+PlotNames[WhichPlot]+"_"+WhichRun+"_"+UBCodeVersion+".pdf"); }
			delete PlotEffCanvas;

		} // End of the loop over the plots

		FileEfficiences->Close();

		std::cout << std::endl << "Efficiency file " << Name << " created" << std::endl << std::endl;
		std::cout << std::endl << "--------------------------------------------------------------------------------------------------------------------------" << std::endl << std::endl;

	} // End of the loop over the samples

} // End of the program 
