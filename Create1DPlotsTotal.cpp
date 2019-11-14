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

	for (int i = 0; i < NCuts; i++) {

		CutExtension = CutExtension + VectorCuts[i];

	}

	// -------------------------------------------------------------------------------------

	int FontStyle = 132;

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

	vector<TCanvas*> PlotCanvas; PlotCanvas.clear();
	vector<TCanvas*> PlotEffCanvas; PlotEffCanvas.clear();
	vector<TH1D*> PlotsTrue; PlotsTrue.clear();
	vector<TH1D*> PlotsTrueReco; PlotsTrueReco.clear();
//	vector<TH1D*> PlotsReco; PlotsReco.clear();
	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();

	vector<TString> LabelsOfSamples;
	vector<TString> NameOfSamples;

	NameOfSamples.push_back("Overlay9");
//	NameOfSamples.push_back("Overlay9_SCE");
//	NameOfSamples.push_back("Overlay9_DLdown");

	TString Name = "myEfficiencies/"+UBCodeVersion+"/FileEfficiences_"+NameOfSamples[0]+"_"+WhichRun+"_"+UBCodeVersion+".root";
	TFile* FileEfficiences = new TFile(Name,"recreate");

	const int NSamples = NameOfSamples.size();
	vector<TFile*> FileSample; FileSample.clear();
	vector<TFile*> TruthFileSample; TruthFileSample.clear();

	vector<TEfficiency*> pEff;
	vector<TH1D*> pEffPlot;

	for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

//		FileSample.push_back(TFile::Open(PathToFiles+"/"+UBCodeVersion+"/CCQEAnalysis_"+NameOfSamples[WhichSample]+"_"+WhichRun+"_"+UBCodeVersion+".root"));
		FileSample.push_back(TFile::Open(PathToFiles+"/"+UBCodeVersion+"/CCQEStudies_"+NameOfSamples[WhichSample]+CutExtension+".root"));
		TruthFileSample.push_back(TFile::Open(PathToFiles+"/"+UBCodeVersion+"/TruthCCQEAnalysis_"+NameOfSamples[WhichSample]+"_"+WhichRun+"_"+UBCodeVersion+".root"));

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

			TH1D* histTrue = (TH1D*)(TruthFileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
			PlotsTrue.push_back(histTrue);

//			TH1D* histTrueReco = (TH1D*)(FileSample[WhichSample]->Get("RecoTrue"+PlotNames[WhichPlot]));
//			TH1D* histTrueReco = (TH1D*)(FileSample[WhichSample]->Get("CC1pRecoTrue"+PlotNames[WhichPlot]));

//			TH1D* histTrueReco = (TH1D*)(FileSample[WhichSample]->Get("CC1pTrue"+PlotNames[WhichPlot]));

			TH1D* histTrueReco = (TH1D*)(FileSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]));

//			TH1D* histTrueReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));

			PlotsTrueReco.push_back(histTrueReco);

//			TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
//			PlotsReco.push_back(histReco);
		
		}

	}

	// Loop over the plots

	for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){
	
		PlotCanvas.push_back(new TCanvas(PlotNames[WhichPlot],PlotNames[WhichPlot],205,34,1024,768));
		PlotCanvas[WhichPlot]->cd();

		TLegend* leg = new TLegend(0.1,0.92,0.9,1.);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.07);
		leg->SetTextFont(FontStyle);
		leg->SetNColumns(2);
		leg->SetMargin(0.15);

		PlotsTrue[WhichPlot]->SetLineColor(kRed);
		PlotsTrue[WhichPlot]->SetLineWidth(3);
		PlotsTrue[WhichPlot]->GetXaxis()->CenterTitle();
		PlotsTrue[WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
		PlotsTrue[WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
		PlotsTrue[WhichPlot]->GetXaxis()->SetTitleSize(0.06);
		PlotsTrue[WhichPlot]->GetXaxis()->SetLabelSize(0.04);
		PlotsTrue[WhichPlot]->GetXaxis()->SetTitleOffset(0.7);
		PlotsTrue[WhichPlot]->GetXaxis()->SetNdivisions(5);

		PlotsTrue[WhichPlot]->GetYaxis()->CenterTitle();
		PlotsTrue[WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
		PlotsTrue[WhichPlot]->GetYaxis()->SetTitleSize(0.12);
		PlotsTrue[WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
		PlotsTrue[WhichPlot]->GetYaxis()->SetRangeUser(0.,1.1*PlotsTrue[WhichPlot]->GetMaximum());
		PlotsTrue[WhichPlot]->GetYaxis()->SetNdivisions(6);
		PlotsTrue[WhichPlot]->GetYaxis()->SetTitleOffset(0.8);
		PlotsTrue[WhichPlot]->GetYaxis()->SetTitleSize(0.06);
		PlotsTrue[WhichPlot]->GetYaxis()->SetLabelSize(0.04);
		PlotsTrue[WhichPlot]->GetYaxis()->SetTitle("# events");

		PlotsTrueReco[WhichPlot]->SetLineColor(kBlue);
		PlotsTrueReco[WhichPlot]->SetLineWidth(3);

//		PlotsReco[WhichPlot]->SetLineColor(kBlack);
//		PlotsReco[WhichPlot]->SetLineWidth(3);

		PlotsTrue[WhichPlot]->Draw();
		PlotsTrueReco[WhichPlot]->Draw("same");
//		PlotsReco[WhichPlot]->Draw("same");

		leg->AddEntry(PlotsTrue[WhichPlot],"True CC1p");
		leg->AddEntry(PlotsTrueReco[WhichPlot],"Candidate Reco CC1p");

//		leg->AddEntry(PlotsReco[WhichPlot],"Reco");
		leg->Draw();		

		//PlotCanvas[WhichPlot]->SaveAs("/home/afroditi/Dropbox/Papers_Analyses/2019/TransverseVariables/Support/"+PlotNames[WhichPlot]+".pdf");
		PlotCanvas[WhichPlot]->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/"+NameOfSamples[0]+"/"+PlotNames[WhichPlot]+"_"+WhichRun+"_"+UBCodeVersion+".pdf");
		PlotCanvas[WhichPlot]->SaveAs("./myPlots/eps/"+UBCodeVersion+"/"+NameOfSamples[0]+"/"+PlotNames[WhichPlot]+"_"+WhichRun+"_"+UBCodeVersion+".eps");
//		delete PlotCanvas[WhichPlot];

		// ---------------------------------------------------------------------------------------------------------------------------


		//if (PlotNames[WhichPlot] == "DeltaPTPlot" || PlotNames[WhichPlot] == "DeltaAlphaTPlot" || PlotNames[WhichPlot] == "DeltaPhiTPlot") {

			PlotEffCanvas.push_back(new TCanvas("Eff"+PlotNames[WhichPlot],"Eff"+PlotNames[WhichPlot],205,34,1024,768));
			PlotEffCanvas[WhichPlot]->cd();

/*			pEff.push_back(new TEfficiency(*PlotsTrueReco[WhichPlot],*PlotsTrue[WhichPlot]));
			pEff[WhichPlot]->SetLineWidth(3);

			pEff[WhichPlot]->GetTotalHistogram()->GetXaxis()->CenterTitle();
			pEff[WhichPlot]->GetTotalHistogram()->GetXaxis()->SetTitleFont(FontStyle);
			pEff[WhichPlot]->GetTotalHistogram()->GetXaxis()->SetLabelFont(FontStyle);
			pEff[WhichPlot]->GetTotalHistogram()->GetXaxis()->SetTitle(PlotsTrue[WhichPlot]->GetXaxis()->GetTitle());
			pEff[WhichPlot]->GetTotalHistogram()->GetXaxis()->SetTitleSize(0.06);
			pEff[WhichPlot]->GetTotalHistogram()->GetXaxis()->SetLabelSize(0.04);
			pEff[WhichPlot]->GetTotalHistogram()->GetXaxis()->SetTitleOffset(0.7);
			pEff[WhichPlot]->GetTotalHistogram()->GetXaxis()->SetNdivisions(5);

			pEff[WhichPlot]->GetTotalHistogram()->GetYaxis()->CenterTitle();
			pEff[WhichPlot]->GetTotalHistogram()->GetYaxis()->SetTitleFont(FontStyle);
			pEff[WhichPlot]->GetTotalHistogram()->GetYaxis()->SetTitleSize(0.12);
			pEff[WhichPlot]->GetTotalHistogram()->GetYaxis()->SetLabelFont(FontStyle);
			pEff[WhichPlot]->GetTotalHistogram()->GetYaxis()->SetRangeUser(0.,25.);
			pEff[WhichPlot]->GetTotalHistogram()->GetYaxis()->SetNdivisions(6);
			pEff[WhichPlot]->GetTotalHistogram()->GetYaxis()->SetTitleOffset(0.45);
			pEff[WhichPlot]->GetTotalHistogram()->GetYaxis()->SetTitleSize(0.08);
			pEff[WhichPlot]->GetTotalHistogram()->GetYaxis()->SetLabelSize(0.04);
			pEff[WhichPlot]->GetTotalHistogram()->GetYaxis()->SetTitle("Efficiency");

			pEff[WhichPlot]->Draw();
			FileEfficiences->cd();
			pEff[WhichPlot]->Write();*/

			PlotsTrueReco[WhichPlot]->Divide(PlotsTrue[WhichPlot]);
			pEffPlot.push_back(PlotsTrueReco[WhichPlot]);
			pEffPlot[WhichPlot]->SetLineWidth(3);
			pEffPlot[WhichPlot]->SetLineColor(kBlack);
			pEffPlot[WhichPlot]->SetMarkerStyle(20);

			pEffPlot[WhichPlot]->GetXaxis()->CenterTitle();
			pEffPlot[WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
			pEffPlot[WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
			pEffPlot[WhichPlot]->GetXaxis()->SetTitle(PlotsTrue[WhichPlot]->GetXaxis()->GetTitle());
			pEffPlot[WhichPlot]->GetXaxis()->SetTitleSize(0.06);
			pEffPlot[WhichPlot]->GetXaxis()->SetLabelSize(0.04);
			pEffPlot[WhichPlot]->GetXaxis()->SetTitleOffset(0.7);
			pEffPlot[WhichPlot]->GetXaxis()->SetNdivisions(5);

			pEffPlot[WhichPlot]->GetYaxis()->CenterTitle();
			pEffPlot[WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
			pEffPlot[WhichPlot]->GetYaxis()->SetTitleSize(0.12);
			pEffPlot[WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
//			pEffPlot[WhichPlot]->GetYaxis()->SetRangeUser(0.,25.);
			pEffPlot[WhichPlot]->GetYaxis()->SetNdivisions(6);
			pEffPlot[WhichPlot]->GetYaxis()->SetTitleOffset(0.45);
			pEffPlot[WhichPlot]->GetYaxis()->SetTitleSize(0.08);
			pEffPlot[WhichPlot]->GetYaxis()->SetLabelSize(0.04);
			pEffPlot[WhichPlot]->GetYaxis()->SetTitle("Efficiency");
			pEffPlot[WhichPlot]->GetYaxis()->SetRangeUser(0.,0.5);

			pEffPlot[WhichPlot]->Draw();
			FileEfficiences->cd();
			pEffPlot[WhichPlot]->Write();

			//PlotEffCanvas[WhichPlot]->SaveAs("/home/afroditi/Dropbox/Papers_Analyses/2019/TransverseVariables/Support/Eff"+PlotNames[WhichPlot]+".pdf");
			PlotEffCanvas[WhichPlot]->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/"+NameOfSamples[0]+"/Eff"+PlotNames[WhichPlot]+"_"+WhichRun+"_"+UBCodeVersion+".pdf");
			PlotEffCanvas[WhichPlot]->SaveAs("./myPlots/eps/"+UBCodeVersion+"/"+NameOfSamples[0]+"/Eff"+PlotNames[WhichPlot]+"_"+WhichRun+"_"+UBCodeVersion+".eps");
//			delete PlotEffCanvas[WhichPlot];

		//}

	} // End of the loop over the plots

	FileEfficiences->Close();

	std::cout << std::endl << "Efficiency file " << Name << " created" << std::endl << std::endl;

} // End of the program 
