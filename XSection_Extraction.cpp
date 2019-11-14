#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLine.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include  "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/ToString.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/MakeMyPlotPretty.cpp"

#include "./Constants.h"

using namespace std;
using namespace Constants;

void XSection_Extraction() {

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;

	TString Subtract = "";
//	TString Subtract = "_BUnsubtracted";

	TString PathToFiles = "../myEvents/OutputFiles/";

	TString PathToExtractedXSec = "myXSec/";

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
//	PlotNames.push_back("ECalPlot"); 
//	PlotNames.push_back("Q2Plot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	vector<TCanvas*> PlotCanvas; PlotCanvas.clear();

	vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
	vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
	vector<vector<TH1D*> > PlotsBkgReco; PlotsBkgReco.clear();
	vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();
	vector<vector<TH1D*> > PlotsTEfficiency; PlotsTEfficiency.clear();

	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();

	vector<TString> LabelsOfSamples;
	vector<TString> NameOfSamples;
	
	NameOfSamples.push_back("Overlay9"); 
//	NameOfSamples.push_back("Overlay9_SCE");
//	NameOfSamples.push_back("Overlay9_DLdown");
	NameOfSamples.push_back("Run1Data9"); NameOfSamples.push_back("ExtBNB9"); NameOfSamples.push_back("OverlayDirt9"); NameOfSamples.push_back("Genie");

	TFile* FileEfficiences = new TFile("myEfficiencies/"+UBCodeVersion+"/FileEfficiences_"+NameOfSamples[0]+"_"+WhichRun+"_"+UBCodeVersion+".root","readonly");

	const int NSamples = NameOfSamples.size();
	vector<TFile*> FileSample; FileSample.clear();

	for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

//		FileSample.push_back(TFile::Open(PathToFiles+UBCodeVersion+"/CCQEAnalysis_"+NameOfSamples[WhichSample]+"_"+WhichRun+"_"+UBCodeVersion+".root"));
		if (NameOfSamples[WhichSample] != "Genie") { FileSample.push_back(TFile::Open(PathToFiles+"/"+UBCodeVersion+"/CCQEStudies_"+NameOfSamples[WhichSample]+CutExtension+".root")); }
		else { FileSample.push_back(TFile::Open("myFiles/"+UBCodeVersion+"/CCQEAnalysis_"+NameOfSamples[WhichSample]+"_"+UBCodeVersion+".root")); }

		vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();
		vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();
		vector<TH1D*> CurrentPlotsBkgReco; CurrentPlotsBkgReco.clear();
		vector<TH1D*> CurrentPlotsCC1pReco; CurrentPlotsCC1pReco.clear();
		vector<TH1D*> CurrentPlotsTEfficiency; CurrentPlotsTEfficiency.clear();

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

			TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
			CurrentPlotsReco.push_back(histReco);

//			TH1D* histBkgReco = (TH1D*)(FileSample[WhichSample]->Get("BkgRecoTrue"+PlotNames[WhichPlot]));
			TH1D* histBkgReco = (TH1D*)(FileSample[WhichSample]->Get("NonCC1pReco"+PlotNames[WhichPlot]));
			CurrentPlotsBkgReco.push_back(histBkgReco);

			TH1D* histCC1pReco = (TH1D*)(FileSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]));
			CurrentPlotsCC1pReco.push_back(histCC1pReco);

			TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
			CurrentPlotsTrue.push_back(histTrue);

			TH1D* histTEfficiency = (TH1D*)(FileEfficiences->Get("CC1pReco"+PlotNames[WhichPlot]));
			CurrentPlotsTEfficiency.push_back(histTEfficiency);
		
		}

		PlotsReco.push_back(CurrentPlotsReco);		
		PlotsTrue.push_back(CurrentPlotsTrue);		
		PlotsBkgReco.push_back(CurrentPlotsBkgReco);
		PlotsCC1pReco.push_back(CurrentPlotsCC1pReco);
		PlotsTEfficiency.push_back(CurrentPlotsTEfficiency);

	}

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Loop over the plots

	for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){
	
		PlotCanvas.push_back(new TCanvas(PlotNames[WhichPlot],PlotNames[WhichPlot],205,34,1024,768));
		PlotCanvas[WhichPlot]->cd();

		TPad *midPad = new TPad("midPad", "", 0.005, 0.3  , 0.995, 0.995);
		TPad *botPad = new TPad("botPad", "", 0.005, 0.005, 0.995, 0.3);
		midPad->SetBottomMargin(0.105);
		midPad->SetTopMargin(0.1);
		botPad->SetTopMargin(0.1);
		botPad->SetBottomMargin(0.05);
		botPad->SetGridy();
		midPad->Draw();
		botPad->Draw();

		TLegend* leg = new TLegend(0.1,0.92,0.65,1.);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.08);
		leg->SetTextFont(FontStyle);
		leg->SetNColumns(3);

		double NOverlay = (double)(PlotsCC1pReco[0][WhichPlot]->GetEntries());
		double NGenie = (double)(PlotsTrue[4][WhichPlot]->GetEntries());

		int NBinsX = PlotsCC1pReco[0][WhichPlot]->GetNbinsX();

		// -------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Apply the relevant weights

		for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

			double BinWidth = PlotsCC1pReco[0][WhichPlot]->GetBinWidth(WhichXBin+1);

			// -------------------------------------------------------------------------------------------------------------------------------------------------------------------
		
			// Effective Efficiencies

			double EffectiveEfficiencyXBin = PlotsTEfficiency[0][WhichPlot]->GetBinContent(WhichXBin+1);
//			double EffectiveEfficiencyXBin = 1.;
			double EffectiveEfficiencyXBinError = PlotsTEfficiency[0][WhichPlot]->GetBinError(WhichXBin+1);
	
			// -------------------------------------------------------------------------------------------------------------------------------------------------------------------

			// Run 1 Data

			double CurrentDataEntry = PlotsReco[1][WhichPlot]->GetBinContent(WhichXBin+1);
			double CurrentDataError = PlotsReco[1][WhichPlot]->GetBinError(WhichXBin+1);
			double DataScaledEntry = 0., DataScaledError = 0.;
			if (EffectiveEfficiencyXBin != 0 ) {  
				DataScaledEntry = CurrentDataEntry / EffectiveEfficiencyXBin * (Units / (Flux * NTargets * BinWidth)); 
				DataScaledError = sqrt( 
							TMath::Power(CurrentDataError / EffectiveEfficiencyXBin,2)  +
							TMath::Power(CurrentDataEntry * EffectiveEfficiencyXBinError / (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2)
						      ) * (Units / (Flux * NTargets * BinWidth));
			}

			PlotsReco[1][WhichPlot]->SetBinContent(WhichXBin+1,DataScaledEntry);
			PlotsReco[1][WhichPlot]->SetBinError(WhichXBin+1,DataScaledError);

			// -------------------------------------------------------------------------------------------------------------------------------------------------------------------

			// ExtBNB

			double CurrentExtBNBEntry = PlotsReco[2][WhichPlot]->GetBinContent(WhichXBin+1);
			double CurrentExtBNBError = PlotsReco[2][WhichPlot]->GetBinError(WhichXBin+1);
			double ExtBNBScaledEntry = 0., ExtBNBScaledError = 0.;
			if (EffectiveEfficiencyXBin != 0 ) {  
				ExtBNBScaledEntry = CurrentExtBNBEntry / EffectiveEfficiencyXBin * (Units / (Flux * NTargets * BinWidth))
							; 
				ExtBNBScaledError = sqrt( 
							TMath::Power(CurrentExtBNBError / EffectiveEfficiencyXBin,2)  +
							TMath::Power(CurrentExtBNBEntry * EffectiveEfficiencyXBinError / (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2)
						      ) * (Units / (Flux * NTargets * BinWidth))
							;
			}

			PlotsReco[2][WhichPlot]->SetBinContent(WhichXBin+1,ExtBNBScaledEntry);
			PlotsReco[2][WhichPlot]->SetBinError(WhichXBin+1,ExtBNBScaledError);

			// -------------------------------------------------------------------------------------------------------------------------------------------------------------------

			// Dirt

			double CurrentDirtEntry = PlotsReco[3][WhichPlot]->GetBinContent(WhichXBin+1);
			double CurrentDirtError = PlotsReco[3][WhichPlot]->GetBinError(WhichXBin+1);
			double DirtScaledEntry = 0., DirtScaledError = 0.;
			if (EffectiveEfficiencyXBin != 0 ) {  
				DirtScaledEntry = CurrentDirtEntry / EffectiveEfficiencyXBin * (Units / (Flux * NTargets * BinWidth))
							; 
				DirtScaledError = sqrt( 
							TMath::Power(CurrentDirtError / EffectiveEfficiencyXBin,2)  +
							TMath::Power(CurrentDirtEntry * EffectiveEfficiencyXBinError / (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2)
						      ) * (Units / (Flux * NTargets * BinWidth))
							;
			}

			PlotsReco[3][WhichPlot]->SetBinContent(WhichXBin+1,DirtScaledEntry);
			PlotsReco[3][WhichPlot]->SetBinError(WhichXBin+1,DirtScaledError);

			// -------------------------------------------------------------------------------------------------------------------------------------------------------------------

			// Overlay

			double CurrentOverlayEntry = PlotsCC1pReco[0][WhichPlot]->GetBinContent(WhichXBin+1);
			double CurrentOverlayError = PlotsCC1pReco[0][WhichPlot]->GetBinError(WhichXBin+1);
			double OverlayScaledEntry = 0., OverlayScaledError = 0.;
			if (EffectiveEfficiencyXBin != 0 ) {
				OverlayScaledEntry = CurrentOverlayEntry / EffectiveEfficiencyXBin * (Units/(Flux * NTargets * BinWidth)) 
							; 
				OverlayScaledError = sqrt( 
							TMath::Power(CurrentOverlayError / EffectiveEfficiencyXBin,2)  +
							TMath::Power(CurrentOverlayEntry * EffectiveEfficiencyXBinError / (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2)
						      ) * (Units / (Flux * NTargets * BinWidth))
							;


			}

			PlotsCC1pReco[0][WhichPlot]->SetBinContent(WhichXBin+1,OverlayScaledEntry);
			PlotsCC1pReco[0][WhichPlot]->SetBinError(WhichXBin+1,OverlayScaledError);

			// -------------------------------------------------------------------------------------------------------------------------------------------------------------------

			// Bkg

			double CurrentBkgEntry = PlotsBkgReco[0][WhichPlot]->GetBinContent(WhichXBin+1);
			double CurrentBkgError = PlotsBkgReco[0][WhichPlot]->GetBinError(WhichXBin+1);
			double BkgScaledEntry = 0., BkgScaledError = 0.;
			if (EffectiveEfficiencyXBin != 0 ) {
				BkgScaledEntry = CurrentBkgEntry / EffectiveEfficiencyXBin * (Units / (Flux * NTargets * BinWidth)) 
							; 
				BkgScaledError = sqrt( 
							TMath::Power(CurrentBkgError / EffectiveEfficiencyXBin,2)  +
							TMath::Power(CurrentBkgEntry * EffectiveEfficiencyXBinError / (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2)
						      ) * (Units / (Flux * NTargets * BinWidth))
							;


			}

			PlotsBkgReco[0][WhichPlot]->SetBinContent(WhichXBin+1,BkgScaledEntry);
			PlotsBkgReco[0][WhichPlot]->SetBinError(WhichXBin+1,BkgScaledError);

			// -------------------------------------------------------------------------------------------------------------------------------------------------------------------

			// Genie

			double CurrentGenieEntry = PlotsTrue[4][WhichPlot]->GetBinContent(WhichXBin+1);
			double CurrentGenieError = PlotsTrue[4][WhichPlot]->GetBinError(WhichXBin+1);
			double GenieScaledEntry = 0., GenieScaledError = 0.;

			GenieScaledEntry = FluxIntegratedXSection * CurrentGenieEntry / BinWidth / NGenieEvents ; 
			GenieScaledError = FluxIntegratedXSection * CurrentGenieError / BinWidth / NGenieEvents ; 

			PlotsTrue[4][WhichPlot]->SetBinContent(WhichXBin+1,GenieScaledEntry);
			PlotsTrue[4][WhichPlot]->SetBinError(WhichXBin+1,GenieScaledError);

			// -------------------------------------------------------------------------------------------------------------------------------------------------------------------

			// Samples for systematics

		} // End of the loop over the bins

		// -------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Plotting the xsections

		//Overlay

		MakeMyPlotPretty(PlotsCC1pReco[0][WhichPlot]);
		PlotsCC1pReco[0][WhichPlot]->SetLineColor(kBlue);
		PlotsCC1pReco[0][WhichPlot]->SetFillColor(kBlue);
		PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitleOffset(0.63);
		PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitle(PlotXAxis[WhichPlot]);
		midPad->cd();

		//Run 1 Data

		MakeMyPlotPretty(PlotsReco[1][WhichPlot]);
		PlotsReco[1][WhichPlot]->SetLineColor(kOrange+7);
		PlotsReco[1][WhichPlot]->SetMarkerColor(kOrange+7);
		PlotsReco[1][WhichPlot]->SetMarkerStyle(20);
		PlotsReco[1][WhichPlot]->SetMarkerSize(1.5);
		// Subtract ExtBNB
		PlotsReco[1][WhichPlot]->Add(PlotsReco[2][WhichPlot],-1); 
		// Subtract Dirt
		PlotsReco[1][WhichPlot]->Add(PlotsReco[3][WhichPlot],-1); 
		// Subtract MC Bkg Subtraction
		if (Subtract == "") { PlotsReco[1][WhichPlot]->Add(PlotsBkgReco[0][WhichPlot],-1); }

		//Genie

		MakeMyPlotPretty(PlotsTrue[4][WhichPlot]);
		PlotsTrue[4][WhichPlot]->SetLineColor(kBlack);
		PlotsTrue[4][WhichPlot]->SetFillColor(kBlack);

		double max = TMath::Max(PlotsCC1pReco[0][WhichPlot]->GetMaximum(),PlotsReco[1][WhichPlot]->GetMaximum());
		PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetRangeUser(0.,1.8 * max);
		PlotsCC1pReco[0][WhichPlot]->Draw("e2");
		PlotsTrue[4][WhichPlot]->Draw("e same");
		PlotsReco[1][WhichPlot]->Draw("e same");

		// -------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Integrated cross-sections & chi2

		double IntegratedGenieXSection = 0., IntegratedGenieXSectionErrorSquared = 0., IntegratedGenieXSectionError = 0.;
		double IntegratedOverlayXSection = 0., IntegratedOverlayXSectionErrorSquared = 0., IntegratedOverlayXSectionError = 0.;
		double IntegratedDataXSection = 0., IntegratedDataXSectionErrorSquared = 0., IntegratedDataXSectionError = 0.;

		double chi2 = 0, num = 0, den = 0;

		for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

			double BinWidth = PlotsCC1pReco[0][WhichPlot]->GetBinWidth(WhichXBin+1);

			double GenieBinEntry = PlotsTrue[4][WhichPlot]->GetBinContent(WhichXBin+1);
			double GenieBinError = PlotsTrue[4][WhichPlot]->GetBinError(WhichXBin+1);
			IntegratedGenieXSection += GenieBinEntry * BinWidth;
			IntegratedGenieXSectionErrorSquared += TMath::Power(GenieBinError * BinWidth,2.);

			double OverlayBinEntry = PlotsCC1pReco[0][WhichPlot]->GetBinContent(WhichXBin+1);
			double OverlayBinError = PlotsCC1pReco[0][WhichPlot]->GetBinError(WhichXBin+1);
			IntegratedOverlayXSection += OverlayBinEntry * BinWidth;
			IntegratedOverlayXSectionErrorSquared += TMath::Power(OverlayBinError * BinWidth,2.);

			double DataBinEntry = PlotsReco[1][WhichPlot]->GetBinContent(WhichXBin+1);
			double DataBinError = PlotsReco[1][WhichPlot]->GetBinError(WhichXBin+1);
			IntegratedDataXSection += DataBinEntry * BinWidth;
			IntegratedDataXSectionErrorSquared += TMath::Power(DataBinError * BinWidth,2.);

			num = TMath::Power(DataBinEntry - OverlayBinEntry,2.);
			den = TMath::Power(DataBinError,2.) + TMath::Power(OverlayBinError,2.);
			if (num != 0 && den != 0) {chi2 += (num / den); }

		}

		int accuracy = 100;

		IntegratedGenieXSection = roundf(IntegratedGenieXSection * accuracy) / accuracy;
		IntegratedGenieXSectionError = sqrt(IntegratedGenieXSectionErrorSquared);
		IntegratedGenieXSectionError = roundf(IntegratedGenieXSectionError * accuracy) / accuracy;

		IntegratedOverlayXSection = roundf(IntegratedOverlayXSection * accuracy) / accuracy;
		IntegratedOverlayXSectionError = sqrt(IntegratedOverlayXSectionErrorSquared);
		IntegratedOverlayXSectionError = roundf(IntegratedOverlayXSectionError * accuracy) / accuracy;

		IntegratedDataXSection = roundf(IntegratedDataXSection * accuracy) / accuracy;
		IntegratedDataXSectionError = sqrt(IntegratedDataXSectionErrorSquared);
		IntegratedDataXSectionError = roundf(IntegratedDataXSectionError * accuracy) / accuracy;

		TLatex latexSigma;
		latexSigma.SetTextFont(FontStyle);
		latexSigma.SetTextSize(0.07);
		TString LabelData = "#splitline{#splitline{#color[807]{#sigma_{Data} = (" +ToString(IntegratedDataXSection)+" #pm "+
				    ToString(IntegratedDataXSectionError)+") #upoint 10^{-38} cm^{2}}}{#color[600]{#sigma_{MC} = (" +
				    ToString(IntegratedOverlayXSection)+" #pm "+ToString(IntegratedOverlayXSectionError)+") #upoint 10^{-38} cm^{2}}}}{#color[1]{#sigma_{Genie} = (" +
				    ToString(IntegratedGenieXSection)+" #pm "+ToString(IntegratedGenieXSectionError)+") #upoint 10^{-38} cm^{2}}}";
		latexSigma.DrawLatexNDC(0.3,0.67, LabelData);

		// -------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Legend & POT Normalization

		leg->AddEntry(PlotsCC1pReco[0][WhichPlot],"CC1p MC","f");
		leg->AddEntry(PlotsTrue[4][WhichPlot],"Genie","l");
		leg->AddEntry(PlotsReco[1][WhichPlot],"Data","lep");
		leg->Draw();	

		TLatex latex;
		latex.SetTextFont(FontStyle);
		latex.SetTextSize(0.07);
		TString Label = ToString(tor860_wcut)+" POT";
		latex.DrawLatexNDC(0.67,0.94, Label);

		// -------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Residual plot

		botPad->cd();
		TH1D* PlotsRecoClone = (TH1D*)(PlotsReco[1][WhichPlot]->Clone());
		PlotsRecoClone->Add(PlotsCC1pReco[0][WhichPlot],-1);
		PlotsRecoClone->Divide(PlotsReco[1][WhichPlot]);

		PlotsRecoClone->SetLineColor(kBlack);
		PlotsRecoClone->SetMarkerColor(kBlack);
		PlotsRecoClone->GetXaxis()->SetTitle();
		PlotsRecoClone->GetXaxis()->SetLabelSize(0);
		PlotsRecoClone->GetYaxis()->SetTitle("#frac{Data-MC}{Data}");
		PlotsRecoClone->GetYaxis()->SetLabelSize(0.1);
		PlotsRecoClone->GetYaxis()->SetRangeUser(-1,1);
		PlotsRecoClone->GetYaxis()->SetTitleSize(0.15);
		PlotsRecoClone->GetYaxis()->SetTitleOffset(0.25);
		PlotsRecoClone->GetYaxis()->SetTitleFont(132);

		PlotsRecoClone->Draw("e same");

		TLatex latexChi2;
		latexChi2.SetTextFont(FontStyle);
		latexChi2.SetTextSize(0.15);
		TString LabelChi2 = "#chi^{2} = " + ToString(chi2) + " / " + ToString(NBinsX-1);
		latexChi2.DrawLatexNDC(0.4,0.3, LabelChi2);

		double RatioMin = PlotsRecoClone->GetXaxis()->GetXmin();
		double RatioMax = PlotsRecoClone->GetXaxis()->GetXmax();
		double YRatioCoord = 0.25;
		TLine* RatioLine = new TLine(RatioMin,YRatioCoord,RatioMax,YRatioCoord);
		RatioLine->SetLineWidth(4);
		RatioLine->SetLineColor(kPink+8);
		RatioLine->SetLineStyle(4);
		RatioLine->Draw("same");

		PlotCanvas[WhichPlot]->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/"+NameOfSamples[0]+"/XSections_"+PlotNames[WhichPlot]+"_"+WhichRun+"_"+UBCodeVersion+Subtract+".pdf");
		PlotCanvas[WhichPlot]->SaveAs("./myPlots/eps/"+UBCodeVersion+"/"+NameOfSamples[0]+"/XSections_"+PlotNames[WhichPlot]+"_"+WhichRun+"_"+UBCodeVersion+Subtract+".eps");
		//delete PlotCanvas[WhichPlot];

	} // End of the loop over the plots

	FileEfficiences->Close();

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Store the extracted xsections

	TString NameExtractedXSec = PathToExtractedXSec+UBCodeVersion+"/ExtractedXSec_"+NameOfSamples[0]+"_"+WhichRun+"_"+UBCodeVersion+".root";
	TFile* ExtractedXSec = TFile::Open(NameExtractedXSec,"recreate");

	for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

		PlotsReco[1][WhichPlot]->Write();
		PlotsCC1pReco[0][WhichPlot]->Write();
		PlotsTrue[4][WhichPlot]->Write();

	}

	std::cout << std::endl << "File " << NameExtractedXSec << " created" << std::endl << std::endl;

} // End of the program 
