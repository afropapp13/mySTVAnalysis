#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
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
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/myFunctions.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/MakeMyPlotPretty.cpp"

#include "./Constants.h"

using namespace std;
using namespace Constants;

void XSection_Extraction(TString OverlaySample) {

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;

	TString Subtract = "";
//	TString Subtract = "_BUnsubtracted";

	int DecimalAccuracy = 2;

	TString PathToFiles = "../myEvents/OutputFiles/";

	TString PathToExtractedXSec = "myXSec/";

	// -------------------------------------------------------------------------------------

	int NEventsPassingSelectionCuts = 0;
	TString CutExtension = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
	VectorCuts.push_back("");
	VectorCuts.push_back("_NuScore");
	VectorCuts.push_back("_ThreePlaneLogChi2");
	VectorCuts.push_back("_Collinearity");

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
	PlotNames.push_back("ECalPlot");
	PlotNames.push_back("EQEPlot"); 
	PlotNames.push_back("Q2Plot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();

	vector<TString> LabelsOfSamples;
	vector<TString> NameOfSamples;

	NameOfSamples.push_back("Overlay9"); 
	NameOfSamples.push_back("BeamOn9"); 
	NameOfSamples.push_back("ExtBNB9"); 
	NameOfSamples.push_back("OverlayDirt9"); 
	NameOfSamples.push_back("Genie");

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		vector<TCanvas*> PlotCanvas; PlotCanvas.clear();

		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
		vector<vector<TH1D*> > PlotsBkgReco; PlotsBkgReco.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();
		vector<vector<TH1D*> > PlotsTEfficiency; PlotsTEfficiency.clear();

		// -----------------------------------------------------------------------------------------------------------------------------------------


		TFile* FileEfficiences = new TFile("myEfficiencies/"+UBCodeVersion+"/FileEfficiences_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root","readonly");

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			if (NameOfSamples[WhichSample] == "BeamOn9" || NameOfSamples[WhichSample] == "ExtBNB9" || NameOfSamples[WhichSample] == "OverlayDirt9") 
				{ FileSample.push_back(TFile::Open(PathToFiles+"/"+UBCodeVersion+"/"+CutExtension+"/STVStudies_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+CutExtension+".root")); }
			if (NameOfSamples[WhichSample] == "Overlay9") 
				{ FileSample.push_back(TFile::Open(PathToFiles+"/"+UBCodeVersion+"/"+CutExtension+"/STVStudies_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+OverlaySample+CutExtension+".root")); }
			if (NameOfSamples[WhichSample] == "Genie") { FileSample.push_back(TFile::Open("myFiles/"+UBCodeVersion+"/CCQEAnalysis_"+NameOfSamples[WhichSample]+"_"+UBCodeVersion+".root")); }

			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();
			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();
			vector<TH1D*> CurrentPlotsBkgReco; CurrentPlotsBkgReco.clear();
			vector<TH1D*> CurrentPlotsCC1pReco; CurrentPlotsCC1pReco.clear();
			vector<TH1D*> CurrentPlotsTEfficiency; CurrentPlotsTEfficiency.clear();

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
				CurrentPlotsReco.push_back(histReco);

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
	
			TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+"_"+Runs[WhichRun]+OverlaySample,PlotNames[WhichPlot]+"_"+Runs[WhichRun]+OverlaySample,205,34,1024,768);
			PlotCanvas->cd();

			TPad *midPad = new TPad("midPad", "", 0.005,0.0, 0.995, 0.995);
			midPad->SetBottomMargin(0.14);
			midPad->SetTopMargin(0.12);
			midPad->SetLeftMargin(0.16);
			midPad->Draw();

			TLegend* leg = new TLegend(0.14,0.89,0.8,0.99);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.06);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(2);

			int NBinsX = PlotsCC1pReco[0][WhichPlot]->GetNbinsX();

			// -----------------------------------------------------------------------------------------------------------------------------------------

			// Apply the relevant weights

			for (int WhichXBin = 0; WhichXBin < NBinsX; WhichXBin++) {

				double BinWidth = PlotsCC1pReco[0][WhichPlot]->GetBinWidth(WhichXBin+1);
// Put it back !
				double ScalingFactor = (Units / (Flux * NTargets * BinWidth)) * (A/Z);
//				double ScalingFactor = (Units / (Flux * NTargets * BinWidth));

				// ----------------------------------------------------------------------------------------------------------------------------------
		
				// Effective Efficiencies
				double EffectiveEfficiencyXBin = PlotsTEfficiency[0][WhichPlot]->GetBinContent(WhichXBin+1);
				double EffectiveEfficiencyXBinError = PlotsTEfficiency[0][WhichPlot]->GetBinError(WhichXBin+1);

//test	
//				double EffectiveEfficiencyXBin = 1.;
//				double EffectiveEfficiencyXBinError = 0.;

				// -----------------------------------------------------------------------------------------------------------------------------------

				// BeamOn

				double CurrentDataEntry = PlotsReco[1][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentDataError = PlotsReco[1][WhichPlot]->GetBinError(WhichXBin+1);
				double DataScaledEntry = 0., DataScaledError = 0.;
				if (EffectiveEfficiencyXBin != 0 ) {  
					DataScaledEntry = CurrentDataEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					DataScaledError = sqrt( 
								TMath::Power(CurrentDataError / EffectiveEfficiencyXBin,2)  +
								TMath::Power(CurrentDataEntry * EffectiveEfficiencyXBinError / (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2)
							      ) * ScalingFactor;
				}

				PlotsReco[1][WhichPlot]->SetBinContent(WhichXBin+1,DataScaledEntry);
				PlotsReco[1][WhichPlot]->SetBinError(WhichXBin+1,DataScaledError);

				// ------------------------------------------------------------------------------------------------------------------------------------------

				// ExtBNB

				double CurrentExtBNBEntry = PlotsReco[2][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentExtBNBError = PlotsReco[2][WhichPlot]->GetBinError(WhichXBin+1);
				double ExtBNBScaledEntry = 0., ExtBNBScaledError = 0.;
				if (EffectiveEfficiencyXBin != 0 ) {  
					ExtBNBScaledEntry = CurrentExtBNBEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					ExtBNBScaledError = sqrt( 
								TMath::Power(CurrentExtBNBError / EffectiveEfficiencyXBin,2)  +
								TMath::Power(CurrentExtBNBEntry * EffectiveEfficiencyXBinError / (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2)
							      ) * ScalingFactor
								;
				}

				PlotsReco[2][WhichPlot]->SetBinContent(WhichXBin+1,ExtBNBScaledEntry);
				PlotsReco[2][WhichPlot]->SetBinError(WhichXBin+1,ExtBNBScaledError);

				// ----------------------------------------------------------------------------------------------------------------------------------------

				// Dirt

				double CurrentDirtEntry = PlotsReco[3][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentDirtError = PlotsReco[3][WhichPlot]->GetBinError(WhichXBin+1);
				double DirtScaledEntry = 0., DirtScaledError = 0.;
				if (EffectiveEfficiencyXBin != 0 ) {  
					DirtScaledEntry = CurrentDirtEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					DirtScaledError = sqrt( 
								TMath::Power(CurrentDirtError / EffectiveEfficiencyXBin,2)  +
								TMath::Power(CurrentDirtEntry * EffectiveEfficiencyXBinError / (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2)
							      ) * ScalingFactor
								;
				}

				PlotsReco[3][WhichPlot]->SetBinContent(WhichXBin+1,DirtScaledEntry);
				PlotsReco[3][WhichPlot]->SetBinError(WhichXBin+1,DirtScaledError);

				// --------------------------------------------------------------------------------------------------------------------------------------------

				// Overlay

				double CurrentOverlayEntry = PlotsCC1pReco[0][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentOverlayError = PlotsCC1pReco[0][WhichPlot]->GetBinError(WhichXBin+1);
				double OverlayScaledEntry = 0., OverlayScaledError = 0.;
				if (EffectiveEfficiencyXBin != 0 ) {
					OverlayScaledEntry = CurrentOverlayEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					OverlayScaledError = sqrt( 
								TMath::Power(CurrentOverlayError / EffectiveEfficiencyXBin,2)  +
								TMath::Power(CurrentOverlayEntry * EffectiveEfficiencyXBinError / (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2)
							      ) * ScalingFactor
								;


				}

				PlotsCC1pReco[0][WhichPlot]->SetBinContent(WhichXBin+1,OverlayScaledEntry);
				PlotsCC1pReco[0][WhichPlot]->SetBinError(WhichXBin+1,OverlayScaledError);

				// ---------------------------------------------------------------------------------------------------------------------------------------

				// Bkg

				double CurrentBkgEntry = PlotsBkgReco[0][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentBkgError = PlotsBkgReco[0][WhichPlot]->GetBinError(WhichXBin+1);
				double BkgScaledEntry = 0., BkgScaledError = 0.;
				if (EffectiveEfficiencyXBin != 0 ) {
					BkgScaledEntry = CurrentBkgEntry / EffectiveEfficiencyXBin * ScalingFactor; 
					BkgScaledError = sqrt( 
								TMath::Power(CurrentBkgError / EffectiveEfficiencyXBin,2)  +
								TMath::Power(CurrentBkgEntry * EffectiveEfficiencyXBinError / (EffectiveEfficiencyXBin * EffectiveEfficiencyXBin),2)
							      ) * ScalingFactor
								;


				}

				PlotsBkgReco[0][WhichPlot]->SetBinContent(WhichXBin+1,BkgScaledEntry);
				PlotsBkgReco[0][WhichPlot]->SetBinError(WhichXBin+1,BkgScaledError);

				// ----------------------------------------------------------------------------------------------------------------------------------------

				// Genie

				double CurrentGenieEntry = PlotsTrue[4][WhichPlot]->GetBinContent(WhichXBin+1);
				double CurrentGenieError = PlotsTrue[4][WhichPlot]->GetBinError(WhichXBin+1);
				double GenieScaledEntry = 0., GenieScaledError = 0.;

				GenieScaledEntry = FluxIntegratedXSection * CurrentGenieEntry / BinWidth / NGenieEvents ; 
				GenieScaledError = FluxIntegratedXSection * CurrentGenieError / BinWidth / NGenieEvents ; 

				PlotsTrue[4][WhichPlot]->SetBinContent(WhichXBin+1,GenieScaledEntry);
				PlotsTrue[4][WhichPlot]->SetBinError(WhichXBin+1,GenieScaledError);

				// -----------------------------------------------------------------------------------------------------------------------------------------

				// Samples for systematics

			} // End of the loop over the bins

			// ---------------------------------------------------------------------------------------------------------------------------------------------------

			// Plotting the xsections

			// Overlay

			MakeMyPlotPretty(PlotsCC1pReco[0][WhichPlot]);
			PlotsCC1pReco[0][WhichPlot]->SetLineColor(OverlayColor);
			PlotsCC1pReco[0][WhichPlot]->SetFillColor(OverlayColor);

			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetTitleOffset(1.05);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetTitleSize(0.06);
			PlotsCC1pReco[0][WhichPlot]->GetXaxis()->SetLabelSize(0.06);

			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitleOffset(1.2);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitle(PlotXAxis[WhichPlot]);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetTitleSize(0.06);
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetLabelSize(0.06);

			midPad->cd();

			// BeamOn

			MakeMyPlotPretty(PlotsReco[1][WhichPlot]);
			PlotsReco[1][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsReco[1][WhichPlot]->SetMarkerColor(BeamOnColor);
			PlotsReco[1][WhichPlot]->SetMarkerStyle(20);
			PlotsReco[1][WhichPlot]->SetMarkerSize(1.5);
			// Subtract ExtBNB
			PlotsReco[1][WhichPlot]->Add(PlotsReco[2][WhichPlot],-1); 
			// Subtract Dirt
			PlotsReco[1][WhichPlot]->Add(PlotsReco[3][WhichPlot],-1); 
			// Subtract MC Bkg Subtraction
			if (Subtract == "") { PlotsReco[1][WhichPlot]->Add(PlotsBkgReco[0][WhichPlot],-1); }

			// GENIE

			MakeMyPlotPretty(PlotsTrue[4][WhichPlot]);
			PlotsTrue[4][WhichPlot]->SetLineColor(GenieColor);
			PlotsTrue[4][WhichPlot]->SetFillColor(GenieColor);

			// Max on plot so that we can include the integrated xsecs
			double max = TMath::Max(PlotsCC1pReco[0][WhichPlot]->GetMaximum(),PlotsReco[1][WhichPlot]->GetMaximum());
			PlotsCC1pReco[0][WhichPlot]->GetYaxis()->SetRangeUser(0.,1.8 * max);

			// Plot MC
			PlotsCC1pReco[0][WhichPlot]->Draw("e2");

			// Plot GENIE
			PlotsTrue[4][WhichPlot]->Draw("e same");

			// Plot BeamOn
			PlotsReco[1][WhichPlot]->Draw("ex0 same");

			// --------------------------------------------------------------------------------------------------------------------------------------------------

			// Integrated cross-sections & chi2

			double IntegratedGenieXSection = round(IntegratedXSec(PlotsTrue[4][WhichPlot]),DecimalAccuracy);
			double IntegratedGenieXSectionError = round(IntegratedXSecError(PlotsTrue[4][WhichPlot]),DecimalAccuracy);

			double IntegratedOverlayXSection = round(IntegratedXSec(PlotsCC1pReco[0][WhichPlot]),DecimalAccuracy);
			double IntegratedOverlayXSectionError = round(IntegratedXSecError(PlotsCC1pReco[0][WhichPlot]),DecimalAccuracy);

			double IntegratedDataXSection = round(IntegratedXSec(PlotsReco[1][WhichPlot]),DecimalAccuracy);
			double IntegratedDataXSectionError = round(IntegratedXSecError(PlotsReco[1][WhichPlot]),DecimalAccuracy);

			double chi2 = Chi2(PlotsReco[1][WhichPlot],PlotsCC1pReco[0][WhichPlot]);

			TString LabelData = "#splitline{#splitline{#color["+ToString(BeamOnColor)+"]{#sigma_{Data} = (" +ToString(IntegratedDataXSection)+" #pm "+
					    ToString(IntegratedDataXSectionError)+") #upoint 10^{-38} cm^{2}}}{#color["+ToString(OverlayColor)+"]{#sigma_{MC} = (" +
					    ToString(IntegratedOverlayXSection)+" #pm "+
					    ToString(IntegratedOverlayXSectionError)+") #upoint 10^{-38} cm^{2}}}}{#color["+ToString(GenieColor)+"]{#sigma_{GENIE} = (" +
					    ToString(IntegratedGenieXSection)+" #pm "+ToString(IntegratedGenieXSectionError)+") #upoint 10^{-38} cm^{2}}}";

			TLatex latexSigma;
			latexSigma.SetTextFont(FontStyle);
			latexSigma.SetTextSize(0.06);
			latexSigma.DrawLatexNDC(0.27,0.67, LabelData);

			// ---------------------------------------------------------------------------------------------------------------------------------------------

			// Legend & POT Normalization

			double tor860_wcut = -99.;
			if (Runs[WhichRun] == "Run1") { tor860_wcut = tor860_wcut_Run1; }

			TString Label = ToString(tor860_wcut)+" POT";
//			latex.DrawLatexNDC(0.47,0.9, Label);

			TLegendEntry* lMC = leg->AddEntry(PlotsCC1pReco[0][WhichPlot],"MC","f");
			lMC->SetTextColor(OverlayColor);

			TLegendEntry* lGenie = leg->AddEntry(PlotsTrue[4][WhichPlot],"GENIE","l");
			lGenie->SetTextColor(GenieColor);

			leg->AddEntry(PlotsReco[1][WhichPlot],"MicroBooNE Data " + Runs[WhichRun] + " " + Label,"ep");
			leg->Draw();	

			// ----------------------------------------------------------------------------------------------------------------------------------------------------

			if (OverlaySample == "") {
				PlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/"+NameOfSamples[0]+"/XSections_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+Subtract+".pdf");
			}

			if (OverlaySample != "") { delete PlotCanvas; }

		} // End of the loop over the plots

		FileEfficiences->Close();

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Store the extracted xsections

		TString NameExtractedXSec = PathToExtractedXSec+UBCodeVersion+"/ExtractedXSec_"+NameOfSamples[0]+"_"+Runs[WhichRun]+OverlaySample+"_"+UBCodeVersion+".root";
		TFile* ExtractedXSec = TFile::Open(NameExtractedXSec,"recreate");

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

			PlotsReco[1][WhichPlot]->Write();
			PlotsCC1pReco[0][WhichPlot]->Write();
			PlotsTrue[4][WhichPlot]->Write();

		}

		ExtractedXSec->Close();

		std::cout << std::endl << "File " << NameExtractedXSec << " created" << std::endl << std::endl;

	} // End of the loop over the runs	

} // End of the program 
