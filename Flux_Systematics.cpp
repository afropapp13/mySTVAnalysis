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
#include <TGaxis.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <numeric>
#include <functional>

#include  "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/MakeMyPlotPretty.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/myFunctions.cpp"

#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

void Flux_Systematics() {

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;
	SetOffsetAndSize();
	TGaxis::SetMaxDigits(3);
	TGaxis::SetExponentOffset(-0.1, 1., "y");	
	
	double TextSize = 0.07;
	
	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(TextSize,"t");	
	
	// ---------------------------------------------------------------------------------------------------------------------------------------	

	TString PathToFiles = "myXSec/";

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

	// ---------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");
//	Runs.push_back("Run2");
//	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");				

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;
	
	// ---------------------------------------------------------------------------------------------------------------------------------------
	
	std::vector<int> NUniverses; NUniverses.clear();	
		
	std::vector<TString> EventWeightLabels; EventWeightLabels.clear();
	
//	int LocalNUniverses = 100;		
	int LocalNUniverses = 80;
	
	EventWeightLabels.push_back("horncurrent_FluxUnisim"); NUniverses.push_back(LocalNUniverses);
	EventWeightLabels.push_back("kminus_PrimaryHadronNormalization"); NUniverses.push_back(LocalNUniverses);
	EventWeightLabels.push_back("kplus_PrimaryHadronFeynmanScaling"); NUniverses.push_back(LocalNUniverses);
	EventWeightLabels.push_back("kzero_PrimaryHadronSanfordWang"); NUniverses.push_back(LocalNUniverses);
	EventWeightLabels.push_back("nucleoninexsec_FluxUnisim"); NUniverses.push_back(LocalNUniverses);
	EventWeightLabels.push_back("nucleonqexsec_FluxUnisim"); NUniverses.push_back(LocalNUniverses);
	EventWeightLabels.push_back("nucleontotxsec_FluxUnisim"); NUniverses.push_back(LocalNUniverses);
	EventWeightLabels.push_back("piminus_PrimaryHadronSWCentralSplineVariation"); NUniverses.push_back(LocalNUniverses);
	EventWeightLabels.push_back("pioninexsec_FluxUnisim"); NUniverses.push_back(LocalNUniverses);
	EventWeightLabels.push_back("pionqexsec_FluxUnisim"); NUniverses.push_back(LocalNUniverses);
	EventWeightLabels.push_back("piontotxsec_FluxUnisim"); NUniverses.push_back(LocalNUniverses);
	EventWeightLabels.push_back("piplus_PrimaryHadronSWCentralSplineVariation"); NUniverses.push_back(LocalNUniverses);
	EventWeightLabels.push_back("expskin_FluxUnisim"); NUniverses.push_back(10);
		
	int NEventWeightLabels = EventWeightLabels.size();

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	cout << endl;

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		TString SystFileName = "mySystematics/"+UBCodeVersion+"/Flux_Systematics_"+Runs[WhichRun]+".root";
		TFile* SystFile = new TFile(SystFileName,"recreate");

		vector<vector<TH1D*> > PlotsReco;

		vector<TString> NameOfSamples;
		std::vector<int> Colors;
		std::vector<int> Markers;
		
		std::vector< std::vector<double> > TotalSystXsecs;
		TotalSystXsecs.resize(N1DPlots);
		
		for (int WhichEventWeightLabel = 0; WhichEventWeightLabel < NEventWeightLabels; WhichEventWeightLabel ++) {
		
			PlotsReco.clear();
			NameOfSamples.clear();
			Colors.clear();
			Markers.clear();
			
			NameOfSamples.push_back(""); // Reference plot
			Colors.push_back(kBlack); 
			Markers.push_back(20);			
		
			for (int i = 0; i < NUniverses[WhichEventWeightLabel]; i ++) {
			
				NameOfSamples.push_back("_"+EventWeightLabels[WhichEventWeightLabel]+"_"+TString(std::to_string(i))); 
				Colors.push_back(kGreen+2); 
				Markers.push_back(22);
			}

			const int NSamples = NameOfSamples.size();
			vector<TFile*> FileSample; FileSample.clear();

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {


				FileSample.push_back(TFile::Open(PathToFiles+UBCodeVersion+"/ExtractedXSec_Overlay9_"+\
						      Runs[WhichRun]+NameOfSamples[WhichSample]+"_"+UBCodeVersion+".root"));

				vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

					TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
					CurrentPlotsReco.push_back(histReco);
			
				}

				PlotsReco.push_back(CurrentPlotsReco);		

			}

			// ----------------------------------------------------------------------------------------------------------------------------

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				// ----------------------------------------------------------------------------------------------------------------

				int NBins = PlotsReco[0][WhichPlot]->GetXaxis()->GetNbins();
				
				if (WhichEventWeightLabel == 0) { 
					
					TotalSystXsecs[WhichPlot].resize(NBins); 
					
					for (int WhichBin = 0; WhichBin < NBins; WhichBin++) {
					
						TotalSystXsecs[WhichPlot][WhichBin] = 0.;
					
					}
					
				}				
				
				std::vector< std::vector<double> > ArrayBinXsecs;
				ArrayBinXsecs.resize(NBins);

				// ----------------------------------------------------------------------------------------------------------------------
		
				TCanvas* PlotCanvas = new TCanvas(EventWeightLabels[WhichEventWeightLabel]+"_"+PlotNames[WhichPlot]+Runs[WhichRun],\
							EventWeightLabels[WhichEventWeightLabel]+"_"+PlotNames[WhichPlot]+Runs[WhichRun],205,34,1024,768);
				PlotCanvas->cd();				

				TPad *midPad = new TPad("midPad", "", 0.005, 0., 0.995, 0.995);
				midPad->SetTopMargin(0.16);
				midPad->SetBottomMargin(0.15);
				midPad->SetLeftMargin(0.2);
				midPad->Draw();

				TLegend* leg = new TLegend(0.18,0.85,0.6,0.98);
				leg->SetBorderSize(0);
				leg->SetTextSize(TextSize);
				leg->SetTextFont(FontStyle);
				leg->SetNColumns(1);

				// ------------------------------------------------------------------------------------------------------------------

				double max = -99.; 
				double min = 1E3; 			

				// Drawing data plots using the efficiencies from nominal & variation samples

				for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

					MakeMyPlotPretty(PlotsReco[WhichSample][WhichPlot]);
					PlotsReco[WhichSample][WhichPlot]->SetLineColor(Colors[WhichSample]);
					PlotsReco[WhichSample][WhichPlot]->SetMarkerStyle(Markers[WhichSample]);
					PlotsReco[WhichSample][WhichPlot]->SetMarkerColor(Colors[WhichSample]);
					PlotsReco[WhichSample][WhichPlot]->SetMarkerSize(2.);

					PlotsReco[WhichSample][WhichPlot]->GetXaxis()->SetTitleSize(TextSize);
					PlotsReco[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(TextSize);
					PlotsReco[WhichSample][WhichPlot]->GetXaxis()->SetTitleOffset(1.);

					PlotsReco[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(1.27);
					PlotsReco[WhichSample][WhichPlot]->GetYaxis()->SetTitle(PlotXAxis[WhichPlot]);
					PlotsReco[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
					PlotsReco[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(TextSize);
					PlotsReco[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(TextSize);
					PlotsReco[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(3);

					double LocalMax = PlotsReco[WhichSample][WhichPlot]->GetMaximum();
					double LocalMin = PlotsReco[WhichSample][WhichPlot]->GetMinimum();				
					max = TMath::Max(LocalMax,max);
					min = TMath::Min(LocalMin,min);				
					PlotsReco[0][WhichPlot]->GetYaxis()->SetRangeUser(0.7*min,1.2*max);

					midPad->cd();
					PlotsReco[WhichSample][WhichPlot]->Draw("hist p0 same");
					
					if (WhichSample == 0) { leg->AddEntry(PlotsReco[WhichSample][WhichPlot],"Nominal","p"); }
					if (WhichSample == 1) 
						{ leg->AddEntry(PlotsReco[WhichSample][WhichPlot],EventWeightLabels[WhichEventWeightLabel],"p"); }

					// 0th element is the reference plot, should not be included in the systematics
					
					if (WhichSample != 0) { 
					
						for (int WhichBin = 0; WhichBin < NBins; WhichBin++){

							ArrayBinXsecs[WhichBin].push_back(PlotsReco[WhichSample][WhichPlot]->GetBinContent(WhichBin+1));

						}
						
					}

				} // End of the loop over the variation samples 

				PlotsReco[0][WhichPlot]->Draw("hist p0 same");
				leg->Draw();	

				TLatex latex;
				latex.SetTextFont(FontStyle);
				latex.SetTextSize(TextSize);

				double tor860_wcut = -99.;

				if (Runs[WhichRun] == "Run1") { tor860_wcut = tor860_wcut_Run1; }
				if (Runs[WhichRun] == "Run2") { tor860_wcut = tor860_wcut_Run2; }
				if (Runs[WhichRun] == "Run3") { tor860_wcut = tor860_wcut_Run3; }
				if (Runs[WhichRun] == "Run4") { tor860_wcut = tor860_wcut_Run4; }
				if (Runs[WhichRun] == "Run5") { tor860_wcut = tor860_wcut_Run5; }

				TString Label = Runs[WhichRun] + " " +ToString(tor860_wcut)+" POT";
				latex.DrawLatexNDC(0.45,0.78, Label);				

				// -------------------------------------------------------------------------------------------------

				// Saving the canvas where the CV & SystVar predictions have been overlaid

				PlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/BeamOn9/Flux_Systematics_"+PlotNames[WhichPlot]+"_"\
						   +Runs[WhichRun]+"_"+EventWeightLabels[WhichEventWeightLabel]+"_"+UBCodeVersion+".pdf");

				delete PlotCanvas;

				// -----------------------------------------------------------------------------------------------------------
				
				// Now plot it using the mean and sigma on a bin by bin basis
				
				TCanvas* MeanStdPlotCanvas = new TCanvas("MeanStd_"+EventWeightLabels[WhichEventWeightLabel]+"_"+\
									   PlotNames[WhichPlot]+Runs[WhichRun],\
									   "MeanStd_"+EventWeightLabels[WhichEventWeightLabel]+"_"+\
									   PlotNames[WhichPlot]+Runs[WhichRun],205,34,1024,768);
				MeanStdPlotCanvas->cd();
				
				TPad *MeanStdmidPad = new TPad("MeanStdmidPad", "", 0.005, 0., 0.995, 0.995);
				MeanStdmidPad->SetTopMargin(0.16);
				MeanStdmidPad->SetBottomMargin(0.15);
				MeanStdmidPad->SetLeftMargin(0.2);
				MeanStdmidPad->Draw();

				TLegend* MeanStdleg = new TLegend(0.16,0.85,0.6,0.98);
				MeanStdleg->SetBorderSize(0);
				MeanStdleg->SetTextSize(TextSize);
				MeanStdleg->SetTextFont(FontStyle);
				MeanStdleg->SetNColumns(1);					
				
				MeanStdmidPad->cd();				
				PlotsReco[0][WhichPlot]->Draw("hist p0 same");
				
				MeanStdleg->AddEntry(PlotsReco[0][WhichPlot],"Nominal","p");
				
				TH1D* MeanStdClone = (TH1D*)(PlotsReco[0][WhichPlot]->Clone());
				
				for (int WhichBin = 0; WhichBin < NBins; WhichBin++) {

					std::vector<double> ArrayBinXsecsBin = ArrayBinXsecs[WhichBin];

					double mean = computeMean(ArrayBinXsecsBin); 
					double std = computeStd(mean,ArrayBinXsecsBin); 					

					MeanStdClone->SetBinContent(WhichBin+1,mean);
					MeanStdClone->SetBinError(WhichBin+1,std);					

				}	

				MeanStdClone->SetLineColor(kGreen+2);
				MeanStdClone->SetMarkerColor(kGreen+2);
				MeanStdClone->SetMarkerStyle(22);				
				MeanStdleg->AddEntry(MeanStdClone,EventWeightLabels[WhichEventWeightLabel],"p");						
//				MeanStdClone->Draw("hist p0 same");
				MeanStdClone->Draw("ex0 same");
				
				MeanStdleg->Draw("same");
				
				MeanStdPlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/BeamOn9/MeanSt_Flux_Systematics_"+PlotNames[WhichPlot]+"_"\
						   +Runs[WhichRun]+"_"+EventWeightLabels[WhichEventWeightLabel]+"_"+UBCodeVersion+".pdf");
						   
				delete MeanStdPlotCanvas;		   
						   
				// -----------------------------------------------------------------------------------------------------------		

				// Store the extracted systematic uncertainty for a given plot & for a given label
				// by using the expression
				// SystUnc = sqrt( (mean - nominal)^2 + std^2 )
				
				std::vector<double> PlotAndLabelSystUnc;
				PlotAndLabelSystUnc.resize(NBins);
				
				for (int WhichBin = 0; WhichBin < NBins; WhichBin++) {
				
					PlotAndLabelSystUnc[WhichBin] = TMath::Sqrt( TMath::Power(PlotsReco[0][WhichPlot]->GetBinContent(WhichBin+1)\
										      - MeanStdClone->GetBinContent(WhichBin+1),2.)\
										      + TMath::Power(MeanStdClone->GetBinError(WhichBin+1),2.) );
										      
					TotalSystXsecs[WhichPlot][WhichBin] = TMath::Sqrt( 
									       TMath::Power(TotalSystXsecs[WhichPlot][WhichBin],2.) 
									       + TMath::Power( PlotAndLabelSystUnc[WhichBin],2.) );


				}

			} // End of the loop over the plots
		
		} // End of the loop over the EventWeight Label
		
		// ----------------------------------------------------------------------------------------------------------------------------
			
		// Loop over the plots to store the relevant uncertainties in the file
		
		SystFile->cd();

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {
		
			TH1D* SystPlot = (TH1D*)PlotsReco[0][WhichPlot]->Clone();
			int NBins = PlotsReco[0][WhichPlot]->GetXaxis()->GetNbins();	

			for (int WhichBin = 1; WhichBin <= NBins; WhichBin++){

				SystPlot->SetBinContent(WhichBin,TotalSystXsecs[WhichPlot][WhichBin]);
				SystPlot->SetBinError(WhichBin,0);			

			}
			
			SystPlot->Write(PlotNames[WhichPlot]);							
		
		}
		
		SystFile->Close();		
		
		cout << endl << "Systematics file " << SystFileName << " has been created" << endl << endl;

	} // End of the loop over the runs	

} // End of the program 
