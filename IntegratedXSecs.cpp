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

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "/home/afroditi/Dropbox/PhD/Secondary_Code/mySimFunctions.cpp"
#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

// --------------------------------------------------------------------------------------------------------------------------------------------

void Multiply(TH1D* True, TH2D* SmearMatrix) {

	int XBins = SmearMatrix->GetXaxis()->GetNbins();
	int YBins = SmearMatrix->GetYaxis()->GetNbins();

	if (XBins != YBins) { std::cout << "Not symmetric matrix" << std::endl; }

	for (int WhichXBin = 0; WhichXBin < XBins; WhichXBin++) {

		double Entry = 0.;

		double TrueInBin = True->GetBinContent(WhichXBin + 1);

		for (int WhichYBin = 0; WhichYBin < YBins; WhichYBin++) {

			double MigrationInBin = SmearMatrix->GetBinContent(WhichXBin + 1,WhichYBin + 1);

			Entry +=  MigrationInBin * TrueInBin;
	
		}

		// Bin entry in reco space

		True->SetBinContent(WhichXBin+1,Entry);

	}

	return;

}

// --------------------------------------------------------------------------------------------------------------------------------------------

void IntegratedXSecs() {

	int DecimalAccuracy = 2;

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;

	TString PathToFiles = "myXSec/";

	// ---------------------------------------------------------------------------------------------------------------------------

//	PlotNames.push_back("DeltaPTPlot"); 
//	PlotNames.push_back("DeltaAlphaTPlot"); 
//	PlotNames.push_back("DeltaPhiTPlot");
//	PlotNames.push_back("MuonMomentumPlot"); 
	PlotNames.push_back("MuonCosThetaPlot"); 
//	PlotNames.push_back("MuonPhiPlot");
//	PlotNames.push_back("ProtonMomentumPlot"); 
//	PlotNames.push_back("ProtonCosThetaPlot");
//	PlotNames.push_back("ProtonPhiPlot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");
//	Runs.push_back("Run2");
	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();

		vector<TString> NameOfSamples; NameOfSamples.clear();
	
		// CV

		NameOfSamples.push_back("Overlay9");

		NameOfSamples.push_back("Genie_v3_0_6_Out_Of_The_Box");
		NameOfSamples.push_back("Genie_v3_0_6_uB_Tune_1");
		NameOfSamples.push_back("SuSav2");
		NameOfSamples.push_back("NuWro");
		NameOfSamples.push_back("GiBUU");
		NameOfSamples.push_back("GENIEv2");
		NameOfSamples.push_back("NEUT");
		NameOfSamples.push_back("GENIEv3_0_4");

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();

		// -------------------------------------------------------------------------------------------------------------------

		// Open the files and grab the relevant plots

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();
			vector<TH1D*> CurrentPlotsCC1pReco; CurrentPlotsCC1pReco.clear();
			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();

			// CV With Statistical Uncertainties

			if (NameOfSamples[WhichSample] == "Overlay9") { // CV with statistical uncertainties only for now

				FileSample.push_back(TFile::Open(PathToFiles+UBCodeVersion+"/WienerSVD_ExtractedXSec_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root")); 

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
					CurrentPlotsReco.push_back(histReco);

//					TH1D* histCC1pReco = (TH1D*)(FileSample[WhichSample]->Get("CC1pReco"+PlotNames[WhichPlot]));
//					CurrentPlotsCC1pReco.push_back(histCC1pReco);

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);
		
				}

			} 

			else {

				if (
					NameOfSamples[WhichSample] == "Genie_v3_0_6_Out_Of_The_Box" || 
					NameOfSamples[WhichSample] == "Genie_v3_0_6_uB_Tune_1" || 
					NameOfSamples[WhichSample] == "SuSav2" ||
					NameOfSamples[WhichSample] == "GENIEv2" ||
					NameOfSamples[WhichSample] == "GENIEv3_0_4"
				) {
					FileSample.push_back(TFile::Open("../myGenieAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root")); 
				}

				if (NameOfSamples[WhichSample] == "NuWro") 
					{ FileSample.push_back(TFile::Open("../myNuWroAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root")); }

				if (NameOfSamples[WhichSample] == "GiBUU") 
					{ FileSample.push_back(TFile::Open("../myGiBUUAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root")); }

				if (NameOfSamples[WhichSample] == "NEUT") 
					{ FileSample.push_back(TFile::Open("../myNEUTAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root")); }

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++){

					TH1D* histReco = nullptr;
					CurrentPlotsReco.push_back(histReco);

					TH1D* histCC1pReco = nullptr;
					CurrentPlotsCC1pReco.push_back(histCC1pReco);

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);
		
				}

			}

			PlotsReco.push_back(CurrentPlotsReco);		
			PlotsCC1pReco.push_back(CurrentPlotsCC1pReco);
			PlotsTrue.push_back(CurrentPlotsTrue);		

		}

		// -----------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

			// -----------------------------------------------------------------------------------------------------------------

			TH2D* Ac = (TH2D*)FileSample[0]->Get("Ac"+PlotNames[WhichPlot]);

			// -----------------------------------------------------------------------------------------------------------------

			cout << PlotNames[WhichPlot] << endl << endl;

			// -----------------------------------------------------------------------------------------------------------------

			// BeamOn Statistical Uncertainty

			cout << "    Data "+ Runs[WhichRun] +" = " << IntegratedXSec(PlotsReco[0][WhichPlot]) << " \\pm " << IntegratedXSecError(PlotsReco[0][WhichPlot]) << endl;

			// -----------------------------------------------------------------------------------------------------------------

//			// Overlay

//			cout << "    Overlay = " << IntegratedXSec(PlotsCC1pReco[0][WhichPlot]) << " \\pm " << IntegratedXSecError(PlotsCC1pReco[0][WhichPlot]) << endl;

			// -----------------------------------------------------------------------------------------------------------------

			// GENIE Overlay

			cout << "    MC "+ Runs[WhichRun] +" = " << IntegratedXSec(PlotsTrue[0][WhichPlot]) << endl;

			// -----------------------------------------------------------------------------------------------------------------

			// GENIE v3.0.6 Out Of The Box

			Multiply(PlotsTrue[1][WhichPlot],Ac);
			cout << "    GENIE v3.0.6 = " << IntegratedXSec(PlotsTrue[1][WhichPlot]) << endl;

			// -----------------------------------------------------------------------------------------------------------------

//			// GENIE v3.0.6 MicroBooNE Tune v1

//			Multiply(PlotsTrue[2][WhichPlot],Ac);
//			cout << "    GENIE v3.0.6 MicroBooNE Tune v1 = " << IntegratedXSec(PlotsTrue[2][WhichPlot]) << endl;

			// -----------------------------------------------------------------------------------------------------------------

			// GENIE SuSav2

			Multiply(PlotsTrue[3][WhichPlot],Ac);
			cout << "    GENIE SuSav2 = " << IntegratedXSec(PlotsTrue[3][WhichPlot]) << endl;

			// -----------------------------------------------------------------------------------------------------------------

			// GENIE v3.0.4

			Multiply(PlotsTrue[8][WhichPlot],Ac);
			cout << "    GENIE v3.0.4 = " << IntegratedXSec(PlotsTrue[8][WhichPlot]) << endl;

			// -----------------------------------------------------------------------------------------------------------------

			// GENIE v2

			Multiply(PlotsTrue[6][WhichPlot],Ac);
			cout << "    GENIE v2.12.10 = " << IntegratedXSec(PlotsTrue[6][WhichPlot]) << endl;

			// -----------------------------------------------------------------------------------------------------------------

			// NuWro

			Multiply(PlotsTrue[4][WhichPlot],Ac);
			cout << "    NuWro = " << IntegratedXSec(PlotsTrue[4][WhichPlot]) << endl;

			// -----------------------------------------------------------------------------------------------------------------

			// GiBUU

			Multiply(PlotsTrue[5][WhichPlot],Ac);
			cout << "    GiBUU = " << IntegratedXSec(PlotsTrue[5][WhichPlot]) << endl;

			// -----------------------------------------------------------------------------------------------------------------

			// NEUT

			Multiply(PlotsTrue[7][WhichPlot],Ac);
			cout << "    NEUT = " << IntegratedXSec(PlotsTrue[7][WhichPlot]) << endl;

			// -----------------------------------------------------------------------------------------------------------------

			cout << endl;
			cout << "----------------------------------------------------------------------------------------------------------------------------" << endl;
			cout << endl << endl ;

		} // End of the loop over the plots

	} // End of the loop over the runs	

} // End of the program 
