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
#include <TVectorD.h>
#include <TMatrixD.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <string>

#include "/home/afroditi/Dropbox/PhD/Secondary_Code/mySimFunctions.cpp"
#include "../myClasses/Constants.h"
#include "../myClasses/Util.h"

using namespace std;
using namespace Constants;

#include "../myClasses/Util.h"

// -----------------------------------------------------------------------------------------------

void WienerSVD_Chi2Covariance(TString Var) {

	int DecimalAccuracy = 2;

	TH1D::SetDefaultSumw2();

	TString PathToFiles = "myXSec/";

	std::cout.precision(3);	

	// ---------------------------------------------------------------------------------------------------------------------------

	vector<TString> PlotNames;

	if (Var == "STV") {

	PlotNames.push_back("DeltaPTPlot"); 
	PlotNames.push_back("DeltaAlphaTPlot"); 
	PlotNames.push_back("DeltaPhiTPlot");

//	PlotNames.push_back("MuonMomentumPlot"); 
//	PlotNames.push_back("MuonCosThetaPlot"); 
//	PlotNames.push_back("ProtonMomentumPlot"); 
//	PlotNames.push_back("ProtonCosThetaPlot");

	} else if (Var == "Long") {

//	PlotNames.push_back("DeltaPLPlot"); 
	PlotNames.push_back("DeltaPnPlot"); 	
	PlotNames.push_back("DeltaPtxPlot"); 
	PlotNames.push_back("DeltaPtyPlot"); 
	PlotNames.push_back("kMissPlot");	
//	PlotNames.push_back("APlot"); 

//	PlotNames.push_back("PMissPlot"); 
//	PlotNames.push_back("PMissMinusPlot");

	} else if (Var == "Kine") {

	PlotNames.push_back("MuonMomentumPlot"); 
	//PlotNames.push_back("MuonCosThetaPlot"); 
//	PlotNames.push_back("MuonPhiPlot");
	//PlotNames.push_back("ProtonMomentumPlot"); 
	//PlotNames.push_back("ProtonCosThetaPlot");
//	PlotNames.push_back("ProtonPhiPlot");

	} else {  

		cout << "What the heck do you want me to compute ?" << endl;
		return;

	}

	const int N1DPlots = PlotNames.size();
	//cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	//Runs.push_back("Run1");
//	Runs.push_back("Run2");
//	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	//cout << "Number of Runs = " << NRuns << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> NameOfSamples; NameOfSamples.clear();
	vector<TString> SampleLabel; SampleLabel.clear();
	
	// CV

//	NameOfSamples.push_back("Overlay9");

/*	NameOfSamples.push_back("GENIEv3_0_4"); SampleLabel.push_back("v3.0.4");
	NameOfSamples.push_back("Genie_v3_0_6_Out_Of_The_Box"); SampleLabel.push_back("v3.0.6");
//	NameOfSamples.push_back("Genie_v3_0_6_uB_Tune_1");
	NameOfSamples.push_back("GENIEv2"); SampleLabel.push_back("v2.12.10");*/
	NameOfSamples.push_back("SuSav2"); SampleLabel.push_back("SuSav2");
/*	NameOfSamples.push_back("GiBUU"); SampleLabel.push_back("GiBUU");
	NameOfSamples.push_back("NEUT"); SampleLabel.push_back("NEUT");
	NameOfSamples.push_back("NuWro"); SampleLabel.push_back("NuWro");
*/
	const int NSamples = NameOfSamples.size();
	
	vector<TFile*> FileSample; FileSample.resize(NSamples);
	
	for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {	
	
		if (
			NameOfSamples[WhichSample] == "Genie_v3_0_6_Out_Of_The_Box" || 
			NameOfSamples[WhichSample] == "Genie_v3_0_6_uB_Tune_1" || 
			NameOfSamples[WhichSample] == "SuSav2" ||
			NameOfSamples[WhichSample] == "GENIEv2" ||
			NameOfSamples[WhichSample] == "GENIEv3_0_4"
		) {
			FileSample[WhichSample] = TFile::Open("../myGenieAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root"); 
		}

		if (NameOfSamples[WhichSample] == "NuWro") 
			{ FileSample[WhichSample] = TFile::Open("../myNuWroAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root"); }

		if (NameOfSamples[WhichSample] == "GiBUU") 
			{ FileSample[WhichSample] = TFile::Open("../myGiBUUAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root"); }

		if (NameOfSamples[WhichSample] == "NEUT") 
			{ FileSample[WhichSample] = TFile::Open("../myNEUTAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root"); }
	
	}		

	// -----------------------------------------------------------------------------------------------------------------------------------------
	
	vector<TFile*> DataFileSample; DataFileSample.resize(NRuns);		
	
	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		DataFileSample[WhichRun] = TFile::Open(PathToFiles+UBCodeVersion+"/WienerSVD_ExtractedXSec_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root");
			
	}		

	// -----------------------------------------------------------------------------------------------------------------------------------------

	double Chi2[NRuns][N1DPlots][NSamples+1];
	int Ndof[NRuns][N1DPlots][NSamples+1];
	double PVal[NRuns][N1DPlots][NSamples+1];	
	
	// -----------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		cout << endl;	

		// -----------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

			// -----------------------------------------------------------------------------------------------------------------

			TH2D* Ac = (TH2D*)DataFileSample[0]->Get("Ac"+PlotNames[WhichPlot]);

			// -----------------------------------------------------------------------------------------------------------------

			//cout << PlotNames[WhichPlot] << endl << endl;

			// -----------------------------------------------------------------------------------------------------------------
			
//TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+"_"+Runs[WhichRun],PlotNames[WhichPlot]+"_"+Runs[WhichRun],205,34,1024,768);
//PlotCanvas->cd();			

			TH1D* DataPlot = (TH1D*)(DataFileSample[WhichRun]->Get("RecoFullUnc"+PlotNames[WhichPlot]));
			TH1D* MCPlot = (TH1D*)(DataFileSample[WhichRun]->Get("True"+PlotNames[WhichPlot]));
			
//DataPlot->Draw("e1x0 same");
//MCPlot->Draw("hist same");

			TH2D* Covariance = (TH2D*)(DataFileSample[WhichRun]->Get("UnfCov"+PlotNames[WhichPlot]));

			TH2D* CovarianceClone = (TH2D*)(Covariance->Clone()); 

			int n = DataPlot->GetXaxis()->GetNbins();

			for (int ix = 1; ix <= n; ix++) {

				for (int iy = 1; iy <= n; iy++) {

					double WidthX = Covariance->GetXaxis()->GetBinWidth(ix);
					double WidthY = Covariance->GetYaxis()->GetBinWidth(iy);

					double TwoDWidth = WidthX * WidthY;
					double BinContent = Covariance->GetBinContent(ix,iy);

					CovarianceClone->SetBinContent(ix,iy,BinContent/TwoDWidth);

				}					

			}			

			// ----------------------------------------------------------------------------------------------------	

			//cout << "MC" << endl;
			CalcChiSquared(MCPlot,DataPlot,CovarianceClone,Chi2[WhichRun][WhichPlot][0],Ndof[WhichRun][WhichPlot][0],PVal[WhichRun][WhichPlot][0]);

			// ----------------------------------------------------------------------------------------------------	

			vector<TH1D*> Plots; Plots.resize(NSamples);
			vector<TH1D*> Clone; Clone.resize(NSamples);

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				//cout << NameOfSamples[WhichSample] << endl;
			
				Plots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot])); 	

				Clone[WhichSample] = Multiply(Plots[WhichSample],Ac);
				CalcChiSquared(Clone[WhichSample],DataPlot,CovarianceClone,Chi2[WhichRun][WhichPlot][WhichSample+1],Ndof[WhichRun][WhichPlot][WhichSample+1],PVal[WhichRun][WhichPlot][WhichSample+1]);
//				CalcChiSquared(Plots[WhichSample],DataPlot,CovarianceClone,Chi2[WhichRun][WhichPlot][WhichSample+1],Ndof[WhichRun][WhichPlot][WhichSample+1],PVal[WhichRun][WhichPlot][WhichSample+1]);						
				
//Clone[WhichSample]->SetLineColor(WhichSample+2);				
//Clone[WhichSample]->Draw("hist same");
//cout << "alternative chi2 = " << Chi2Func(Clone[WhichSample],MCPlot) << endl;

			}

			// ----------------------------------------------------------------------------------------------------			

		} // End of the loop over the plots

		cout << Runs[WhichRun] << endl << endl;
		cout << "MC uB Tune" ;

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

			cout << " & " << Chi2[WhichRun][WhichPlot][0] << "/" << Ndof[WhichRun][WhichPlot][0] ;
			if (WhichPlot == N1DPlots-1) { cout << "  \\\\"; }

		} // End of the loop over the plots

		cout << endl;

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			cout << SampleLabel[WhichSample] ;

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				cout << " & " << Chi2[WhichRun][WhichPlot][WhichSample+1] << "/" << Ndof[WhichRun][WhichPlot][WhichSample+1] ;
				if (WhichPlot == N1DPlots-1) { cout << "  \\\\"; }

			} // End of the loop over the plots

			cout << endl;

		} // End of the loop over the samples

		cout << endl;

	} // End of the loop over the runs	

} // End of the program 
