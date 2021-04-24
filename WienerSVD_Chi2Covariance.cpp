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
#include <string>

#include "/home/afroditi/Dropbox/PhD/Secondary_Code/mySimFunctions.cpp"
#include "../myClasses/Constants.h"
#include "../myClasses/Util.h"

using namespace std;
using namespace Constants;

// -----------------------------------------------------------------------------------------------

void CalcChiSquared(TH1D* h_model, TH1D* h_data, TH2D* cov, double &chi, int &ndof, double &pval) {

    // Clone them so we can scale them 

    TH1D* h_model_clone = (TH1D*)h_model->Clone();
    TH1D* h_data_clone  = (TH1D*)h_data->Clone();
    TH2D* h_cov_clone   = (TH2D*)cov->Clone();

    // Getting covariance matrix in TMatrix form

    TMatrix cov_m;
    cov_m.Clear();
    cov_m.ResizeTo(h_cov_clone->GetNbinsX(), h_cov_clone->GetNbinsX());

    // loop over rows

    for (int i = 0; i < h_cov_clone->GetNbinsX(); i++) {

        // loop over columns

        for (int j = 0; j < h_cov_clone->GetNbinsY(); j++) {

            cov_m[i][j] = h_cov_clone->GetBinContent(i+1, j+1);
        
		}
    
    }

    // Inverting the covariance matrix

    TMatrix inverse_cov_m = cov_m.Invert();

    // Calculating the chi2 = Summation_ij{ (x_i - mu_j)*E_ij^(-1)*(x_j - mu_j)  }
    // x = data, mu = model, E^(-1) = inverted covariance matrix 

    chi = 0;
    
	for (int i = 0; i < h_cov_clone->GetNbinsX(); i++) {

        for (int j = 0; j < h_cov_clone->GetNbinsY(); j++) {

            chi += (h_data_clone->GetBinContent(i+1) - h_model_clone->GetBinContent(i+1)) * inverse_cov_m[i][j] * (h_data_clone->GetBinContent(j+1) - h_model_clone->GetBinContent(j+1));
 
        }
 
    }

    ndof = h_data_clone->GetNbinsX();
    pval = TMath::Prob(chi, ndof);

    std::cout << "Chi2/dof: " << chi << "/" << h_data_clone->GetNbinsX() << " = " << chi/double(ndof) <<  std::endl;
    std::cout << "p-value: " <<  pval << "\n" <<  std::endl;

    delete h_model_clone;
    delete h_data_clone;
    delete h_cov_clone;

}

// -----------------------------------------------------------------------------------------------

void WienerSVD_Chi2Covariance() {

	int DecimalAccuracy = 2;

	TH1D::SetDefaultSumw2();
	vector<TString> PlotNames;

	TString PathToFiles = "myXSec/";
	TString PathToCovFiles = "myMigrationMatrices/";	

	// ---------------------------------------------------------------------------------------------------------------------------

	PlotNames.push_back("DeltaPTPlot"); 
//	PlotNames.push_back("DeltaAlphaTPlot"); 
//	PlotNames.push_back("DeltaPhiTPlot");
//	PlotNames.push_back("MuonMomentumPlot"); 
//	PlotNames.push_back("MuonCosThetaPlot"); 
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
//	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	// -----------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> NameOfSamples; NameOfSamples.clear();
	
	// CV

//	NameOfSamples.push_back("Overlay9");

//	NameOfSamples.push_back("Genie_v3_0_6_Out_Of_The_Box");
//	NameOfSamples.push_back("Genie_v3_0_6_uB_Tune_1");
	NameOfSamples.push_back("SuSav2");
	// NameOfSamples.push_back("NuWro");
	// NameOfSamples.push_back("GiBUU");
	// NameOfSamples.push_back("GENIEv2");
	// NameOfSamples.push_back("NEUT");
	// NameOfSamples.push_back("GENIEv3_0_4");

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
	vector<TFile*> CovarianceFiles; CovarianceFiles.resize(NRuns);			
	
	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		DataFileSample[WhichRun] = TFile::Open(PathToFiles+UBCodeVersion+"/WienerSVD_ExtractedXSec_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root");
			
		CovarianceFiles[WhichRun] = TFile::Open(PathToCovFiles+UBCodeVersion+"/WienerSVD_Total_CovarianceMatrices_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root");
			
	}		
	
	// -----------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// -----------------------------------------------------------------------------------------------------------------------------

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

			// -----------------------------------------------------------------------------------------------------------------

			cout << PlotNames[WhichPlot] << endl << endl;

			// -----------------------------------------------------------------------------------------------------------------
			
			TH1D* DataPlot = (TH1D*)(DataFileSample[WhichRun]->Get("Reco"+PlotNames[WhichPlot]));
			TH1D* MCPlot = (TH1D*)(DataFileSample[WhichRun]->Get("True"+PlotNames[WhichPlot]));
			
			TH2D* Covariance = (TH2D*)(CovarianceFiles[WhichRun]->Get("TotalCovariance_"+PlotNames[WhichPlot]));

			TH2D* CovarianceClone = (TH2D*)(Covariance->Clone()); 

			int n = DataPlot->GetXaxis()->GetNbins();

			for (int ix = 1; ix <= n; ix++) {

				double WidthX = Covariance->GetXaxis()->GetBinWidth(ix);

				for (int iy = 1; iy <= n; iy++) {

					double WidthY = Covariance->GetYaxis()->GetBinWidth(iy);

					double TwoDWidth = WidthX * WidthY;
					double BinContent = Covariance->GetBinContent(ix,iy);

					CovarianceClone->SetBinContent(ix,iy,BinContent/TwoDWidth);

				}					

			}			

			// ----------------------------------------------------------------------------------------------------
					
			vector<TH1D*> Plots; Plots.resize(NSamples);					
					
			for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

				cout << NameOfSamples[WhichSample] << endl;
			
				Plots[WhichSample] = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot])); 

				double chi = -99.;
				int ndof = -99;
				double pval = -99.;	

				CalcChiSquared(Plots[WhichSample],DataPlot,CovarianceClone,chi,ndof,pval);							
				
			}

			// ----------------------------------------------------------------------------------------------------


			// ----------------------------------------------------------------------------------------------------			

		} // End of the loop over the plots

	} // End of the loop over the runs	

} // End of the program 
