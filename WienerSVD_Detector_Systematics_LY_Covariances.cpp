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

#include "../Secondary_Code/CenterAxisTitle.cpp"
#include "../Secondary_Code/SetOffsetAndSize.cpp"
#include "../Secondary_Code/MakeMyPlotPretty.cpp"
#include "../Secondary_Code/myFunctions.cpp"

#include "../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

void WienerSVD_Detector_Systematics_LY_Covariances(TString NomOverlaySample = "Overlay9",TString BeamOnSample = "BeamOn9",TString BeamOffSample = "ExtBNB9",TString DirtSample = "OverlayDirt9") {

	TH2D::SetDefaultSumw2();
	TH1D::SetDefaultSumw2();
	
//	TGaxis::SetMaxDigits(3);
//	TGaxis::SetExponentOffset(-0.1, 1., "y");	
	
	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();
	
	// -------------------------------------------------------------------------------------

	TString CutExtension = "_NoCuts";

	vector<TString> VectorCuts; VectorCuts.clear();
	// v52
	VectorCuts.push_back("");
	VectorCuts.push_back("_PID");
	VectorCuts.push_back("_NuScore");

	int NCuts = (int)(VectorCuts.size());	

	for (int i = 0; i < NCuts; i++) {

		CutExtension = CutExtension + VectorCuts[i];

	}
	
	TString ExactFileLocation = PathToFiles+CutExtension;			

	// ---------------------------------------------------------------------------------------------------------------------	

	vector<TString> PlotNames;

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
//	PlotNames.push_back("EQEPlot"); 
//	PlotNames.push_back("Q2Plot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	// ---------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs;
	Runs.push_back("Run1");
//	Runs.push_back("Run2");
	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");				

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;	
	
	// -----------------------------------------------------------------------------------------				

	vector<TString> NameOfSamples;
	NameOfSamples.push_back("CV"); // Reference plot

	// LY Detector Variations

	NameOfSamples.push_back("LYDown");
	NameOfSamples.push_back("LYRayleigh");
	NameOfSamples.push_back("LYAttenuation");

	const int NSamples = NameOfSamples.size();	

	// ----------------------------------------------------------------------------------------------------------------------

	// Covariance matrices for each run / detector variation sample / plot

	std::vector< std::vector< std::vector<TH2D*> > > DetLYCovarianceMatrix; DetLYCovarianceMatrix.clear(); 
	
	// -------------------------------------------------------------------------------------

	// Files for each run / detector variation sample / plot

	vector< vector<TFile*>> MCFileSample;
	MCFileSample.resize(NRuns, vector<TFile*>(NSamples));

	vector< vector<TFile*>> BeamOnFileSample;
	BeamOnFileSample.resize(NRuns, vector<TFile*>(NSamples));

	vector< vector<TFile*>> BeamOffFileSample;
	BeamOffFileSample.resize(NRuns, vector<TFile*>(NSamples));

	vector< vector<TFile*>> DirtFileSample;
	DirtFileSample.resize(NRuns, vector<TFile*>(NSamples));	

	// -----------------------------------------------------------------------------------------------------------------------

	cout << endl;

	// -----------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// ------------------------------------------------------------------------------------------------------------------

		// Plots for detector variation samples / plots

		vector< vector <TH1D*> > CC1pPlots;
		CC1pPlots.resize(NSamples, vector<TH1D*>(N1DPlots));

		vector< vector <TH1D*> > NonCC1pPlots;
		NonCC1pPlots.resize(NSamples, vector<TH1D*>(N1DPlots));

		vector< vector <TH1D*> > BeamOnPlots;
		BeamOnPlots.resize(NSamples, vector<TH1D*>(N1DPlots));
		
		vector< vector <TH1D*> > BkgSubBeamOnPlots;
		BkgSubBeamOnPlots.resize(NSamples, vector<TH1D*>(N1DPlots));		

		vector< vector <TH1D*> > BeamOffPlots;
		BeamOffPlots.resize(NSamples, vector<TH1D*>(N1DPlots));

		vector< vector <TH1D*> > DirtPlots;
		DirtPlots.resize(NSamples, vector<TH1D*>(N1DPlots));

		vector< vector <TH2D*> > Covariances;
		Covariances.resize(NSamples, vector<TH2D*>(N1DPlots));
		
		vector <TH2D*> TotalCovariances;
		TotalCovariances.resize(N1DPlots);				
		
		// ------------------------------------------------------------------------------------------------------------------

		TString FileName = MigrationMatrixPath+"WienerSVD_Detector_Systematics_LY_"+NomOverlaySample+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
		TFile* SystFile = new TFile(FileName,"recreate");

		// -----------------------------------------------------------------------------------------				

		TFile* BeamOnFile = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOnSample+"_"+Runs[WhichRun]+CutExtension+".root");
		TFile* BeamOffFile = TFile::Open(ExactFileLocation+"/STVStudies_"+BeamOffSample+"_"+Runs[WhichRun]+CutExtension+".root");
		TFile* DirtFile = TFile::Open(ExactFileLocation+"/STVStudies_"+DirtSample+"_"+Runs[WhichRun]+CutExtension+".root");

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {
		
			TFile* MCFile = TFile::Open(ExactFileLocation+"/STVStudies_"+NomOverlaySample+"_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+CutExtension+".root");			
		
			BeamOnFileSample[WhichRun][WhichSample] = BeamOnFile;
			BeamOffFileSample[WhichRun][WhichSample] = BeamOffFile;
			DirtFileSample[WhichRun][WhichSample] = DirtFile;
			MCFileSample[WhichRun][WhichSample] = MCFile;			

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot++) {

				BeamOnPlots[WhichSample][WhichPlot] = (TH1D*)(BeamOnFileSample[WhichRun][WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
				BeamOffPlots[WhichSample][WhichPlot] = (TH1D*)(BeamOffFileSample[WhichRun][WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
				CC1pPlots[WhichSample][WhichPlot] = (TH1D*)(MCFileSample[WhichRun][WhichSample]->Get("CC1p"+PlotNames[WhichPlot]));
				NonCC1pPlots[WhichSample][WhichPlot] = (TH1D*)(MCFileSample[WhichRun][WhichSample]->Get("NonCC1p"+PlotNames[WhichPlot]));
				DirtPlots[WhichSample][WhichPlot] = (TH1D*)(DirtFileSample[WhichRun][WhichSample]->Get("Reco"+PlotNames[WhichPlot]));	
				
				if (BeamOnSample == "BeamOn9") {  
				
					BkgSubBeamOnPlots[WhichSample][WhichPlot] = (TH1D*)(BeamOnPlots[WhichSample][WhichPlot]->Clone());
					BkgSubBeamOnPlots[WhichSample][WhichPlot]->Add(BeamOffPlots[WhichSample][WhichPlot]);
					BkgSubBeamOnPlots[WhichSample][WhichPlot]->Add(NonCC1pPlots[WhichSample][WhichPlot]);				
					BkgSubBeamOnPlots[WhichSample][WhichPlot]->Add(DirtPlots[WhichSample][WhichPlot]);				
				
				} else {

					BkgSubBeamOnPlots[WhichSample][WhichPlot] = (TH1D*)(CC1pPlots[WhichSample][WhichPlot]->Clone());			
				
				}	
	

				// -------------------------------------------------------------------------------------------------------
				
				int NBins = BkgSubBeamOnPlots[WhichSample][WhichPlot]->GetXaxis()->GetNbins();
				
				for (int WhichXBin = 1; WhichXBin <= NBins; WhichXBin++) {
				
					double NomDataEntryX = BkgSubBeamOnPlots[0][WhichPlot]->GetBinContent(WhichXBin);
					double NomDataErrorX = BkgSubBeamOnPlots[0][WhichPlot]->GetBinError(WhichXBin);
					
					double VarDataEntryX = BkgSubBeamOnPlots[WhichSample][WhichPlot]->GetBinContent(WhichXBin);
					double VarDataErrorX = BkgSubBeamOnPlots[WhichSample][WhichPlot]->GetBinError(WhichXBin);
				
					for (int WhichYBin = 1; WhichYBin <= NBins; WhichYBin++) {	
					
						double NomDataEntryY = BkgSubBeamOnPlots[0][WhichPlot]->GetBinContent(WhichYBin);
						double NomDataErrorY = BkgSubBeamOnPlots[0][WhichPlot]->GetBinError(WhichYBin);
						
						double VarDataEntryY = BkgSubBeamOnPlots[WhichSample][WhichPlot]->GetBinContent(WhichYBin);
						double VarDataErrorY = BkgSubBeamOnPlots[WhichSample][WhichPlot]->GetBinError(WhichYBin);
				
						// CV

						double CovEntry = TMath::Max((NomDataEntryX - VarDataEntryX) * (NomDataEntryY - VarDataEntryY),1E-8);
						Covariances[WhichSample][WhichPlot]->SetBinContent(WhichXBin+1,WhichYBin+1,CovEntry);
						
						// Error

	//					double CovError = TMath::Sqrt( 
	//					TMath::Power(DataEntryY - POTDataEntryY,2.) * ( TMath::Power(DataErrorX,2.) + TMath::Power(POTDataErrorX,2.) ) +
	//					TMath::Power(DataEntryX - POTDataEntryX,2.) * ( TMath::Power(DataErrorY,2.) + TMath::Power(POTDataErrorY,2.) )
	//								);

						double CovError = 0.000001;

						Covariances[WhichSample][WhichPlot]->SetBinError(WhichXBin+1,WhichYBin+1,CovError);				
				
					} // End of the loop over the Y bins
				
				} // End of the loop over the X bins
			
				// -------------------------------------------------------------------------------------------------------			
				
				SystFile->cd();
				TString XTitle = BkgSubBeamOnPlots[0][WhichPlot]->GetXaxis()->GetTitle();
				Covariances[WhichSample][WhichPlot]->GetXaxis()->SetTitle("i bin "+XTitle);
				Covariances[WhichSample][WhichPlot]->GetYaxis()->SetTitle("i bin "+XTitle);			
				Covariances[WhichSample][WhichPlot]->Write("Covariance_"+PlotNames[WhichPlot]+NameOfSamples[WhichSample]);
				
				if (WhichSample == 0) { TotalCovariances[WhichPlot] = Covariances[WhichSample][WhichPlot]; }
				else { TotalCovariances[WhichPlot]->Add(Covariances[WhichSample][WhichPlot]); }
			
			} // End of the loop over the plots		

		} // End of the loop over the detector variation samples
		
		// ----------------------------------------------------------------------------------------------------		
		
		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot++) {	
		
			SystFile->cd();
			TotalCovariances[WhichPlot]->Write("Covariance_"+PlotNames[WhichPlot]);			
				
			TString TStringPlotCanvas = "Total_"+PlotNames[WhichPlot]+Runs[WhichRun];		
			TCanvas* PlotCanvas = new TCanvas(TStringPlotCanvas,TStringPlotCanvas,205,34,1024,768);
			PlotCanvas->cd();
			TotalCovariances[WhichPlot]->Draw("coltz");		

		} // End of the loop over the plots

		// ----------------------------------------------------------------------------------------------------

		cout << endl << "LY Detector Systematics file " << FileName << " has been created" << endl << endl;

		// ----------------------------------------------------------------------------------------------------

	} // End of the loop over the runs	

} // End of the program 
