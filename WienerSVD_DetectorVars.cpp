#include <TFile.h>
#include <TF1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "ubana/AnalysisCode/Secondary_Code/GlobalSettings.cpp"
#include "ubana/AnalysisCode/Secondary_Code/myFunctions.cpp"

#include "ubana/myClasses/Util.h"

// --------------------------------------------------------------------------------------------------------------------------------------------

void FracUnc(TH1D* CV, TH1D* Vars, int index, TLegend* leg, TString LegLabel) {

	std::vector<int> Colors = {kBlack,kRed+1,kGreen+2,kOrange+1,kBlue-3,kMagenta, kCyan, kGreen, kGray+1};

	TH1D* CVClone = (TH1D*)(CV->Clone());
	CVClone->Add(Vars,-1);
	CVClone->Divide(CV);
	CVClone->Scale(100.); // %

	int NBins = CVClone->GetXaxis()->GetNbins();

	// Loop over the bins and, if there are negative entries, flip them

	for (int i = 1; i < NBins + 1; i++) {

		double BinContent = CVClone->GetBinContent(i);

		if (BinContent < 0) { CVClone->SetBinContent(i,-BinContent); }

	}

	CVClone->GetXaxis()->CenterTitle();
	CVClone->GetXaxis()->SetTitleFont(FontStyle);
	CVClone->GetXaxis()->SetTitleSize(TextSize);
	CVClone->GetXaxis()->SetLabelFont(FontStyle);
	CVClone->GetXaxis()->SetLabelSize(TextSize);
	CVClone->GetXaxis()->SetNdivisions(8);

	CVClone->GetYaxis()->CenterTitle();
	CVClone->GetYaxis()->SetTitleFont(FontStyle);
	CVClone->GetYaxis()->SetTitleSize(TextSize);
	CVClone->GetYaxis()->SetLabelFont(FontStyle);
	CVClone->GetYaxis()->SetLabelSize(TextSize);
	CVClone->GetYaxis()->SetNdivisions(6);
	CVClone->GetYaxis()->SetTitleOffset(0.95);
	CVClone->GetYaxis()->SetTitle("Uncertainty [%]");

	CVClone->GetYaxis()->SetRangeUser(0.,19);
	TString PlotName = CVClone->GetName();
	if (string(PlotName).find("SingleBin") != std::string::npos) { CVClone->GetYaxis()->SetRangeUser(0.,4.); }

	CVClone->SetMarkerColor(Colors[index]);
	CVClone->SetMarkerStyle(20);
	CVClone->SetMarkerSize(2.);
	CVClone->SetLineColor(Colors[index]);
	CVClone->SetLineWidth(3);

	CVClone->Draw("hist same");

	leg->AddEntry(CVClone,LegLabel,"l");

}

// --------------------------------------------------------------------------------------------------------------------------------------------

TH1D* Multiply(TH1D* True, TH2D* SmearMatrix) {

	TH1D* TrueClone = (TH1D*)(True->Clone());

	int XBins = SmearMatrix->GetXaxis()->GetNbins();
	int YBins = SmearMatrix->GetYaxis()->GetNbins();

	if (XBins != YBins) { std::cout << "Not symmetric matrix" << std::endl; }

	TVectorD signal(XBins);
	TMatrixD response(XBins,XBins);

	H2V(True, signal);
	H2M(SmearMatrix, response, kFALSE); // X axis: Reco, Y axis: True

	TVectorD RecoSpace = response * signal;
	V2H(RecoSpace, TrueClone);	

	return TrueClone;

}

// --------------------------------------------------------------------------------------------------------------------------------------------

// TString Syst = "Stat" "POT" "NTarget" "LY" "TPC" "SCERecomb2" "XSec" "G4" "Flux" "Dirt" "MC_Stat"

void WienerSVD_DetectorVars() {

	// -------------------------------------------------------------------------------------

	GlobalSettings();
	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	TGaxis::SetMaxDigits(3);

	TString ExactFileLocation = PathToFiles+CutExtension;

	// -------------------------------------------------------------------------------------

//	vector<TString> PlotNames;
//	PlotNames.push_back("DeltaPTPlot"); 
//	PlotNames.push_back("DeltaAlphaTPlot"); 
//	PlotNames.push_back("DeltaPhiTPlot"); 
//	PlotNames.push_back("MuonMomentumPlot"); 
//	PlotNames.push_back("MuonPhiPlot"); 
//	PlotNames.push_back("MuonCosThetaPlot");
//	PlotNames.push_back("MuonCosThetaSingleBinPlot");
//	PlotNames.push_back("ProtonMomentumPlot"); 
//	PlotNames.push_back("ProtonPhiPlot"); 
//	PlotNames.push_back("ProtonCosThetaPlot");

//	PlotNames.push_back("CCQEMuonMomentumPlot"); 
//	PlotNames.push_back("CCQEMuonCosThetaPlot"); 
//	PlotNames.push_back("CCQEProtonMomentumPlot"); 
//	PlotNames.push_back("CCQEProtonCosThetaPlot");

	const int NPlots = PlotNames.size();

	// -------------------------------------------------------------------------------------------------------------------------------------

	vector<TString> Runs; Runs.clear();

	Runs.push_back("Run3");

	const int NRuns = Runs.size();

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	// Base Plots

	vector <TH1D*> TrueCC1pPlots; TrueCC1pPlots.resize(NPlots); 
	vector <TH2D*> CC1pResponseMatrix; CC1pResponseMatrix.resize(NPlots); 
	vector <TH1D*> ClosureTestCC1pPlots; ClosureTestCC1pPlots.resize(NPlots); 
	vector <TH1D*> CC1pPlots; CC1pPlots.resize(NPlots); 
	vector <TH1D*> NonCC1pPlots; NonCC1pPlots.resize(NPlots);
	vector<vector <TH2D*> > FracCovariances; FracCovariances.resize(NRuns,vector<TH2D*>(NPlots));

	// -------------------------------------------------------------------------------------------------------------------------------------

	// Alternative MC Models

	std::vector<TString> CV;
	std::vector<TString> Vars;
	std::vector<int> Colors;

	Vars.push_back("LYDown"); Colors.push_back(kRed+1); CV.push_back("CV");
	Vars.push_back("LYRayleigh"); Colors.push_back(kGreen+2); CV.push_back("CV");
	Vars.push_back("LYAttenuation"); Colors.push_back(kOrange+1); CV.push_back("CV");

	Vars.push_back("X"); Colors.push_back(kRed+1); CV.push_back("CV");
	Vars.push_back("YZ"); Colors.push_back(kGreen+2); CV.push_back("CV");
	Vars.push_back("ThetaXZ"); Colors.push_back(kOrange+1); CV.push_back("CV");
	Vars.push_back("ThetaYZ"); Colors.push_back(kBlue-3); CV.push_back("CV");

	Vars.push_back("SCE"); Colors.push_back(kBlue); CV.push_back("CVextra");
	Vars.push_back("Recombination2"); Colors.push_back(kMagenta); CV.push_back("CVextra");

	int NVars = Vars.size();

	// -------------------------------------------------------------------------------------------------------------------------------------

	// Base Samples

	TFile* CVTrueMCFile[NRuns][NVars];
	TFile* CVMCFile[NRuns][NVars];
	TFile* CVResponseMatrixFile[NRuns][NVars];

	// Alternative Samples

	TFile* VarsMCFile[NRuns][NVars];
	TFile* VarsResponseMatrixFile[NRuns][NVars];

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	// Base Plots

	TH2D* CVCC1pResponseMatrix[NRuns][NPlots][NVars];
	TH1D* CVForwardFoldedCC1pPlots[NRuns][NPlots][NVars];
//	TH1D* CVCC1pPlots[NRuns][NPlots][NVars];
	TH1D* CVNonCC1pPlots[NRuns][NPlots][NVars];
	TH1D* CVTrueCC1pPlots[NRuns][NPlots][NVars]; 

	// Alternative Plots

	TH2D* VarsCC1pResponseMatrix[NRuns][NPlots][NVars];
	TH1D* VarsForwardFoldedCC1pPlots[NRuns][NPlots][NVars];
//	TH1D* VarsCC1pPlots[NRuns][NPlots][NVars];
	TH1D* VarsNonCC1pPlots[NRuns][NPlots][NVars];

	// ---------------------------------------------------------------------------------------------------------------------------------------------

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Loop over the samples and open them

		for (int WhichVar = 0; WhichVar < NVars; WhichVar++) {

			// Corresponding CV

			TString CVTStringBaseMC = ExactFileLocation+"/STVStudies_Overlay9_"+Runs[WhichRun]+"_"+CV[WhichVar]+CutExtension+".root"; 
			TString CVTrueTStringBaseMC = PathToFiles+"/TruthSTVAnalysis_Overlay9_"+Runs[WhichRun]+"_"+CV[WhichVar]+"_"+UBCodeVersion+".root"; 
			TString CVResponseFileName = MigrationMatrixPath+"FileResponseMatrices_Overlay9_"+Runs[WhichRun]+"_"+CV[WhichVar]+"_"+UBCodeVersion+".root";

			CVResponseMatrixFile[WhichRun][WhichVar] = new TFile(CVResponseFileName,"readonly");
			CVTrueMCFile[WhichRun][WhichVar] = TFile::Open(CVTrueTStringBaseMC,"readonly");
			CVMCFile[WhichRun][WhichVar] = TFile::Open(CVTStringBaseMC,"readonly");

			// Detector Variations

			TString VarsTStringBaseMC = ExactFileLocation+"/STVStudies_Overlay9_"+Runs[WhichRun]+"_"+Vars[WhichVar]+CutExtension+".root"; 
			TString VarsResponseFileName = MigrationMatrixPath+"FileResponseMatrices_Overlay9_"+Runs[WhichRun]+"_"+Vars[WhichVar]+"_"+UBCodeVersion+".root";

			VarsResponseMatrixFile[WhichRun][WhichVar] = new TFile(VarsResponseFileName,"readonly");
			VarsMCFile[WhichRun][WhichVar] = TFile::Open(VarsTStringBaseMC,"readonly");

		}

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------

		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot++) {

			// --------------------------------------------------------------------------------------------------------------------------------------------------------------

			// Grab the plots from each CV / variation

			for (int WhichVar = 0; WhichVar < NVars; WhichVar++) {

				// Corresponding CV

				CVCC1pResponseMatrix[WhichRun][WhichPlot][WhichVar] = (TH2D*)(CVResponseMatrixFile[WhichRun][WhichVar]->Get("POTScaledCC1pReco"+PlotNames[WhichPlot]+"2D"));
				CVTrueCC1pPlots[WhichRun][WhichPlot][WhichVar] = (TH1D*)(CVTrueMCFile[WhichRun][WhichVar]->Get("True"+PlotNames[WhichPlot]));
//				CVCC1pPlots[WhichRun][WhichPlot][WhichVar] = (TH1D*)(CVMCFile[WhichRun][WhichVar]->Get("CC1pReco"+PlotNames[WhichPlot]));
				CVNonCC1pPlots[WhichRun][WhichPlot][WhichVar] = (TH1D*)(CVMCFile[WhichRun][WhichVar]->Get("NonCC1pReco"+PlotNames[WhichPlot]));
	
				CVForwardFoldedCC1pPlots[WhichRun][WhichPlot][WhichVar] = Multiply(CVTrueCC1pPlots[WhichRun][WhichPlot][WhichVar],CVCC1pResponseMatrix[WhichRun][WhichPlot][WhichVar]);
				CVForwardFoldedCC1pPlots[WhichRun][WhichPlot][WhichVar]->Add(CVNonCC1pPlots[WhichRun][WhichPlot][WhichVar]);

				// Detector Variations	

				VarsCC1pResponseMatrix[WhichRun][WhichPlot][WhichVar] = (TH2D*)(VarsResponseMatrixFile[WhichRun][WhichVar]->Get("POTScaledCC1pReco"+PlotNames[WhichPlot]+"2D"));
//				VarsCC1pPlots[WhichRun][WhichPlot][WhichVar] = (TH1D*)(VarsMCFile[WhichRun][WhichVar]->Get("CC1pReco"+PlotNames[WhichPlot]));
				VarsNonCC1pPlots[WhichRun][WhichPlot][WhichVar] = (TH1D*)(VarsMCFile[WhichRun][WhichVar]->Get("NonCC1pReco"+PlotNames[WhichPlot]));
	
				VarsForwardFoldedCC1pPlots[WhichRun][WhichPlot][WhichVar] = Multiply(CVTrueCC1pPlots[WhichRun][WhichPlot][WhichVar],VarsCC1pResponseMatrix[WhichRun][WhichPlot][WhichVar]);
				VarsForwardFoldedCC1pPlots[WhichRun][WhichPlot][WhichVar]->Add(VarsNonCC1pPlots[WhichRun][WhichPlot][WhichVar]);			

			}

			// --------------------------------------------------------------------------------------------------------------------------------------------------------------

			// Create canvas & legend for a given plot

			TString CanvasName = Runs[WhichRun]+"_"+PlotNames[WhichPlot];
			TCanvas* Canvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
			Canvas->SetBottomMargin(0.16);
			Canvas->SetTopMargin(0.12);
			Canvas->SetLeftMargin(0.15);

			TLegend* leg = new TLegend(0.15,0.89,0.9,0.99);
			leg->SetBorderSize(0);
			leg->SetTextFont(FontStyle);
			leg->SetTextSize(TextSize-0.03);
			leg->SetNColumns(3);

			// --------------------------------------------------------------------------------------------------------------------------------------------------------------

			// Plot in the canvas the fractional contribution for each variation

			for (int WhichVar = 0; WhichVar < NVars; WhichVar++) {

				FracUnc(CVForwardFoldedCC1pPlots[WhichRun][WhichPlot][WhichVar],VarsForwardFoldedCC1pPlots[WhichRun][WhichPlot][WhichVar], WhichVar,leg,Vars[WhichVar]);

			}

			leg->Draw();

			Canvas->SaveAs(PlotPath+"Overlay9/WienerSVD_DetVars_FracUnc_"+PlotNames[WhichPlot]+"Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
			delete Canvas;

			// --------------------------------------------------------------------------------------------------------------------------------------------------------------

		} // End of the loop over the plots
		
		// ---------------------------------------------------------------------------------------	

		// Loop over the samples and open them

		for (int WhichVar = 0; WhichVar < NVars; WhichVar++) {

			CVResponseMatrixFile[WhichRun][WhichVar]->Close();
			CVTrueMCFile[WhichRun][WhichVar]->Close();
			CVMCFile[WhichRun][WhichVar]->Close();

			VarsResponseMatrixFile[WhichRun][WhichVar]->Close();
			VarsMCFile[WhichRun][WhichVar]->Close();

		}

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------

		cout << endl << "Detector variation contributions assessed!" << endl << endl;

		// --------------------------------------------------------------------------------------------------------------------------------------------------------------

	} // End of the loop over the runs	

} // End of the program
