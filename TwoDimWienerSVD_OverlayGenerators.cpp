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
#include <TMatrixD.h>
#include <TVectorD.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include  "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/mySimFunctions.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/MakeMyPlotPretty.cpp"

#include "../myClasses/Tools.h"
#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "../myClasses/Util.h"

//----------------------------------------//

int LocateBinWithValue(TH1D* h, double Value) {

	int NBins = h->GetXaxis()->GetNbins();

	for (int i = 1; i <= NBins; i++) {

		double CurrentEntry = h->GetBinContent(i);
		if (CurrentEntry == Value) { return i; } 

	}

	return -99;

}

//----------------------------------------//

void PrettyPlot(TH1D* h) {

	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetTitleSize(0.06);
	h->GetXaxis()->SetLabelSize(0.06);
	h->GetXaxis()->SetTitleOffset(1.05);
	h->GetXaxis()->SetNdivisions(8);
	h->GetXaxis()->CenterTitle();

	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetNdivisions(8);
	h->GetYaxis()->SetTitleOffset(1.35);
	h->GetYaxis()->SetTitleSize(0.06);
	h->GetYaxis()->SetLabelSize(0.06);
	h->GetYaxis()->CenterTitle();

}

//----------------------------------------//

void TwoDimWienerSVD_OverlayGenerators(bool PlotGENIE = true, bool PlotGen = false, 
								 	   bool PlotGENIEFSITweaks = false, bool PlotGENIEFlagTweaks = false, 
								 	   bool PlotGENIECT = false, bool PlotNuclModels = false, 
								 	   bool PlotNuWro = false, bool PlotNominal = false, 
									   bool GiBUUComp = false, bool All = false
									) {

	//----------------------------------------//

	Tools tools;

	//----------------------------------------//

	int DecimalAccuracy = 2;

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(4);		

	TString PathToFiles = "myXSec/";

	TString Extra = "";
	if (!PlotGENIE && PlotGen) { Extra = "OtherGen"; }
	if (PlotGENIE && PlotGen) { Extra = "All"; }
	if (!PlotGENIE && !PlotGen && PlotGENIEFSITweaks) { Extra = "GENIEFSITweaks"; }
	if (!PlotGENIE && !PlotGen && PlotGENIEFlagTweaks) { Extra = "GENIEFlagTweaks"; }
	if (!PlotGENIE && !PlotGen && PlotGENIECT) { Extra = "GENIEClosureTest"; }
	if (!PlotGENIE && !PlotGen && PlotNuclModels) { Extra = "GENIENuclModels"; }
	if (!PlotGENIE && !PlotGen && PlotNuWro) { Extra = "NuWro"; }
	if (!PlotGENIE && !PlotGen && PlotNominal) { Extra = "Nominal"; }
	if (!PlotGENIE && !PlotGen && GiBUUComp) { Extra = "GiBUUComp"; }		
	if (!PlotGENIE && !PlotGen && All) { Extra = "AllGens"; }			

	//----------------------------------------//

	vector<TString> PlotNames;
	vector< vector<double> > SliceDiscriminators;
	vector< vector< vector<double> > > SliceBinning;

	//----------------------------------------//		

	// 2D analysis

	PlotNames.push_back("SerialDeltaPT_MuonCosThetaPlot"); 
	PlotNames.push_back("SerialDeltaPT_ProtonCosThetaPlot");
	PlotNames.push_back("SerialDeltaPT_DeltaAlphaTPlot");	
	PlotNames.push_back("SerialMuonMomentum_MuonCosThetaPlot");
	PlotNames.push_back("SerialProtonMomentum_ProtonCosThetaPlot");
	PlotNames.push_back("SerialDeltaAlphaT_MuonCosThetaPlot");
	PlotNames.push_back("SerialDeltaAlphaT_ProtonCosThetaPlot");
	PlotNames.push_back("SerialDeltaAlphaT_DeltaPTPlot");		
	PlotNames.push_back("SerialDeltaPhiT_DeltaPTPlot");
	PlotNames.push_back("SerialDeltaPn_DeltaPTPlot");	
	PlotNames.push_back("SerialProtonCosTheta_MuonCosThetaPlot");
	PlotNames.push_back("SerialDeltaPty_DeltaPtxPlot");	
	PlotNames.push_back("SerialDeltaPtx_DeltaPtyPlot");
	PlotNames.push_back("SerialECal_DeltaPTPlot");
	PlotNames.push_back("SerialECal_DeltaAlphaTPlot");

	//----------------------------------------//	

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	//----------------------------------------//

	vector<TString> Runs;
//	Runs.push_back("Run1");
//	Runs.push_back("Run2");	
//	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run4a");
//	Runs.push_back("Run5");
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	//----------------------------------------//

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		vector<vector<TH1D*> > PlotsTotalReco; PlotsTotalReco.clear();
		vector<vector<TH1D*> > PlotsNormOnly; PlotsNormOnly.clear();		
		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
		vector<vector<TH1D*> > QEPlotsTrue; QEPlotsTrue.clear();
		vector<vector<TH1D*> > MECPlotsTrue; MECPlotsTrue.clear();
		vector<vector<TH1D*> > RESPlotsTrue; RESPlotsTrue.clear();
		vector<vector<TH1D*> > DISPlotsTrue; DISPlotsTrue.clear();	
		vector<vector<TH1D*> > COHPlotsTrue; COHPlotsTrue.clear();		

		gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();

		vector<TString> NameOfSamples; NameOfSamples.clear();
		vector<int> Colors; Colors.clear();		
		vector<TString> Labels; Labels.clear();

		// CV

		NameOfSamples.push_back("Overlay9"); Colors.push_back(OverlayColor); Labels.push_back("GENIE v3 G18 Tune"); //Labels.push_back("MC uB Tune");                     

		//----------------------------------------//

		if (PlotGENIE) {

			NameOfSamples.push_back("GENIEv2");	Colors.push_back(GENIEv2Color); Labels.push_back("GENIE v2");
			NameOfSamples.push_back("Genie_v3_0_6_Out_Of_The_Box");	Colors.push_back(Geniev3OutOfTheBoxColor); Labels.push_back("GENIE v3.0.6 G18");					
			NameOfSamples.push_back("SuSav2"); Colors.push_back(SuSav2Color); Labels.push_back("GENIE v3 G21");

		}

		//----------------------------------------//	

		if (PlotGen) {

			NameOfSamples.push_back("Overlay9NuWro"); Colors.push_back(NuWroColor); Labels.push_back("NuWro 19.02.1");			
			NameOfSamples.push_back("GiBUU"); Colors.push_back(GiBUUColor); Labels.push_back("GiBUU 2021");
			NameOfSamples.push_back("NEUT"); Colors.push_back(NEUTColor); Labels.push_back("NEUT v5.4.0");

		}	

		//----------------------------------------//			

		if (PlotGENIECT) {

			NameOfSamples.push_back("Genie_v3_0_6_Nominal"); Colors.push_back(NEUTColor); Labels.push_back("GENIE v3 G18");

		}

		//----------------------------------------//

		if (PlotGENIEFSITweaks) {

			NameOfSamples.push_back("Genie_v3_0_6_NoFSI"); Colors.push_back(GiBUUColor); Labels.push_back("GENIE v3 G18 No FSI Tune");			
			NameOfSamples.push_back("Genie_v3_0_6_hN2018"); Colors.push_back(GENIEv2Color); Labels.push_back("GENIE v3 G18 hN Tune");

		}

		//----------------------------------------//

		if (PlotGENIEFlagTweaks) {

			NameOfSamples.push_back("Genie_v3_0_6_NoRPA"); Colors.push_back(NuWroColor); Labels.push_back("GENIE v3 G18 No RPA Tune");
			NameOfSamples.push_back("Genie_v3_0_6_NoCoulomb"); Colors.push_back(GENIEv3_0_4_Color); Labels.push_back("GENIE v3 G18 No Coulomb Tune");

		}

		//----------------------------------------//

		if (PlotNuclModels) {

			NameOfSamples.push_back("Genie_v3_0_6_RFG"); Colors.push_back(GiBUUColor); Labels.push_back("GENIE v3 G18 RFG Tune");			

		}  

		//----------------------------------------//

		if (All) {

			NameOfSamples.push_back("GENIEv2");	Colors.push_back(GENIEv2Color); Labels.push_back("GENIE v2");
			NameOfSamples.push_back("Genie_v3_0_6_Out_Of_The_Box");	Colors.push_back(Geniev3OutOfTheBoxColor); Labels.push_back("GENIE v3.0.6 G18");					
			NameOfSamples.push_back("SuSav2"); Colors.push_back(SuSav2Color); Labels.push_back("GENIE v3 G21");
			NameOfSamples.push_back("Overlay9NuWro"); Colors.push_back(NuWroColor); Labels.push_back("NuWro 19.02.1");			
			NameOfSamples.push_back("GiBUU"); Colors.push_back(GiBUUColor); Labels.push_back("GiBUU 2021");
			NameOfSamples.push_back("NEUT"); Colors.push_back(NEUTColor); Labels.push_back("NEUT v5.4.0");
			NameOfSamples.push_back("Genie_v3_0_6_Nominal"); Colors.push_back(NEUTColor); Labels.push_back("GENIE v3 G18");
			//NameOfSamples.push_back("Genie_v3_0_6_NoFSI"); Colors.push_back(GiBUUColor); Labels.push_back("GENIE v3 G18 No FSI Tune");			
			NameOfSamples.push_back("Genie_v3_0_6_hN2018"); Colors.push_back(GENIEv2Color); Labels.push_back("GENIE v3 G18 hN Tune");
			//NameOfSamples.push_back("Genie_v3_0_6_NoRPA"); Colors.push_back(NuWroColor); Labels.push_back("GENIE v3 G18 No RPA Tune");
			//NameOfSamples.push_back("Genie_v3_0_6_NoCoulomb"); Colors.push_back(GENIEv3_0_4_Color); Labels.push_back("GENIE v3 G18 No Coulomb Tune");			
			//NameOfSamples.push_back("Genie_v3_0_6_RFG"); Colors.push_back(GiBUUColor); Labels.push_back("GENIE v3 G18 RFG Tune");			

		}			             

		//----------------------------------------//

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();

		//----------------------------------------//

		// Unfolding uncertainty

		TFile* fUnc = TFile::Open(PathToFiles+UBCodeVersion+"/WienerSVD_UnfoldingUnc_Combined_"+UBCodeVersion+".root","readonly");

		//----------------------------------------//

		// File to store all the generator predictions

		TFile* fGenXSec = nullptr;
		TString GenXSecName = PathToFiles+UBCodeVersion+"/GenXSec/All_XSecs_2D_" + Runs[WhichRun] + "_"+UBCodeVersion+".root";

		// If All == true	

		if (All) {

			fGenXSec = TFile::Open(GenXSecName,"recreate");

		}			

		//----------------------------------------//

		// Open the files and grap the relevant plots

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			vector<TH1D*> CurrentPlotsTotalReco; CurrentPlotsTotalReco.clear();
			vector<TH1D*> CurrentPlotsNormOnly; CurrentPlotsNormOnly.clear();			
			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();
			vector<TH1D*> CurrentPlotsCC1pReco; CurrentPlotsCC1pReco.clear();
			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();
			vector<TH1D*> QECurrentPlotsTrue; QECurrentPlotsTrue.clear();
			vector<TH1D*> MECCurrentPlotsTrue; MECCurrentPlotsTrue.clear();	
			vector<TH1D*> RESCurrentPlotsTrue; RESCurrentPlotsTrue.clear();
			vector<TH1D*> DISCurrentPlotsTrue; DISCurrentPlotsTrue.clear();	
			vector<TH1D*> COHCurrentPlotsTrue; COHCurrentPlotsTrue.clear();			

			// CV With Statistical Uncertainties

			if (NameOfSamples[WhichSample] == "Overlay9") { // CV with statistical uncertainties only for now

				TString FileSampleName = PathToFiles+UBCodeVersion+"/WienerSVD_ExtractedXSec_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; 
				FileSample.push_back(TFile::Open(FileSampleName,"readonly")); 

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TH1D* histTotalReco = (TH1D*)(FileSample[WhichSample]->Get("StatReco"+PlotNames[WhichPlot]));
					CurrentPlotsTotalReco.push_back(histTotalReco);

					TH1D* histNormOnly = (TH1D*)(FileSample[WhichSample]->Get("NormOnlyReco"+PlotNames[WhichPlot]));
					CurrentPlotsNormOnly.push_back(histNormOnly);					

					TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
					if (PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot") { histReco = (TH1D*)(FileSample[WhichSample]->Get("RecoFullUnc"+PlotNames[WhichPlot])); }
					CurrentPlotsReco.push_back(histReco);

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);

					TH1D* QEhistTrue = (TH1D*)(FileSample[WhichSample]->Get("QETrue"+PlotNames[WhichPlot]));
					QECurrentPlotsTrue.push_back(QEhistTrue);
					TH1D* MEChistTrue = (TH1D*)(FileSample[WhichSample]->Get("MECTrue"+PlotNames[WhichPlot]));
					MECCurrentPlotsTrue.push_back(MEChistTrue);
					TH1D* REShistTrue = (TH1D*)(FileSample[WhichSample]->Get("RESTrue"+PlotNames[WhichPlot]));
					RESCurrentPlotsTrue.push_back(REShistTrue);
					TH1D* DIShistTrue = (TH1D*)(FileSample[WhichSample]->Get("DISTrue"+PlotNames[WhichPlot]));
					DISCurrentPlotsTrue.push_back(DIShistTrue);
					TH1D* COHhistTrue = (TH1D*)(FileSample[WhichSample]->Get("COHTrue"+PlotNames[WhichPlot]));
					COHCurrentPlotsTrue.push_back(COHhistTrue);					
		
				}

			} else if (NameOfSamples[WhichSample] == "Overlay9NuWro") {

				TString FileSampleName = PathToFiles+UBCodeVersion+"/"+NameOfSamples[WhichSample]+"WienerSVD_ExtractedXSec_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; 
				FileSample.push_back(TFile::Open(FileSampleName,"readonly")); 

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					//TH1D* histTotalReco = (TH1D*)(FileSample[WhichSample]->Get("StatReco"+PlotNames[WhichPlot]));
					CurrentPlotsTotalReco.push_back(nullptr);
					CurrentPlotsNormOnly.push_back(nullptr);					

					TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
					CurrentPlotsReco.push_back(histReco);

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("NoSmearAltTrue"+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);

					TH1D* QEhistTrue = (TH1D*)(FileSample[WhichSample]->Get("QENoSmearAltTrue"+PlotNames[WhichPlot]));
					QECurrentPlotsTrue.push_back(QEhistTrue);
					TH1D* MEChistTrue = (TH1D*)(FileSample[WhichSample]->Get("MECNoSmearAltTrue"+PlotNames[WhichPlot]));
					MECCurrentPlotsTrue.push_back(MEChistTrue);
					TH1D* REShistTrue = (TH1D*)(FileSample[WhichSample]->Get("RESNoSmearAltTrue"+PlotNames[WhichPlot]));
					RESCurrentPlotsTrue.push_back(REShistTrue);
					TH1D* DIShistTrue = (TH1D*)(FileSample[WhichSample]->Get("DISNoSmearAltTrue"+PlotNames[WhichPlot]));
					DISCurrentPlotsTrue.push_back(DIShistTrue);
					TH1D* COHhistTrue = (TH1D*)(FileSample[WhichSample]->Get("COHNoSmearAltTrue"+PlotNames[WhichPlot]));
					COHCurrentPlotsTrue.push_back(COHhistTrue);					
		
				}

			}

			else {

				if (
					NameOfSamples[WhichSample] == "Genie_v3_0_6_Out_Of_The_Box" || 
					NameOfSamples[WhichSample] == "Genie_v3_0_6_uB_Tune_1" || 
					NameOfSamples[WhichSample] == "Genie_v3_0_6_Nominal" || 
					NameOfSamples[WhichSample] == "Genie_v3_0_6_NoFSI" || 
					NameOfSamples[WhichSample] == "Genie_v3_0_6_NoRPA" || 
					NameOfSamples[WhichSample] == "Genie_v3_0_6_NoCoulomb" || 
					NameOfSamples[WhichSample] == "Genie_v3_0_6_hN2018" ||
					NameOfSamples[WhichSample] == "Genie_v3_0_6_RFG" ||  
					NameOfSamples[WhichSample] == "Genie_v3_0_6_EffSF" ||  
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

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TH1D* histTotalReco = nullptr;
					CurrentPlotsTotalReco.push_back(histTotalReco);

					TH1D* histNormOnly = nullptr;
					CurrentPlotsNormOnly.push_back(histNormOnly);					

					TH1D* histReco = nullptr;
					CurrentPlotsReco.push_back(histReco);

					TH1D* histCC1pReco = nullptr;
					CurrentPlotsCC1pReco.push_back(histCC1pReco);

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);

					TH1D* QEhistTrue = (TH1D*)(FileSample[WhichSample]->Get("QETrue"+PlotNames[WhichPlot]));
					QECurrentPlotsTrue.push_back(QEhistTrue);
					TH1D* MEChistTrue = (TH1D*)(FileSample[WhichSample]->Get("MECTrue"+PlotNames[WhichPlot]));
					MECCurrentPlotsTrue.push_back(MEChistTrue);
					TH1D* REShistTrue = (TH1D*)(FileSample[WhichSample]->Get("RESTrue"+PlotNames[WhichPlot]));
					RESCurrentPlotsTrue.push_back(REShistTrue);
					TH1D* DIShistTrue = (TH1D*)(FileSample[WhichSample]->Get("DISTrue"+PlotNames[WhichPlot]));
					DISCurrentPlotsTrue.push_back(DIShistTrue);
					TH1D* COHhistTrue = (TH1D*)(FileSample[WhichSample]->Get("COHTrue"+PlotNames[WhichPlot]));
					COHCurrentPlotsTrue.push_back(COHhistTrue);					
		
				}

			}

			PlotsTotalReco.push_back(CurrentPlotsTotalReco);
			PlotsNormOnly.push_back(CurrentPlotsNormOnly);					
			PlotsReco.push_back(CurrentPlotsReco);		
			PlotsCC1pReco.push_back(CurrentPlotsCC1pReco);
			PlotsTrue.push_back(CurrentPlotsTrue);	

			QEPlotsTrue.push_back(QECurrentPlotsTrue);			
			MECPlotsTrue.push_back(MECCurrentPlotsTrue);
			RESPlotsTrue.push_back(RESCurrentPlotsTrue);
			DISPlotsTrue.push_back(DISCurrentPlotsTrue);
			COHPlotsTrue.push_back(COHCurrentPlotsTrue);				

		}

		//----------------------------------------//

		// Loop over the plots

		vector< vector<TH1D*> > BeamOnStatShape;
		vector< vector<TH1D*> > BeamOnStatOnly;
		vector< vector<TH1D*> > BeamOnNormOnly;
		vector< vector< vector<TH1D*> > > MC;
		vector< vector< vector<TH1D*> > > QEMC;
		vector< vector< vector<TH1D*> > > MECMC;
		vector< vector< vector<TH1D*> > > RESMC;
		vector< vector< vector<TH1D*> > > DISMC;
		vector< vector< vector<TH1D*> > > COHMC;													

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {		

			//----------------------------------------//

			// Setting up the relevant discriminators

			SliceDiscriminators.clear();
			SliceBinning.clear();

			if (PlotNames[WhichPlot] == "SerialDeltaPT_MuonCosThetaPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsMuonCosTheta); 
				SliceBinning.push_back(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices);

			}

			if (PlotNames[WhichPlot] == "SerialDeltaPT_ProtonCosThetaPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsProtonCosTheta); 
				SliceBinning.push_back(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices);

			}

			if (PlotNames[WhichPlot] == "SerialDeltaPT_DeltaAlphaTPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsDeltaAlphaT); 
				SliceBinning.push_back(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices);

			}			

			if (PlotNames[WhichPlot] == "SerialMuonMomentum_MuonCosThetaPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsMuonCosTheta); 
				SliceBinning.push_back(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices);

			}	

			if (PlotNames[WhichPlot] == "SerialProtonMomentum_ProtonCosThetaPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsProtonCosTheta); 
				SliceBinning.push_back(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices);

			}	

			if (PlotNames[WhichPlot] == "SerialDeltaAlphaT_MuonCosThetaPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsMuonCosTheta); 
				SliceBinning.push_back(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices);

			}	

			if (PlotNames[WhichPlot] == "SerialDeltaAlphaT_ProtonCosThetaPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsProtonCosTheta); 
				SliceBinning.push_back(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices);

			}	

			if (PlotNames[WhichPlot] == "SerialDeltaAlphaT_DeltaPTPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsDeltaPT); 
				SliceBinning.push_back(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices);

			}	

			if (PlotNames[WhichPlot] == "SerialDeltaPhiT_DeltaPTPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsDeltaPT); 
				SliceBinning.push_back(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices);

			}	

			if (PlotNames[WhichPlot] == "SerialDeltaPn_DeltaPTPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsDeltaPT); 
				SliceBinning.push_back(TwoDArrayNBinsDeltaPnInDeltaPTSlices);

			}	

			if (PlotNames[WhichPlot] == "SerialECal_DeltaPTPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsDeltaPT); 
				SliceBinning.push_back(TwoDArrayNBinsECalInDeltaPTSlices);

			}	

			if (PlotNames[WhichPlot] == "SerialECal_DeltaAlphaTPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsDeltaAlphaT); 
				SliceBinning.push_back(TwoDArrayNBinsECalInDeltaAlphaTSlices);

			}				

			if (PlotNames[WhichPlot] == "SerialProtonCosTheta_MuonCosThetaPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsMuonCosTheta); 
				SliceBinning.push_back(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices);

			}	

			if (PlotNames[WhichPlot] == "SerialDeltaPtx_DeltaPtyPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsDeltaPty); 
				SliceBinning.push_back(TwoDArrayNBinsDeltaPtxInDeltaPtySlices);

			}	

			if (PlotNames[WhichPlot] == "SerialDeltaPty_DeltaPtxPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsDeltaPtx); 
				SliceBinning.push_back(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices);

			}	

			//----------------------------------------//

			BeamOnStatShape.resize(N1DPlots);
			BeamOnStatOnly.resize(N1DPlots);
			BeamOnNormOnly.resize(N1DPlots);							
			MC.resize(N1DPlots);
			QEMC.resize(N1DPlots);
			MECMC.resize(N1DPlots);
			RESMC.resize(N1DPlots);
			DISMC.resize(N1DPlots);
			COHMC.resize(N1DPlots);															

			//----------------------------------------//

			TH2D* Ac = (TH2D*)FileSample[0]->Get("Ac"+PlotNames[WhichPlot]);
			TH2D* Cov = (TH2D*)FileSample[0]->Get("UnfCov"+PlotNames[WhichPlot]);	

			//----------------------------------------//

			TH1D* UncHist = (TH1D*)(fUnc->Get("UnfUnc_"+PlotNames[WhichPlot]));

			//----------------------------------------//

			// The covariance matrix needs to be scaled by the 2D bin width

			TH2D* CovClone = (TH2D*)(Cov->Clone()); 

			int n = Cov->GetXaxis()->GetNbins();

			for (int ix = 1; ix <= n; ix++) {

				for (int iy = 1; iy <= n; iy++) {

					double WidthX = Cov->GetXaxis()->GetBinWidth(ix);
					double WidthY = Cov->GetYaxis()->GetBinWidth(iy);

					double TwoDWidth = WidthX * WidthY;

					double BinContent = Cov->GetBinContent(ix,iy);
					// Division by bin width already included
					// Still need to include scaling due to slice range
					// That is done in Tools::Get2DHistoBins
					double NewBinContent = BinContent/TwoDWidth;

					// Only for the diagonal elements
					// Add the unfolding uncertainty
					// On top of everything else
					// That is done both for the final xsec result and for the unfolded covariance
					if (ix == iy) { 
						
						// unfolded covariance matrix
						double UnfUncBin = UncHist->GetBinContent(ix);
//						double UnfUncBin = 0.;
						NewBinContent = NewBinContent + TMath::Power(UnfUncBin,2.) ; 

						// xsec uncertainty
						double CurrentUnc = PlotsReco[0][WhichPlot]->GetBinError(ix);
						double NewError = TMath::Sqrt( TMath::Power(CurrentUnc,2.) + TMath::Power(UnfUncBin,2.) ) ;
						PlotsReco[0][WhichPlot]->SetBinError(ix,NewError);
						
					}

					CovClone->SetBinContent(ix,iy,NewBinContent);

				}					

			}	

			//CovClone->Draw("coltz text");

			//------------------------------------//

			// Number of N-dimensional slices

			int NSlices = 1;
			vector<double> SerialVectorRanges;
			vector<int> SerialVectorBins;
			vector<int> SerialVectorLowBin;
			vector<int> SerialVectorHighBin;						

			int BinCounter = 1;

			// vector< vector<double> > SliceDiscriminators;
			// 1st index = how many disciminators
			// 2nd index = values / ranges of discriminators

			// How many discriminators do we have ? e.g. cos theta mu, Pmu et al
			for (int islice = 0; islice < (int)(SliceDiscriminators.size()); islice++) { 
				
				// For a given discriminator, how many slices do we have ? SliceDiscrimSize - 1
				int SliceDiscrimSize = SliceDiscriminators.at(islice).size()-1;
				NSlices *= SliceDiscrimSize; 

				for (int iSliceDiscrimSize = 0; iSliceDiscrimSize < SliceDiscrimSize; iSliceDiscrimSize++) {

					// Accessing the vector<double> with the bin ranges
					int SliceDiscrimValue = SliceBinning.at(0).at(iSliceDiscrimSize).size();

					// Storing the number of bins for a specific slice					
					SerialVectorBins.push_back(SliceDiscrimValue-1);
					for (int iBin = 0; iBin < SliceDiscrimValue; iBin++) {

						double BinValue = SliceBinning.at(0).at(iSliceDiscrimSize).at(iBin);
						// First bin number for a given slice
						if (iBin == 0) { SerialVectorLowBin.push_back(BinCounter); }
						// Last bin number for a given slice
						if (iBin == SliceDiscrimValue-2) { SerialVectorHighBin.push_back(BinCounter); }	

						// Storing the binning for a specific slice
						SerialVectorRanges.push_back(BinValue);

						// Increase the global bin counter
						// But not for the last bin
						if (iBin != SliceDiscrimValue-1) { BinCounter++; }

					} // End of the loop over the bins of a given slice
					//cout << "End of the loop over the bins of a given slice" << endl;

				} // End of the loop over the slices for a given discriminator

			} // End of the loop over the number of discriminators	

			//------------------------------------//

			BeamOnStatShape[WhichPlot].resize(NSlices);
			BeamOnStatOnly[WhichPlot].resize(NSlices);
			BeamOnNormOnly[WhichPlot].resize(NSlices);
			MC[WhichPlot].resize(NSlices);
			QEMC[WhichPlot].resize(NSlices);
			MECMC[WhichPlot].resize(NSlices);
			RESMC[WhichPlot].resize(NSlices);
			DISMC[WhichPlot].resize(NSlices);	
			COHMC[WhichPlot].resize(NSlices);																	

			//------------------------------------//

			int StartIndex = 0;
			int BinStartIndex = 0;			

			//------------------------------------//

			// MC multiplication by the additional matrix Ac
			// NOT for index 0 = overlay, that has already been done

			for (int WhichSample = 1; WhichSample < NSamples; WhichSample++) {

				PlotsTrue[WhichSample][WhichPlot] = Multiply(PlotsTrue[WhichSample][WhichPlot],Ac);
				QEPlotsTrue[WhichSample][WhichPlot] = Multiply(QEPlotsTrue[WhichSample][WhichPlot],Ac);
				MECPlotsTrue[WhichSample][WhichPlot] = Multiply(MECPlotsTrue[WhichSample][WhichPlot],Ac);
				RESPlotsTrue[WhichSample][WhichPlot] = Multiply(RESPlotsTrue[WhichSample][WhichPlot],Ac);
				DISPlotsTrue[WhichSample][WhichPlot] = Multiply(DISPlotsTrue[WhichSample][WhichPlot],Ac);
				COHPlotsTrue[WhichSample][WhichPlot] = Multiply(COHPlotsTrue[WhichSample][WhichPlot],Ac);																					

			}							

			//------------------------------------//

			// Loop over the N-dimensional slices

			for (int NDimSlice = 0; NDimSlice < NSlices; NDimSlice++) {	

				//------------------------------------//

				MC[WhichPlot][NDimSlice].resize(NSamples);
				QEMC[WhichPlot][NDimSlice].resize(NSamples);
				MECMC[WhichPlot][NDimSlice].resize(NSamples);	
				RESMC[WhichPlot][NDimSlice].resize(NSamples);
				DISMC[WhichPlot][NDimSlice].resize(NSamples);	
				COHMC[WhichPlot][NDimSlice].resize(NSamples);																		

				//------------------------------------//

				TString NameCopy = PlotNames[WhichPlot];

				NameCopy.ReplaceAll("unf_","");
				NameCopy.ReplaceAll("TrueUnf_","");
				NameCopy.ReplaceAll("True_","");
				NameCopy.ReplaceAll("True","");
				NameCopy.ReplaceAll("NoSmearAlt","");			

				NameCopy.ReplaceAll("_Run1","");
				NameCopy.ReplaceAll("_Run2","");
				NameCopy.ReplaceAll("_Run3","");
				NameCopy.ReplaceAll("_Run4","");
				NameCopy.ReplaceAll("_Run4a","");	
				NameCopy.ReplaceAll("_Run5","");
				NameCopy.ReplaceAll("_Combined","");	

				NameCopy = NameCopy + "_" + TString(std::to_string(NDimSlice));	

				//------------------------------------//		

				// Get the number of bins and the bin ranges for the specific slice	

				int SliceNBins = SerialVectorBins.at(NDimSlice);
				std::vector<double> SerialSliceBinning;		

				for (int iBin = 0; iBin < SliceNBins+1; iBin++) { 

					double value = SerialVectorRanges.at(StartIndex+iBin);
					SerialSliceBinning.push_back(value);

				} // End of the number of bins and the bin ranges declaration	

				//------------------------------------//

				// Canvas, pads & legend				

				TString CanvasName = PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_Slice_"+TString(std::to_string(NDimSlice));
				TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
				PlotCanvas->cd();
				PlotCanvas->SetBottomMargin(0.14);
				PlotCanvas->SetTopMargin(0.12);
				PlotCanvas->SetLeftMargin(0.18);

				TLegend* leg = new TLegend(0.6,0.58,0.7,0.85);
				TLegend* legChi2 = new TLegend(0.78,0.72,0.88,0.85);			
				if (
					PlotNames[WhichPlot] == "SerialDeltaAlphaT_DeltaPTPlot" ||																				
					PlotNames[WhichPlot] == "SerialDeltaPtx_DeltaPtyPlot" ||
					PlotNames[WhichPlot] == "SerialDeltaPty_DeltaPtxPlot" ||
					PlotNames[WhichPlot] == "SerialProtonCosTheta_MuonCosThetaPlot" ||
					PlotNames[WhichPlot] == "SerialDeltaAlphaT_MuonCosThetaPlot" ||
					(PlotNames[WhichPlot] == "SerialMuonMomentum_MuonCosThetaPlot" && NDimSlice == 3) ||
					(PlotNames[WhichPlot] == "SerialDeltaPn_DeltaPTPlot" && NDimSlice == 2) ||										
					PlotNames[WhichPlot] == "SerialDeltaAlphaT_ProtonCosThetaPlot"
					) { 
						
						leg = new TLegend(0.22,0.58,0.32,0.85); 
						legChi2 = new TLegend(0.4,0.72,0.5,0.85);					
						
				}

				if (PlotNominal) { leg = new TLegend(0.6,0.68,0.71,0.85); }
				if (PlotNominal && PlotNames[WhichPlot] == "DeltaPtyPlot") { leg = new TLegend(0.2,0.58,0.5,0.85); }			

				leg->SetBorderSize(0);
				leg->SetTextSize(0.03);
				leg->SetTextFont(FontStyle);
				leg->SetNColumns(1);
				leg->SetMargin(0.15);

				legChi2->SetBorderSize(0);
				legChi2->SetTextSize(0.03);
				legChi2->SetTextFont(FontStyle);
				legChi2->SetNColumns(1);
				legChi2->SetMargin(0.15);

				//------------------------------------//

				// Corresponding covariance matrix

				TH2D* SliceCovMatrix = tools.Get2DHistoBins(CovClone,SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning);				

				//------------------------------------//

				// BeamOn Total Uncertainty																				

				BeamOnStatShape[WhichPlot][NDimSlice] = tools.GetHistoBins(PlotsReco[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"StatShape");										

				PrettyPlot(BeamOnStatShape[WhichPlot][NDimSlice]);

				double MaxValue = BeamOnStatShape[WhichPlot][NDimSlice]->GetMaximum();
				int MaxValueBin = LocateBinWithValue(BeamOnStatShape[WhichPlot][NDimSlice],MaxValue);
				double MaxValueError = BeamOnStatShape[WhichPlot][NDimSlice]->GetBinError(MaxValueBin);

				double MinValue = BeamOnStatShape[WhichPlot][NDimSlice]->GetMinimum();
														
				BeamOnStatShape[WhichPlot][NDimSlice]->GetYaxis()->SetRangeUser(XSecRange[ MapUncorCor[ NameCopy ] ].first,XSecRange[ MapUncorCor[ NameCopy ] ].second);

				BeamOnStatShape[WhichPlot][NDimSlice]->SetLineColor(BeamOnColor);
				BeamOnStatShape[WhichPlot][NDimSlice]->SetMarkerColor(BeamOnColor);
				BeamOnStatShape[WhichPlot][NDimSlice]->SetMarkerSize(1.);
				BeamOnStatShape[WhichPlot][NDimSlice]->SetMarkerStyle(20);
				BeamOnStatShape[WhichPlot][NDimSlice]->SetLineWidth(1);	
				BeamOnStatShape[WhichPlot][NDimSlice]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);							

				BeamOnStatShape[WhichPlot][NDimSlice]->Draw("e1x0 same"); // Total Unc (Shape + Stat)

				//------------------------------//

				// arrays for MC NSamples

				double Chi2[NSamples];			
				int Ndof[NSamples];
				double pval[NSamples];

				//------------------------------//									
		
				// Alternative MC

				for (int WhichSample = 1; WhichSample < NSamples; WhichSample++) {

					MC[WhichPlot][NDimSlice][WhichSample] = tools.GetHistoBins(PlotsTrue[WhichSample][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning, NameOfSamples[WhichSample]);
					MC[WhichPlot][NDimSlice][WhichSample]->SetLineColor(Colors[WhichSample]);
					MC[WhichPlot][NDimSlice][WhichSample]->SetMarkerColor(Colors[WhichSample]);
					MC[WhichPlot][NDimSlice][WhichSample]->Draw("hist same");

					QEMC[WhichPlot][NDimSlice][WhichSample] = tools.GetHistoBins(QEPlotsTrue[WhichSample][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning, NameOfSamples[WhichSample]);
					MECMC[WhichPlot][NDimSlice][WhichSample] = tools.GetHistoBins(MECPlotsTrue[WhichSample][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning, NameOfSamples[WhichSample]);										
					RESMC[WhichPlot][NDimSlice][WhichSample] = tools.GetHistoBins(RESPlotsTrue[WhichSample][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning, NameOfSamples[WhichSample]);
					DISMC[WhichPlot][NDimSlice][WhichSample] = tools.GetHistoBins(DISPlotsTrue[WhichSample][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning, NameOfSamples[WhichSample]);
					COHMC[WhichPlot][NDimSlice][WhichSample] = tools.GetHistoBins(COHPlotsTrue[WhichSample][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning, NameOfSamples[WhichSample]);

					CalcChiSquared(MC[WhichPlot][NDimSlice][WhichSample],BeamOnStatShape[WhichPlot][NDimSlice],SliceCovMatrix,Chi2[WhichSample],Ndof[WhichSample],pval[WhichSample]);
					TString Chi2NdofAlt = " (" + to_string_with_precision(Chi2[WhichSample],2) + "/" + TString(std::to_string(Ndof[WhichSample])) +")";
					TLegendEntry* lGenie = leg->AddEntry(MC[WhichPlot][NDimSlice][WhichSample],Labels[WhichSample],"l");
					lGenie->SetTextColor(Colors[WhichSample]); 										
					TLegendEntry* lGenieChi2 = legChi2->AddEntry(MC[WhichPlot][NDimSlice][WhichSample],Chi2NdofAlt,"");
					lGenieChi2->SetTextColor(Colors[WhichSample]);					

				}		

				//------------------------------//

				// Overlay

				MC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(PlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay");
				PrettyPlot(MC[WhichPlot][NDimSlice][0]);
				MC[WhichPlot][NDimSlice][0]->SetLineColor(Colors[0]);
				MC[WhichPlot][NDimSlice][0]->SetMarkerColor(Colors[0]);	
				MC[WhichPlot][NDimSlice][0]->Draw("hist same");	

				QEMC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(QEPlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay");
				MECMC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(MECPlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay");								
				RESMC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(RESPlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay");
				DISMC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(DISPlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay");
				COHMC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(COHPlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay");

				CalcChiSquared(MC[WhichPlot][NDimSlice][0],BeamOnStatShape[WhichPlot][NDimSlice],SliceCovMatrix,Chi2[0],Ndof[0],pval[0]);
				TString Chi2NdofAlt = " (" + to_string_with_precision(Chi2[0],2) + "/" + TString(std::to_string(Ndof[0])) +")";
				TLegendEntry* lGenie = leg->AddEntry(MC[WhichPlot][NDimSlice][0],Labels[0],"l");
				lGenie->SetTextColor(Colors[0]); 										
				TLegendEntry* lGenieChi2 = legChi2->AddEntry(MC[WhichPlot][NDimSlice][0],Chi2NdofAlt,"");
				lGenieChi2->SetTextColor(Colors[0]);					

				//------------------------------//	

				// Stat Unc Only	

				BeamOnStatOnly[WhichPlot][NDimSlice] = tools.GetHistoBins(PlotsTotalReco[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"StatOnly");											
				PrettyPlot(BeamOnStatOnly[WhichPlot][NDimSlice]);
				BeamOnStatOnly[WhichPlot][NDimSlice]->SetLineColor(BeamOnColor);
				BeamOnStatOnly[WhichPlot][NDimSlice]->SetMarkerColor(BeamOnColor);
				BeamOnStatOnly[WhichPlot][NDimSlice]->SetLineWidth(1);			
				BeamOnStatOnly[WhichPlot][NDimSlice]->Draw("e1x0 same"); // Stat Only

				// Plot again on top
				BeamOnStatShape[WhichPlot][NDimSlice]->Draw("e1x0 same"); // Total Unc (Shape + Stat)				
				
				//------------------------------//

				// Norm Unc Only

				BeamOnNormOnly[WhichPlot][NDimSlice] = tools.GetHistoBins(PlotsNormOnly[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"NormOnly");
				PrettyPlot(BeamOnNormOnly[WhichPlot][NDimSlice]); // includes scaling factor for multi dimensional analysis			
				BeamOnNormOnly[WhichPlot][NDimSlice]->SetFillColorAlpha(kGray+1, 0.45);	
				BeamOnNormOnly[WhichPlot][NDimSlice]->SetLineColor(kGray+1);
				BeamOnNormOnly[WhichPlot][NDimSlice]->SetMarkerColor(kGray+1);
				BeamOnNormOnly[WhichPlot][NDimSlice]->Draw("e2 same");				

				// -----------------------------------------------------------------------------------------------------------------			

				// Legend & Run / POT

				double tor860_wcut = -99.;
				if (Runs[WhichRun] == "Run1") { tor860_wcut = Fulltor860_wcut_Run1; }
				if (Runs[WhichRun] == "Run2") { tor860_wcut = Fulltor860_wcut_Run2; }
				if (Runs[WhichRun] == "Run3") { tor860_wcut = Fulltor860_wcut_Run3; }
				if (Runs[WhichRun] == "Run4") { tor860_wcut = Fulltor860_wcut_Run4; }
				if (Runs[WhichRun] == "Run4a") { tor860_wcut = Fulltor860_wcut_Run4a; }				
				if (Runs[WhichRun] == "Run5") { tor860_wcut = Fulltor860_wcut_Run5; }
				if (Runs[WhichRun] == "Combined") { tor860_wcut = Fulltor860_wcut_Combined; }
				TString Label = ToString(tor860_wcut)+" POT";			

				// ---------------------------------------------------------------------------------------------------------
				// ---------------------------------------------------------------------------------------------------------

				leg->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],"MicroBooNE Data","ep");
				leg->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],"(Stat #oplus Shape Unc)","");
				leg->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],Label,"");
				leg->AddEntry(BeamOnNormOnly[WhichPlot][NDimSlice],"Norm Unc","f");

				legChi2->Draw();
				leg->Draw();			

				TLatex *textSlice = new TLatex();
				textSlice->SetTextFont(FontStyle);
				textSlice->SetTextSize(0.06);
				TString PlotNameDuplicate = NameCopy;
				TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("Reco","") ;
				textSlice->DrawLatexNDC(0.24, 0.92, LatexLabel[ MapUncorCor[ReducedPlotName] ]);	

				//----------------------------------------//

				// If the option All == true is activated, store all the relevant xsecs in a file

				if (All) {

					fGenXSec->cd();

					// Unfolded covariance matrix
					SliceCovMatrix->Write("UnfCov_" + NameCopy);					

					// Data
					BeamOnStatShape[WhichPlot][NDimSlice]->Write("StatShape_" + NameCopy);
					BeamOnNormOnly[WhichPlot][NDimSlice]->Write("NormOnly_" + NameCopy);
					BeamOnStatOnly[WhichPlot][NDimSlice]->Write("StatOnly_" + NameCopy);

					// Overlay GENIE
					MC[WhichPlot][NDimSlice][0]->Write("OverlayGENIE_" + NameCopy);
					QEMC[WhichPlot][NDimSlice][0]->Write("QEOverlayGENIE_" + NameCopy);
					MECMC[WhichPlot][NDimSlice][0]->Write("MECOverlayGENIE_" + NameCopy);
					RESMC[WhichPlot][NDimSlice][0]->Write("RESOverlayGENIE_" + NameCopy);
					DISMC[WhichPlot][NDimSlice][0]->Write("DISOverlayGENIE_" + NameCopy);
					COHMC[WhichPlot][NDimSlice][0]->Write("COHOverlayGENIE_" + NameCopy);																									

					// Store the remaining generator xsecs
					// Start from 1, the GENIE Overlay = 0 has been stored above

					for (int igen = 1; igen < NSamples; igen++ ) {

						TString GenName = NameOfSamples[igen] + "_" + NameCopy;
						MC[WhichPlot][NDimSlice][igen]->Write(GenName);
						QEMC[WhichPlot][NDimSlice][igen]->Write("QE" + GenName);
						MECMC[WhichPlot][NDimSlice][igen]->Write("MEC" + GenName);
						RESMC[WhichPlot][NDimSlice][igen]->Write("RES" + GenName);
						DISMC[WhichPlot][NDimSlice][igen]->Write("DIS" + GenName);
						COHMC[WhichPlot][NDimSlice][igen]->Write("COH" + GenName);																														

					}

				}							

				//------------------------------------//

				// Saving the canvas with the data (total uncertainties) vs overlay & generator predictions

				PlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/BeamOn9/"+Extra+"MultiDimWienerSVD_Generator_TotalUnc_Data_2DXSections_"+PlotNames[WhichPlot]+"_Slice_"+TString(std::to_string(NDimSlice))+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");

				delete PlotCanvas;

				//------------------------------------//

				// Update the starting index to move to the next slice

				StartIndex += (SliceNBins+1);
				BinStartIndex += SliceNBins;

				//------------------------------------//

			} // End of the loop over the discriminators slices

			//----------------------------------------//

		} // End of the loop over the plots

		//----------------------------------------//

		if (All) {

			fGenXSec->Close();
			cout << endl << GenXSecName << " file created" << endl << endl;	

		}		

		//----------------------------------------//

	} // End of the loop over the runs	

} // End of the program 
