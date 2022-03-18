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

void ThreeDimWienerSVD_OverlayGenerators(TString PlotName = "", int FirstDiscrIndex = 0, int SecondDiscrIndex = 0, 
								 		 bool PlotGENIE = true, bool PlotGen = false, 
								 		 bool PlotGENIEFSITweaks = false, bool PlotGENIEFlagTweaks = false, 
								 		 bool PlotGENIECT = false, bool PlotNuclModels = false, 
								 		 bool PlotNuWro = false, bool PlotNominal = false, 
										 bool GiBUUComp = false, bool All = false, bool NoFSIPlusGiBUU = false,
										 bool G21FSI = false
										) {

	//----------------------------------------//

	Tools tools;

	//----------------------------------------//

	int DecimalAccuracy = 2;

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(6);		

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
	if (!PlotGENIE && !PlotGen && NoFSIPlusGiBUU) { Extra = "NoFSIPlusGiBUU"; }
	if (!PlotGENIE && !PlotGen && G21FSI) { Extra = "G21FSI"; }				

	//----------------------------------------//	

	vector<TString> PlotNames;
	vector< vector<double> > SliceDiscriminators;
	vector< vector< vector<double> > > SliceBinning;
	
	//----------------------------------------//

	// 3D names

	PlotNames.push_back(PlotName);

	//----------------------------------------//	

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	//----------------------------------------//

	vector<TString> Runs;
//	Runs.push_back("Run1");
//	Runs.push_back("Run2");	
//	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run5");
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	//----------------------------------------//

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		vector<vector<TH1D*> > PlotsFullUncReco; PlotsFullUncReco.clear();
		vector<vector<TH1D*> > PlotsXSecReco; PlotsXSecReco.clear();
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

		NameOfSamples.push_back("Overlay9"); Colors.push_back(OverlayColor); Labels.push_back("G18              "); //Labels.push_back("MC uB Tune");                     

		//----------------------------------------//	

		if (PlotGENIE) {

			NameOfSamples.push_back("GENIEv2");	Colors.push_back(GENIEv2Color); Labels.push_back("Gv2       ");
			NameOfSamples.push_back("Genie_v3_0_6_Out_Of_The_Box");	Colors.push_back(Geniev3OutOfTheBoxColor); Labels.push_back("G18 No Tune");					
			NameOfSamples.push_back("SuSav2"); Colors.push_back(SuSav2Color); Labels.push_back("G21hN      ");

		}

		//----------------------------------------//		

		if (PlotGen) {

			NameOfSamples.push_back("Overlay9NuWro"); Colors.push_back(NuWroColor); Labels.push_back("NuWro           ");			
			NameOfSamples.push_back("GiBUU"); Colors.push_back(GiBUUColor); Labels.push_back("GiB              ");
			NameOfSamples.push_back("NEUT"); Colors.push_back(NEUTColor); Labels.push_back("NEUT            ");

		}	

		//----------------------------------------//			

		if (PlotGENIECT) {

			NameOfSamples.push_back("Genie_v3_0_6_Nominal"); Colors.push_back(NEUTColor); Labels.push_back("G18 Nom    ");

		}

		//----------------------------------------//

		if (PlotGENIEFSITweaks) {

			NameOfSamples.push_back("Genie_v3_0_6_NoFSI"); Colors.push_back(GiBUUColor); Labels.push_back("G18 No FSI ");			
			NameOfSamples.push_back("Genie_v3_0_6_hN2018"); Colors.push_back(GENIEv2Color); Labels.push_back("G18 hN Tune");

		}

		//----------------------------------------//

		if (PlotGENIEFlagTweaks) {

			NameOfSamples.push_back("Genie_v3_0_6_NoRPA"); Colors.push_back(NuWroColor); Labels.push_back("G18 No RPA");
			NameOfSamples.push_back("Genie_v3_0_6_NoCoulomb"); Colors.push_back(GENIEv3_0_4_Color); Labels.push_back("G18 No Coul");

		}

		//----------------------------------------//

		if (PlotNuclModels) {

//			NameOfSamples.push_back("Genie_v3_0_6_RFG"); Colors.push_back(GiBUUColor); Labels.push_back("G18 RFG Tune");
			NameOfSamples.push_back("GENIEv2");	Colors.push_back(GENIEv2Color); Labels.push_back("Gv2        ");
			NameOfSamples.push_back("GENIEv2LFG"); Colors.push_back(GiBUUColor); Labels.push_back("Gv2 LFG    ");	
			NameOfSamples.push_back("GENIEv2EffSF"); Colors.push_back(NEUTColor); Labels.push_back("Gv2 EffSF  ");					

		}  

		//----------------------------------------//

		if (NoFSIPlusGiBUU) {

			NameOfSamples.push_back("GiBUUNoFSI"); Colors.push_back(GiBUUColor); Labels.push_back("GiB No FSI ");
			NameOfSamples.push_back("Genie_v3_0_6_NoFSI"); Colors.push_back(NuWroColor); Labels.push_back("G18 No FSI ");			
			NameOfSamples.push_back("GiBUU"); Colors.push_back(NEUTColor); Labels.push_back("GiB              ");		

		}

		//----------------------------------------//

		if (G21FSI) {

			NameOfSamples.push_back("SuSav2"); Colors.push_back(GiBUUColor); Labels.push_back("G21hN         ");
			NameOfSamples.push_back("G21hA"); Colors.push_back(NuWroColor); Labels.push_back("G21hA         ");
			NameOfSamples.push_back("G21G4"); Colors.push_back(NEUTColor); Labels.push_back("G21G4         ");
			//NameOfSamples.push_back("G21NoFSI"); Colors.push_back(NEUTColor); Labels.push_back("G21 No FSI   ");

		}		

		//----------------------------------------//

		if (All) {

			NameOfSamples.push_back("GENIEv2");	Colors.push_back(GENIEv2Color); Labels.push_back("Gv2        ");
			NameOfSamples.push_back("GENIEv2LFG"); Colors.push_back(GiBUUColor); Labels.push_back("Gv2 LFG    ");	
			NameOfSamples.push_back("GENIEv2EffSF"); Colors.push_back(NEUTColor); Labels.push_back("Gv2 EffSF  ");			
			NameOfSamples.push_back("Genie_v3_0_6_Out_Of_The_Box");	Colors.push_back(Geniev3OutOfTheBoxColor); Labels.push_back("G18 No Tune");					
			NameOfSamples.push_back("SuSav2"); Colors.push_back(SuSav2Color); Labels.push_back("G21hN      ");
			NameOfSamples.push_back("G21hA"); Colors.push_back(SuSav2Color); Labels.push_back("G21hA      ");
			NameOfSamples.push_back("G21G4"); Colors.push_back(SuSav2Color); Labels.push_back("G21G4      ");
			NameOfSamples.push_back("G21NoFSI"); Colors.push_back(SuSav2Color); Labels.push_back("G21 No FSI ");									
			NameOfSamples.push_back("Overlay9NuWro"); Colors.push_back(NuWroColor); Labels.push_back("NuWro      ");			
			NameOfSamples.push_back("GiBUU"); Colors.push_back(GiBUUColor); Labels.push_back("GiB        ");
			NameOfSamples.push_back("GiBUUNoFSI"); Colors.push_back(GiBUUColor); Labels.push_back("GiB No FSI ");
			NameOfSamples.push_back("GiBUUTscaling"); Colors.push_back(GiBUUColor); Labels.push_back("GiB Tscal  ");						
			NameOfSamples.push_back("NEUT"); Colors.push_back(NEUTColor); Labels.push_back("NEUT       ");
			NameOfSamples.push_back("NEUTv5401_RFG"); Colors.push_back(NEUTColor); Labels.push_back("NEUTv540RFG");			
			NameOfSamples.push_back("Genie_v3_0_6_Nominal"); Colors.push_back(NEUTColor); Labels.push_back("G18 Nom    ");
			NameOfSamples.push_back("Genie_v3_0_6_NoFSI"); Colors.push_back(GiBUUColor); Labels.push_back("G18 No FSI ");			
			NameOfSamples.push_back("Genie_v3_0_6_hN2018"); Colors.push_back(GENIEv2Color); Labels.push_back("G18 hN     ");
			NameOfSamples.push_back("Genie_v3_0_6_NoRPA"); Colors.push_back(NuWroColor); Labels.push_back("G18 No RPA ");
			//NameOfSamples.push_back("Genie_v3_0_6_NoCoulomb"); Colors.push_back(GENIEv3_0_4_Color); Labels.push_back("G18 No Coul");			
			NameOfSamples.push_back("Genie_v3_0_6_RFG"); Colors.push_back(GiBUUColor); Labels.push_back("G18 RFG    ");
			NameOfSamples.push_back("Genie_v3_0_6_EffSF"); Colors.push_back(GiBUUColor); Labels.push_back("G18 EffSF  ");						

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
		TString GenXSecName = PathToFiles+UBCodeVersion+"/GenXSec/All_XSecs_3D_" + PlotNames[0] + "_" + TString(std::to_string(FirstDiscrIndex)) + "_" + Runs[WhichRun] + "_"+UBCodeVersion+".root";

		// If All == true	

		if (All) {

			fGenXSec = TFile::Open(GenXSecName,"recreate");

		}			

		//----------------------------------------//

		// Open the files and grap the relevant plots

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			vector<TH1D*> CurrentPlotsFullUncReco; CurrentPlotsFullUncReco.clear();
			vector<TH1D*> CurrentPlotsXSecReco; CurrentPlotsXSecReco.clear();
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

					TH1D* histXSecReco = (TH1D*)(FileSample[WhichSample]->Get("XSecReco"+PlotNames[WhichPlot]));
					CurrentPlotsXSecReco.push_back(histXSecReco);

					TH1D* histFullUnco = (TH1D*)(FileSample[WhichSample]->Get("RecoFullUnc"+PlotNames[WhichPlot]));
					CurrentPlotsFullUncReco.push_back(histFullUncReco);										

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
					CurrentPlotsXSecReco.push_back(nullptr);
					CurrentPlotsFullUncReco.push_back(nullptr);															

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
					NameOfSamples[WhichSample] == "G21hA" ||
					NameOfSamples[WhichSample] == "G21G4" ||
					NameOfSamples[WhichSample] == "G21NoFSI" ||										
					NameOfSamples[WhichSample] == "GENIEv2" ||
					NameOfSamples[WhichSample] == "GENIEv2LFG" ||
					NameOfSamples[WhichSample] == "GENIEv2EffSF" ||					
					NameOfSamples[WhichSample] == "GENIEv3_0_4"
				) {
					FileSample.push_back(TFile::Open("../myGenieAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root")); 
				}

				if (NameOfSamples[WhichSample] == "NuWro") 
					{ FileSample.push_back(TFile::Open("../myNuWroAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root")); }

				if (NameOfSamples[WhichSample] == "GiBUU") 
					{ FileSample.push_back(TFile::Open("../myGiBUUAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root")); }

				if (NameOfSamples[WhichSample] == "GiBUUNoFSI" || NameOfSamples[WhichSample] == "GiBUUTscaling") 
					{ FileSample.push_back(TFile::Open("../myROOTGiBUU/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root")); }

				if (NameOfSamples[WhichSample] == "NEUT"  || NameOfSamples[WhichSample] == "NEUTv5401_RFG") 
					{ FileSample.push_back(TFile::Open("../myNEUTAnalysis/OutputFiles/STVAnalysis_"+NameOfSamples[WhichSample]+".root")); }

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TH1D* histTotalReco = nullptr;
					CurrentPlotsTotalReco.push_back(histTotalReco);

					TH1D* histXSecReco = nullptr;
					CurrentPlotsXSecReco.push_back(histXSecReco);

					TH1D* histFullUncReco = nullptr;
					CurrentPlotsFullUncReco.push_back(histFullUncReco);										

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

			PlotsFullUncReco.push_back(CurrentPlotsFullUncReco);
			PlotsXSecReco.push_back(CurrentPlotsXSecReco);
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

		vector< vector<TH1D*> > BeamOnXSec;
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

			int StartIndex = 0;	
			// Determine the global start bin number
			int GlobalIndex = 0;				

			if (PlotNames[WhichPlot] == "SerialECal_DeltaPTDeltaAlphaTPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsDeltaAlphaT);
				SliceBinning.push_back(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[FirstDiscrIndex]);

				for (int ig = 0; ig < FirstDiscrIndex; ig++) {

					for (int icolumn = 0; icolumn < TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[ig].size(); icolumn++) {

						GlobalIndex += TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[ig][icolumn].size()-1;

					}

				}				

			}	

			if (PlotNames[WhichPlot] == "SerialECal_DeltaPtxDeltaPtyPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsDeltaPtx);
				SliceBinning.push_back(TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices[FirstDiscrIndex]);

				for (int ig = 0; ig < FirstDiscrIndex; ig++) {

					for (int icolumn = 0; icolumn < TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices[ig].size(); icolumn++) {

						GlobalIndex += TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices[ig][icolumn].size()-1;

					}

				}				

			}			

			if (PlotNames[WhichPlot] == "SerialECal_MuonCosThetaMuonMomentumPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsMuonMomentum);
				SliceBinning.push_back(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[FirstDiscrIndex]);

				for (int ig = 0; ig < FirstDiscrIndex; ig++) {

					for (int icolumn = 0; icolumn < TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[ig].size(); icolumn++) {

						GlobalIndex += TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[ig][icolumn].size()-1;

					}

				}

			}

			if (PlotNames[WhichPlot] == "SerialECal_ProtonCosThetaProtonMomentumPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsProtonMomentum);
				SliceBinning.push_back(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[FirstDiscrIndex]);

				for (int ig = 0; ig < FirstDiscrIndex; ig++) {

					for (int icolumn = 0; icolumn < TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[ig].size(); icolumn++) {

						GlobalIndex += TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[ig][icolumn].size()-1;

					}

				}

			}																																			

			//----------------------------------------//

			BeamOnXSec.resize(N1DPlots);
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
			Cov->Scale(1./TMath::Power(MultiDimScaleFactor[PlotNames[WhichPlot]],2.)); // includes scaling factor for multi dimensional analysis

			TH2D* NormCov = (TH2D*)FileSample[0]->Get("NormUnfCov"+PlotNames[WhichPlot]);	
			NormCov->Scale(1./TMath::Power(MultiDimScaleFactor[PlotNames[WhichPlot]],2.)); // includes scaling factor for multi dimensional analysis

			TH2D* ShapeCov = (TH2D*)FileSample[0]->Get("ShapeUnfCov"+PlotNames[WhichPlot]);	
			ShapeCov->Scale(1./TMath::Power(MultiDimScaleFactor[PlotNames[WhichPlot]],2.)); // includes scaling factor for multi dimensional analysis

			//----------------------------------------//

			TH1D* UncHist = (TH1D*)(fUnc->Get("UnfUnc_"+PlotNames[WhichPlot]));

			//----------------------------------------//

			// The covariance matrix needs to be scaled by the 2D bin width

			TH2D* CovClone = (TH2D*)(Cov->Clone()); 
			TH2D* NormCovClone = (TH2D*)(NormCov->Clone());	
			TH2D* ShapeCovClone = (TH2D*)(ShapeCov->Clone());			

			int n = Cov->GetXaxis()->GetNbins();

			for (int ix = 1; ix <= n; ix++) {

				for (int iy = 1; iy <= n; iy++) {

					double WidthX = Cov->GetXaxis()->GetBinWidth(ix);
					double WidthY = Cov->GetYaxis()->GetBinWidth(iy);

					double TwoDWidth = WidthX * WidthY;

					double BinContent = Cov->GetBinContent(ix,iy);
					double NewBinContent = BinContent/TwoDWidth;

					double NormBinContent = NormCov->GetBinContent(ix,iy);
					double NormNewBinContent = NormBinContent/TwoDWidth;

					double ShapeBinContent = ShapeCov->GetBinContent(ix,iy);
					double ShapeNewBinContent = ShapeBinContent/TwoDWidth;					

					// Only for the diagonal elements
					// Add the unfolding uncertainty
					// On top of everything else
					// That is done both for the final xsec result and for the unfolded covariance
					if (ix == iy) { 
						
						// unfolded covariance matrix
						double UnfUncBin = UncHist->GetBinContent(ix);
//						double UnfUncBin = 0.;
						NewBinContent = NewBinContent + TMath::Power(UnfUncBin,2.); 
						ShapeNewBinContent = ShapeNewBinContent + TMath::Power(UnfUncBin,2.) ;						

						// xsec uncertainty
						double CurrentUnc = PlotsReco[0][WhichPlot]->GetBinError(ix);
						double NewError = TMath::Sqrt( TMath::Power(CurrentUnc,2.) + TMath::Power(UnfUncBin,2.) ) ;
						PlotsReco[0][WhichPlot]->SetBinError(ix,NewError);
						
					}

					CovClone->SetBinContent(ix,iy,NewBinContent);
					ShapeCovClone->SetBinContent(ix,iy,ShapeNewBinContent);
					NormCovClone->SetBinContent(ix,iy,NormNewBinContent);					

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
						if (iBin == 0) { SerialVectorLowBin.push_back(GlobalIndex + BinCounter); }
						// Last bin number for a given slice
						if (iBin == SliceDiscrimValue-2) { SerialVectorHighBin.push_back(GlobalIndex + BinCounter); }	

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

			BeamOnXSec[WhichPlot].resize(NSlices);
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

				NameCopy = NameCopy + "_" + TString(std::to_string(SecondDiscrIndex + NDimSlice));	

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

				TString CanvasName = Extra + "_" + PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_Slice_"+TString(std::to_string(SecondDiscrIndex + NDimSlice));
				TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
				PlotCanvas->cd();
				PlotCanvas->SetBottomMargin(0.14);
				PlotCanvas->SetTopMargin(0.12);
				PlotCanvas->SetLeftMargin(0.19);
				PlotCanvas->SetRightMargin(0.02);				

				TLegend* leg = new TLegend(0.62,0.52,0.72,0.85);
				//TLegend* legChi2 = new TLegend(0.8,0.72,0.9,0.85);			
				if (
					PlotNames[WhichPlot] == "SerialDeltaPT_DeltaAlphaTPlot" ||
					PlotNames[WhichPlot] == "SerialDeltaAlphaT_DeltaPTPlot" ||																				
					PlotNames[WhichPlot] == "SerialDeltaPtx_DeltaPtyPlot" ||
					PlotNames[WhichPlot] == "SerialDeltaPty_DeltaPtxPlot" ||
					PlotNames[WhichPlot] == "SerialProtonCosTheta_MuonCosThetaPlot" ||
					PlotNames[WhichPlot] == "SerialDeltaAlphaT_MuonCosThetaPlot" ||
					(PlotNames[WhichPlot] == "SerialMuonMomentum_MuonCosThetaPlot" && NDimSlice == 3) ||
					(PlotNames[WhichPlot] == "SerialECal_ProtonCosThetaProtonMomentumPlot" && NDimSlice == 2) ||
					(PlotNames[WhichPlot] == "SerialECal_MuonCosThetaMuonMomentumPlot" && NDimSlice == 2) ||															
					PlotNames[WhichPlot] == "SerialDeltaAlphaT_ProtonCosThetaPlot"
					) { 
						
						leg = new TLegend(0.22,0.52,0.32,0.85); 				
						
				}

				if (PlotNominal) { leg = new TLegend(0.6,0.68,0.71,0.85); }
				if (PlotNominal && PlotNames[WhichPlot] == "DeltaPtyPlot") { leg = new TLegend(0.2,0.58,0.5,0.85); }			

				leg->SetBorderSize(0);
				leg->SetTextSize(0.05);
				leg->SetTextFont(FontStyle);
				leg->SetNColumns(1);
				leg->SetMargin(0.15);

				//------------------------------------//

				// Corresponding covariance matrix

				TH2D* SliceCovMatrix = tools.Get2DHistoBins(CovClone,SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning);				

				//------------------------------------//

				// BeamOn Total Uncertainty																				

				BeamOnStatShape[WhichPlot][NDimSlice] = tools.GetHistoBins(PlotsReco[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"StatShape");										

				PrettyPlot(BeamOnStatShape[WhichPlot][NDimSlice]);														
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
					TString Chi2NdofAlt = " (" + to_string_with_precision(Chi2[WhichSample],1) + "/" + TString(std::to_string(Ndof[WhichSample])) +")";
					TLegendEntry* lGenie = leg->AddEntry(MC[WhichPlot][NDimSlice][WhichSample],Labels[WhichSample] + Chi2NdofAlt,"l");
					lGenie->SetTextColor(Colors[WhichSample]); 										
					//TLegendEntry* lGenieChi2 = legChi2->AddEntry(MC[WhichPlot][NDimSlice][WhichSample],Chi2NdofAlt,"");
					//lGenieChi2->SetTextColor(Colors[WhichSample]);					

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
				TString Chi2NdofAlt = " (" + to_string_with_precision(Chi2[0],1) + "/" + TString(std::to_string(Ndof[0])) +")";
				TLegendEntry* lGenie = leg->AddEntry(MC[WhichPlot][NDimSlice][0],Labels[0] + Chi2NdofAlt,"l");
				lGenie->SetTextColor(Colors[0]); 										
				//TLegendEntry* lGenieChi2 = legChi2->AddEntry(MC[WhichPlot][NDimSlice][0],Chi2NdofAlt,"");
				//lGenieChi2->SetTextColor(Colors[0]);					

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

				//------------------------------//

				// XSec Only

				BeamOnXSec[WhichPlot][NDimSlice] = tools.GetHistoBins(PlotsXSecReco[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"XSec");
				PrettyPlot(BeamOnXSec[WhichPlot][NDimSlice]); // includes scaling factor for multi dimensional analysis		
				BeamOnXSec[WhichPlot][NDimSlice]->SetLineColor(BeamOnColor);
				BeamOnXSec[WhichPlot][NDimSlice]->SetMarkerColor(BeamOnColor);
				BeamOnXSec[WhichPlot][NDimSlice]->SetLineWidth(1);	
				BeamOnXSec[WhichPlot][NDimSlice]->SetMarkerSize(1.);
				BeamOnXSec[WhichPlot][NDimSlice]->SetMarkerStyle(20);
				BeamOnXSec[WhichPlot][NDimSlice]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);	
				BeamOnXSec[WhichPlot][NDimSlice]->GetYaxis()->SetRangeUser(XSecRange[ MapUncorCor[ NameCopy ] ].first,XSecRange[ MapUncorCor[ NameCopy ] ].second);																			

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

				//legChi2->Draw();
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
					BeamOnXSec[WhichPlot][NDimSlice]->Write("XSecOnly_" + NameCopy);					

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

				PlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/BeamOn9/"+Extra+"MultiDimWienerSVD_Generator_TotalUnc_Data_2DXSections_"+PlotNames[WhichPlot]+"_Slice_"+TString(std::to_string(SecondDiscrIndex + NDimSlice))+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");

				//delete PlotCanvas;

				//------------------------------------//

				// Update the starting index to move to the next slice

				StartIndex += (SliceNBins+1);
				BinStartIndex += SliceNBins;
//SliceCovMatrix->Draw("coltz");
				//------------------------------------//

			} // End of the loop over the discriminators slices

			// ----------------------------------------------------------------------------------------------

		} // End of the loop over the plots

		//----------------------------------------//

		if (All) {

			fGenXSec->Close();
			cout << endl << GenXSecName << " file created" << endl << endl;	

		}			

	} // End of the loop over the runs	

} // End of the program 
