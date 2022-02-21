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
//#include <math>

#include  "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/mySimFunctions.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/MakeMyPlotPretty.cpp"

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


	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetNdivisions(8);
	h->GetYaxis()->SetTitleOffset(1.35);
	h->GetYaxis()->SetTitleSize(0.06);
	h->GetYaxis()->SetLabelSize(0.06);

	TString PlotName = h->GetName();

	PlotName.ReplaceAll("unf_","");
	PlotName.ReplaceAll("TrueUnf_","");
	PlotName.ReplaceAll("True_","");
	PlotName.ReplaceAll("True","");
	PlotName.ReplaceAll("NoSmearAlt","");			

	PlotName.ReplaceAll("_Run1","");
	PlotName.ReplaceAll("_Run2","");
	PlotName.ReplaceAll("_Run3","");
	PlotName.ReplaceAll("_Run4","");
	PlotName.ReplaceAll("_Run4a","");	
	PlotName.ReplaceAll("_Run5","");
	PlotName.ReplaceAll("_Combined","");		

	h->Scale(1./MultiDimScaleFactor[PlotName]);

}

//----------------------------------------//

void WienerSVD_OverlayGenerators(bool PlotGENIE = true, bool PlotGen = false, 
								 bool PlotGENIEFSITweaks = false, bool PlotGENIEFlagTweaks = false, 
								 bool PlotGENIECT = false, bool PlotNuclModels = false, 
								 bool PlotNuWro = false, bool PlotNominal = false, 
								 bool GiBUUComp = false, bool All = false
								) {

	//----------------------------------------//

	int DecimalAccuracy = 2;

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(4);		

	TString PathToFiles = "myXSec/";

	//----------------------------------------//

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

//	vector<TString> PlotNames;
//	PlotNames.push_back("DeltaPTPlot");
//	PlotNames.push_back("DeltaAlphaTPlot");
//	PlotNames.push_back("DeltaPhiTPlot");		 
	if (All) { PlotNames = OneDimXSec; }

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
		TString GenXSecName = PathToFiles+UBCodeVersion+"/GenXSec/All_XSecs_1D_" + Runs[WhichRun] + "_"+UBCodeVersion+".root";

		// If All == true	

		if (All) {

			fGenXSec = TFile::Open(GenXSecName,"recreate");

		}	

		//----------------------------------------//

		// Open the files and grap the relevant plots

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

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

					TH1D* histXSecReco = nullptr;
					CurrentPlotsXSecReco.push_back(histXSecReco);					

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

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

			//----------------------------------------//

			TH2D* Ac = (TH2D*)FileSample[0]->Get("Ac"+PlotNames[WhichPlot]);
			TH2D* Cov = (TH2D*)FileSample[0]->Get("UnfCov"+PlotNames[WhichPlot]);	
			Cov->Scale(1./TMath::Power(MultiDimScaleFactor[PlotNames[WhichPlot]],2.)); // includes scaling factor for multi dimensional analysis		

			//----------------------------------------//

			TH1D* UncHist = (TH1D*)(fUnc->Get("UnfUnc_"+PlotNames[WhichPlot]));
			UncHist->Scale(1./MultiDimScaleFactor[PlotNames[WhichPlot]]); // includes scaling factor for multi dimensional analysis

			//----------------------------------------//

			// The covariance matrix needs to be scaled by the 2D bin width

			TH2D* CovClone = (TH2D*)(Cov->Clone()); 

			int n = Cov->GetXaxis()->GetNbins();

			for (int ix = 1; ix <= n; ix++) {

				for (int iy = 1; iy <= n; iy++) {

					double WidthX = Cov->GetXaxis()->GetBinWidth(ix);
					double WidthY = Cov->GetYaxis()->GetBinWidth(iy);

					double TwoDWidth = WidthX * WidthY;
//					double TwoDWidth = 1.;

					double BinContent = Cov->GetBinContent(ix,iy);
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

			// -----------------------------------------------------------------------------------------------------------------------------			

			TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+"_"+Runs[WhichRun],PlotNames[WhichPlot]+"_"+Runs[WhichRun],205,34,1024,768);
			PlotCanvas->cd();

//			TPad *midPad = new TPad("midPad", "",0.005, 0.2, 0.995, 0.995);
			TPad *midPad = new TPad("midPad", "",0.005, 0., 0.995, 0.995);
			midPad->SetBottomMargin(0.14);
			midPad->SetTopMargin(0.12);
			midPad->SetLeftMargin(0.19);
			midPad->Draw();

			TLegend* leg = new TLegend(0.6,0.58,0.7,0.85);
			TLegend* legChi2 = new TLegend(0.78,0.72,0.88,0.85);			
			if (PlotNames[WhichPlot] == "MuonPhiPlot" || PlotNames[WhichPlot] == "ProtonPhiPlot" || PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot") { leg = new TLegend(0.25,0.2,0.8,0.45); }
			if (
				PlotNames[WhichPlot] == "DeltaAlphaTPlot" || PlotNames[WhichPlot] == "MuonCosThetaPlot" || 
				PlotNames[WhichPlot] == "ProtonCosThetaPlot" || PlotNames[WhichPlot] == "DeltaPtyPlot" ||
				PlotNames[WhichPlot] == "DeltaAlphaT_DeltaPT_0_20To0_40Plot" || PlotNames[WhichPlot] == "DeltaAlphaT_DeltaPT_0_40To1_00Plot" ||
				PlotNames[WhichPlot] == "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_70To1_00Plot" ||
				PlotNames[WhichPlot] == "DeltaPn_DeltaPT_0_20To0_40Plot" ||
				PlotNames[WhichPlot] == "DeltaAlphaT_MuonCosTheta_Minus1_00To0_00Plot" || 
				PlotNames[WhichPlot] == "DeltaAlphaT_MuonCosTheta_0_00To0_50Plot" ||
				PlotNames[WhichPlot] == "DeltaAlphaT_MuonCosTheta_0_50To0_75Plot" ||
				PlotNames[WhichPlot] == "MuonMomentum_MuonCosTheta_0_75To1_00Plot" ||
				PlotNames[WhichPlot] == "DeltaPty_DeltaPtx_Minus0_55ToMinus0_15Plot" ||
				PlotNames[WhichPlot] == "DeltaPty_DeltaPtx_Minus0_15To0_15Plot" || 
				PlotNames[WhichPlot] == "DeltaPty_DeltaPtx_0_15To0_55Plot" ||
				PlotNames[WhichPlot] == "DeltaPtx_DeltaPty_Minus0_75ToMinus0_15Plot" ||
				PlotNames[WhichPlot] == "DeltaPtx_DeltaPty_Minus0_15To0_15Plot" ||		
				PlotNames[WhichPlot] == "DeltaPtx_DeltaPty_0_15To0_45Plot" ||	
				PlotNames[WhichPlot] == "SerialECal_MuonCosThetaMuonMomentumPlot" ||
				PlotNames[WhichPlot] == "SerialECal_ProtonCosThetaProtonMomentumPlot" ||																	
				PlotNames[WhichPlot] == "SerialDeltaPT_MuonCosThetaPlot" ||	
				PlotNames[WhichPlot] == "SerialMuonMomentum_MuonCosThetaPlot" ||
				PlotNames[WhichPlot] == "SerialProtonMomentum_ProtonCosThetaPlot" ||
				PlotNames[WhichPlot] == "SerialDeltaPtx_DeltaPtyPlot" ||
				PlotNames[WhichPlot] == "SerialDeltaPty_DeltaPtxPlot" ||
				PlotNames[WhichPlot] == "SerialDeltaAlphaT_MuonCosThetaPlot" ||
				PlotNames[WhichPlot] == "SerialDeltaAlphaT_ProtonCosThetaPlot" ||	
				PlotNames[WhichPlot] == "DeltaAlphaT_ProtonCosTheta_Minus1_00To0_00Plot" ||
				PlotNames[WhichPlot] == "DeltaAlphaT_ProtonCosTheta_0_00To0_50Plot" ||
				PlotNames[WhichPlot] == "DeltaAlphaT_ProtonCosTheta_0_50To0_75Plot" ||				
				PlotNames[WhichPlot] == "DeltaAlphaT_ProtonCosTheta_0_75To1_00Plot" ||		
				PlotNames[WhichPlot] == "ProtonCosTheta_MuonCosTheta_Minus1_00To0_00Plot" ||
				PlotNames[WhichPlot] == "ProtonCosTheta_MuonCosTheta_0_00To0_50Plot" ||
				PlotNames[WhichPlot] == "ProtonCosTheta_MuonCosTheta_0_50To0_75Plot" ||
				PlotNames[WhichPlot] == "ProtonCosTheta_MuonCosTheta_0_75To1_00Plot" ||	
				PlotNames[WhichPlot] == "DeltaPn_DeltaPT_0_40To1_00Plot" ||																																															
				PlotNames[WhichPlot] == "SerialDeltaPT_ProtonCosThetaPlot"	
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
			if (PlotNames[WhichPlot] == "MuonPhiPlot" || PlotNames[WhichPlot] == "ProtonPhiPlot") { leg->SetNColumns(2); }
			leg->SetMargin(0.15);

			legChi2->SetBorderSize(0);
			legChi2->SetTextSize(0.03);
			legChi2->SetTextFont(FontStyle);
			legChi2->SetNColumns(1);
			if (PlotNames[WhichPlot] == "MuonPhiPlot" || PlotNames[WhichPlot] == "ProtonPhiPlot") { legChi2->SetNColumns(2); }
			legChi2->SetMargin(0.15);			

			// ------------------------------------------------------------------------------------------------------------------

			// BeamOn Total Uncertainty

			PrettyPlot(PlotsReco[0][WhichPlot]); // includes scaling factor for multi dimensional analysis

			double MaxValue = PlotsReco[0][WhichPlot]->GetMaximum();
			int MaxValueBin = LocateBinWithValue(PlotsReco[0][WhichPlot],MaxValue);
			double MaxValueError = PlotsReco[0][WhichPlot]->GetBinError(MaxValueBin);

			double MinValue = PlotsReco[0][WhichPlot]->GetMinimum();
													
			PlotsReco[0][WhichPlot]->GetYaxis()->SetRangeUser(XSecRange[PlotNames[WhichPlot]].first,XSecRange[PlotNames[WhichPlot]].second);

			PlotsReco[0][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsReco[0][WhichPlot]->SetMarkerColor(BeamOnColor);
			PlotsReco[0][WhichPlot]->SetMarkerSize(1.);
			PlotsReco[0][WhichPlot]->SetMarkerStyle(20);
			PlotsReco[0][WhichPlot]->SetLineWidth(1);	
			PlotsReco[0][WhichPlot]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);					

			midPad->cd();

			if (PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot") { 

				PlotsReco[0][WhichPlot]->GetXaxis()->SetTitle("");
				PlotsReco[0][WhichPlot]->GetXaxis()->SetLabelSize(0.);
				PlotsReco[0][WhichPlot]->GetXaxis()->SetTickSize(0.);								

				PlotsReco[0][WhichPlot]->GetYaxis()->SetTitle("#sigma #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

			}
			
			//------------------------------//

			// The N-dimensional analysis has been developed based on the bin number, not the actual range

			if (string(PlotNames[WhichPlot]).find("Serial") != std::string::npos) {	

				TString XaxisTitle = PlotsReco[0][WhichPlot]->GetXaxis()->GetTitle();
				XaxisTitle.ReplaceAll("deg","bin #");
				XaxisTitle.ReplaceAll("GeV/c","bin #");
				XaxisTitle.ReplaceAll("GeV","bin #");				
				PlotsReco[0][WhichPlot]->GetXaxis()->SetTitle(XaxisTitle);

				TString YaxisTitle = VarLabel[PlotNames[WhichPlot]];
				YaxisTitle.ReplaceAll("deg","");
				YaxisTitle.ReplaceAll("GeV/c","");
				YaxisTitle.ReplaceAll("GeV","");
				YaxisTitle.ReplaceAll("/c","");
				YaxisTitle.ReplaceAll("^{2}","");												
				PlotsReco[0][WhichPlot]->GetYaxis()->SetTitle(YaxisTitle);				

			}			

			//------------------------------//

			PlotsReco[0][WhichPlot]->Draw("e1x0 same"); // Total Unc (Shape + Stat)

			PrettyPlot(PlotsTotalReco[0][WhichPlot]); // includes scaling factor for multi dimensional analysis
			PlotsTotalReco[0][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsTotalReco[0][WhichPlot]->SetMarkerColor(BeamOnColor);
			PlotsTotalReco[0][WhichPlot]->SetLineWidth(1);			
			PlotsTotalReco[0][WhichPlot]->Draw("e1x0 same"); // Stat Only

			PrettyPlot(PlotsXSecReco[0][WhichPlot]); // includes scaling factor for multi dimensional analysis
			PlotsXSecReco[0][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsXSecReco[0][WhichPlot]->SetMarkerColor(BeamOnColor);
			PlotsXSecReco[0][WhichPlot]->SetLineWidth(1);
			PlotsXSecReco[0][WhichPlot]->SetMarkerSize(1.);
			PlotsXSecReco[0][WhichPlot]->SetMarkerStyle(20);
			PlotsXSecReco[0][WhichPlot]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);									
			//PlotsXSecReco[0][WhichPlot]->Draw("e1x0 same"); // XSec Only			
			
			PrettyPlot(PlotsNormOnly[0][WhichPlot]); // includes scaling factor for multi dimensional analysis			
			PlotsNormOnly[0][WhichPlot]->SetFillColorAlpha(kGray+1, 0.45);	
			PlotsNormOnly[0][WhichPlot]->SetLineColor(kGray+1);
			PlotsNormOnly[0][WhichPlot]->SetMarkerColor(kGray+1);
			if (PlotNames[WhichPlot] != "MuonCosThetaSingleBinPlot") { PlotsNormOnly[0][WhichPlot]->Draw("e2 hist same"); } // Norm unc Only					

			// -----------------------------------------------------------------------------------------------------------------

			// Overlay GENIE v3 + Tune

			PrettyPlot(PlotsTrue[0][WhichPlot]); // includes scaling factor for multi dimensional analysis
			PlotsTrue[0][WhichPlot]->SetLineColor(Colors[0]);
			PlotsTrue[0][WhichPlot]->SetMarkerColor(Colors[0]);

			// -----------------------------------------------------------------------------------------------------------------

			// arrays for NSamples

			double Chi2[NSamples];			
			int Ndof[NSamples];
			double pval[NSamples];

			// Clones for the NSamples-1 model predictions
			// index 0 corresponds to nominal overlay / CV

			TH1D* Clone[NSamples-1];			

			for (int WhichSample = 1; WhichSample < NSamples; WhichSample++) {

				// Remove the bin width multiplication
				UnReweight(PlotsTrue[WhichSample][WhichPlot]);
				UnReweight(QEPlotsTrue[WhichSample][WhichPlot]);	
				UnReweight(MECPlotsTrue[WhichSample][WhichPlot]);
				UnReweight(RESPlotsTrue[WhichSample][WhichPlot]);	
				UnReweight(DISPlotsTrue[WhichSample][WhichPlot]);	
				UnReweight(COHPlotsTrue[WhichSample][WhichPlot]);																	
				// Apply the additional smearing matrix Ac
				Clone[WhichSample-1] = Multiply(PlotsTrue[WhichSample][WhichPlot],Ac);
				QEPlotsTrue[WhichSample][WhichPlot] = Multiply(QEPlotsTrue[WhichSample][WhichPlot],Ac);		
				MECPlotsTrue[WhichSample][WhichPlot] = Multiply(MECPlotsTrue[WhichSample][WhichPlot],Ac);	
				RESPlotsTrue[WhichSample][WhichPlot] = Multiply(RESPlotsTrue[WhichSample][WhichPlot],Ac);	
				DISPlotsTrue[WhichSample][WhichPlot] = Multiply(DISPlotsTrue[WhichSample][WhichPlot],Ac);
				COHPlotsTrue[WhichSample][WhichPlot] = Multiply(COHPlotsTrue[WhichSample][WhichPlot],Ac);												
				// Divide again by the bin width
				Reweight(Clone[WhichSample-1]);
				Reweight(QEPlotsTrue[WhichSample][WhichPlot]);
				Reweight(MECPlotsTrue[WhichSample][WhichPlot]);
				Reweight(RESPlotsTrue[WhichSample][WhichPlot]);
				Reweight(DISPlotsTrue[WhichSample][WhichPlot]);
				Reweight(COHPlotsTrue[WhichSample][WhichPlot]);																		

				//Clone[WhichSample-1] = PlotsTrue[WhichSample][WhichPlot];				
				Clone[WhichSample-1]->SetLineColor(Colors[WhichSample]);
				Clone[WhichSample-1]->SetMarkerColor(Colors[WhichSample]);

				PrettyPlot(Clone[WhichSample-1]); // includes scaling factor for multi dimensional analysis
				Clone[WhichSample-1]->Draw("hist same");		

//				CalcChiSquared(PlotsTrue[WhichSample][WhichPlot],PlotsReco[0][WhichPlot],CovClone,Chi2[WhichSample],Ndof[WhichSample],pval[WhichSample]);
				CalcChiSquared(Clone[WhichSample-1],PlotsReco[0][WhichPlot],CovClone,Chi2[WhichSample],Ndof[WhichSample],pval[WhichSample]);
				TString Chi2NdofAlt = " (" + to_string_with_precision(Chi2[WhichSample],2) + "/" + TString(std::to_string(Ndof[WhichSample])) +")";
				if (PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot") { Chi2NdofAlt = ""; } 
//				TLegendEntry* lGenie = leg->AddEntry(Clone[WhichSample-1],Labels[WhichSample] + Chi2NdofAlt,"l");
				TLegendEntry* lGenie = leg->AddEntry(Clone[WhichSample-1],Labels[WhichSample],"l");
				lGenie->SetTextColor(Colors[WhichSample]); 										
				TLegendEntry* lGenieChi2 = legChi2->AddEntry(Clone[WhichSample-1],Chi2NdofAlt,"");
				lGenieChi2->SetTextColor(Colors[WhichSample]);

			}

			//----------------------------------------//

			// Legend & Run / POT

			double tor860_wcut = -99.;
			if (Runs[WhichRun] == "Run1") { tor860_wcut = Fulltor860_wcut_Run1; }
			if (Runs[WhichRun] == "Run2") { tor860_wcut = Fulltor860_wcut_Run2; }
			if (Runs[WhichRun] == "Run3") { tor860_wcut = Fulltor860_wcut_Run3; }
			if (Runs[WhichRun] == "Run4") { tor860_wcut = Fulltor860_wcut_Run4; }
			if (Runs[WhichRun] == "Run5") { tor860_wcut = Fulltor860_wcut_Run5; }
			if (Runs[WhichRun] == "Combined") { tor860_wcut = Fulltor860_wcut_Combined; }
			TString Label = ToString(tor860_wcut)+" POT";

			// ---------------------------------------------------------------------------------------------------------

			CalcChiSquared(PlotsTrue[0][WhichPlot],PlotsReco[0][WhichPlot],CovClone,Chi2[0],Ndof[0],pval[0]);
			TString Chi2NdofNom = " (" + to_string_with_precision(Chi2[0],2) + "/" + TString(std::to_string(Ndof[0])) +")";
			if (PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot") { Chi2NdofNom = ""; }			
//			TLegendEntry* lGenie_GenieOverlay = leg->AddEntry(PlotsTrue[0][WhichPlot],Labels[0]+Chi2NdofNom,"l");
			TLegendEntry* lGenie_GenieOverlay = leg->AddEntry(PlotsTrue[0][WhichPlot],Labels[0],"l");
			PlotsTrue[0][WhichPlot]->Draw("hist same"); lGenie_GenieOverlay->SetTextColor(Colors[0]); 
			TLegendEntry* lGenie_GenieOverlayChi2 = legChi2->AddEntry(PlotsTrue[0][WhichPlot],Chi2NdofNom,"");
			lGenie_GenieOverlayChi2->SetTextColor(Colors[0]);				

			// ---------------------------------------------------------------------------------------------------------
			// ---------------------------------------------------------------------------------------------------------

			PlotsReco[0][WhichPlot]->Draw("e1x0 same"); // BeamOn Stat Total

			leg->AddEntry(PlotsReco[0][WhichPlot],"MicroBooNE Data","ep");

			if (PlotNames[WhichPlot] != "MuonCosThetaSingleBinPlot") { 
				
				leg->AddEntry(PlotsTotalReco[0][WhichPlot],"(Stat #oplus Shape Unc)","");
			
			} else {

				leg->AddEntry(PlotsTotalReco[0][WhichPlot],"(Total Unc)","");

			}

			leg->AddEntry(PlotsReco[0][WhichPlot],Label,"");

			if (PlotNames[WhichPlot] != "MuonCosThetaSingleBinPlot") { 
				
				leg->AddEntry(PlotsNormOnly[0][WhichPlot],"Norm Unc","f"); 
			
			}			

			if (PlotNames[WhichPlot] != "MuonCosThetaSingleBinPlot") { legChi2->Draw(); }
			leg->Draw();			

			TLatex *textSlice = new TLatex();
			textSlice->SetTextFont(FontStyle);
			textSlice->SetTextSize(0.06);
			TString PlotNameDuplicate = PlotNames[WhichPlot];
			TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("Reco","") ;
			textSlice->DrawLatexNDC(0.24, 0.92, LatexLabel[ ReducedPlotName ] );	

			//----------------------------------------//

			// If the option All == true is activated, store all the relevant xsecs in a file

			if (All) {

				fGenXSec->cd();

				// Unfolded covariance matrix
				CovClone->Write("UnfCov_" + PlotNames[WhichPlot]);

				// Data
				PlotsReco[0][WhichPlot]->Write("StatShape_" + PlotNames[WhichPlot]); // Stat + Shape
				PlotsNormOnly[0][WhichPlot]->Write("NormOnly_" + PlotNames[WhichPlot]); // Norm only
				PlotsTotalReco[0][WhichPlot]->Write("StatOnly_" + PlotNames[WhichPlot]); // Stat only
				PlotsXSecReco[0][WhichPlot]->Write("XSecOnly_" + PlotNames[WhichPlot]); // Stat only				

				// Overlay GENIE
				PlotsTrue[0][WhichPlot]->Write("OverlayGENIE_" + PlotNames[WhichPlot]);
				QEPlotsTrue[0][WhichPlot]->Write("QEOverlayGENIE_" + PlotNames[WhichPlot]);
				MECPlotsTrue[0][WhichPlot]->Write("MECOverlayGENIE_" + PlotNames[WhichPlot]);	
				RESPlotsTrue[0][WhichPlot]->Write("RESOverlayGENIE_" + PlotNames[WhichPlot]);
				DISPlotsTrue[0][WhichPlot]->Write("DISOverlayGENIE_" + PlotNames[WhichPlot]);	
				COHPlotsTrue[0][WhichPlot]->Write("COHOverlayGENIE_" + PlotNames[WhichPlot]);																		

				// Store the remaining generator xsecs
				// 0 is the overlay that has been stored above

				for (int igen = 1; igen < NSamples; igen++ ) {

					TString GenName = NameOfSamples[igen] + "_" + PlotNames[WhichPlot];
					Clone[igen-1]->Write(GenName);
					QEPlotsTrue[igen][WhichPlot]->Write("QE" + GenName);
					MECPlotsTrue[igen][WhichPlot]->Write("MEC" + GenName);
					RESPlotsTrue[igen][WhichPlot]->Write("RES" + GenName);
					DISPlotsTrue[igen][WhichPlot]->Write("DIS" + GenName);										
					COHPlotsTrue[igen][WhichPlot]->Write("COH" + GenName);															

				}

			}		

			//----------------------------------------//

			// Saving the canvas with the data (total uncertainties) vs overlay & generator predictions

			PlotCanvas->SaveAs("./myPlots/pdf/"+UBCodeVersion+"/BeamOn9/"+Extra+"WienerSVD_Generator_TotalUnc_Data_XSections_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");

			delete PlotCanvas;

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
