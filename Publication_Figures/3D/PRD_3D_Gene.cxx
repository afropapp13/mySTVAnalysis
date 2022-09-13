#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TLine.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../../../myClasses/Constants.h"
#include "../../../myClasses/Util.h"

using namespace std;
using namespace Constants;

#include "/home/afroditi/Dropbox/PhD/Secondary_Code/mySimFunctions.cpp"

//----------------------------------------//

void PRD_3D_Gene() {

	//----------------------------------------//

	int DecimalAccuracy = 2;
	int FontStyle = 132;
	double TextSize = 0.08;	
	TString Runs = "Combined";	

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(4);		

	double XMinCanvas = 0.1;
	double XMaxCanvas = 0.99;
	double YMinCanvas = 0.1;
	double YMaxCanvas = 0.99;	

	double YaxisRangeMin = 0.1;	
	double YaxisRangeMax = 2.9;

	double XaxisRangeMin = ArrayNBinsECal[0];	
	double XaxisRangeMax = ArrayNBinsECal[NBinsECal];			

	//----------------------------------------//

//	vector<int> MarkerStyles{20,21,22,33,34};
	vector<int> MarkerStyles{20,20,20,20,20};	

	//----------------------------------------//

	// MC Samples to compare

	vector<TString> MCSample; vector<TString> Label; vector<int> MCColors;  vector<int> LineStyle;

//	MCSample.push_back("GiBUUNoFSI"); Label.push_back("GiB No FSI"); MCColors.push_back(NuWroColor); LineStyle.push_back(kDashed);
//	MCSample.push_back("Genie_v3_0_6_NoFSI"); Label.push_back("G18 No FSI"); MCColors.push_back(OverlayColor); LineStyle.push_back(kDashed);	
	MCSample.push_back("OverlayGENIE"); Label.push_back("G18    "); MCColors.push_back(OverlayColor); LineStyle.push_back(kSolid);
	MCSample.push_back("GiBUU"); Label.push_back("GiBUU"); MCColors.push_back(GiBUUColor); LineStyle.push_back(kSolid);	
//	MCSample.push_back("GiBUUTscaling"); Label.push_back("GiBUUTscaling");	
	MCSample.push_back("NEUT");  Label.push_back("NEUT"); MCColors.push_back(NEUTColor); LineStyle.push_back(kSolid);
	MCSample.push_back("Overlay9NuWro");  Label.push_back("NuWro"); MCColors.push_back(NuWroColor); LineStyle.push_back(kSolid);	
//	MCSample.push_back("NEUTv5401_RFG");  Label.push_back("NEUTv5401_RFG");	
//	MCSample.push_back("Overlay9NuWro"); Label.push_back("NuWro");
//	MCSample.push_back("GENIEv2"); Label.push_back("Gv2");
//	MCSample.push_back("GENIEv2LFG"); Label.push_back("Gv2 LFG");
//	MCSample.push_back("GENIEv2EffSF"); Label.push_back("Gv2 EffSF");		
//	MCSample.push_back("Genie_v3_0_6_Out_Of_The_Box"); Label.push_back("G18 No Tune");					
//	MCSample.push_back("SuSav2"); Label.push_back("G21hN");
//	MCSample.push_back("G21hA"); Label.push_back("G21hA");	
//	MCSample.push_back("G21G4"); Label.push_back("G21G4");
//	MCSample.push_back("G21NoFSI"); Label.push_back("G21NoFSI");		
//	MCSample.push_back("Genie_v3_0_6_hN2018"); Label.push_back("G18 hN Tune");
//	MCSample.push_back("Genie_v3_0_6_NoRPA"); Label.push_back("G18 No RPA Tune");
//	MCSample.push_back("Genie_v3_0_6_RFG"); Label.push_back("G18 RFG Tune");
//	MCSample.push_back("Genie_v3_0_6_EffSF"); Label.push_back("G18 EffSF Tune");			

	int NMC = MCSample.size();	

	//----------------------------------------//

	// File with all the xsecs

	TString XSecFileName = "../../myXSec/" + UBCodeVersion + "/GenXSec/All_XSecs_Combined_" + UBCodeVersion + ".root";
	TFile* fXSec = new TFile(XSecFileName,"readonly");	

	//----------------------------------------//

	std::vector< std::vector< std::vector<TString> > > ThreeDimCanvases;

	ThreeDimCanvases.push_back( { 
												// ECal in DeltaPT & DeltaAlphaT slices
												// DeltaPT < 0.2 GeV/c
												{"SerialECal_DeltaPTDeltaAlphaTPlot_0","SerialECal_DeltaPTDeltaAlphaTPlot_1","SerialECal_DeltaPTDeltaAlphaTPlot_2","SerialECal_DeltaPTDeltaAlphaTPlot_3"},
												// 0.2 < DeltaPT < 0.4 GeV/c
												{"SerialECal_DeltaPTDeltaAlphaTPlot_4","SerialECal_DeltaPTDeltaAlphaTPlot_5","SerialECal_DeltaPTDeltaAlphaTPlot_6","SerialECal_DeltaPTDeltaAlphaTPlot_7"},
												// DeltaPT > 0.4 GeV/c
												{"SerialECal_DeltaPTDeltaAlphaTPlot_8","SerialECal_DeltaPTDeltaAlphaTPlot_9","SerialECal_DeltaPTDeltaAlphaTPlot_10","SerialECal_DeltaPTDeltaAlphaTPlot_11"}																							
											}
	);

	ThreeDimCanvases.push_back( { 
												// ECal in DeltaPtx & DeltaPty slices
												// DeltaPty < -0.15 GeV/c
												{"SerialECal_DeltaPtxDeltaPtyPlot_0","SerialECal_DeltaPtxDeltaPtyPlot_3","SerialECal_DeltaPtxDeltaPtyPlot_6"},
												// -0.15 < DeltaPty < 0.15 GeV/c
												{"SerialECal_DeltaPtxDeltaPtyPlot_1","SerialECal_DeltaPtxDeltaPtyPlot_4","SerialECal_DeltaPtxDeltaPtyPlot_7"},
												// DeltaPty > 0.15 GeV/c
												{"SerialECal_DeltaPtxDeltaPtyPlot_2","SerialECal_DeltaPtxDeltaPtyPlot_5","SerialECal_DeltaPtxDeltaPtyPlot_8"}																							
											}
	);	

	ThreeDimCanvases.push_back( { 
												// ECal in MuonCosTheta & MuonMomentum slices
												// -1 < MuonCosTheta < 0
												{"SerialECal_MuonCosThetaMuonMomentumPlot_0","SerialECal_MuonCosThetaMuonMomentumPlot_1","SerialECal_MuonCosThetaMuonMomentumPlot_2"},
												// 0 < MuonCosTheta < 0.5
												{"SerialECal_MuonCosThetaMuonMomentumPlot_3","SerialECal_MuonCosThetaMuonMomentumPlot_4","SerialECal_MuonCosThetaMuonMomentumPlot_5"},
												// 0.5 < MuonCosTheta < 0.75
												{"SerialECal_MuonCosThetaMuonMomentumPlot_6","SerialECal_MuonCosThetaMuonMomentumPlot_7","SerialECal_MuonCosThetaMuonMomentumPlot_8"},
												// 0.75 < MuonCosTheta < 1
												{"SerialECal_MuonCosThetaMuonMomentumPlot_9","SerialECal_MuonCosThetaMuonMomentumPlot_10","SerialECal_MuonCosThetaMuonMomentumPlot_11"},																																			
											}
	);	

	ThreeDimCanvases.push_back( { 
												// ECal in ProtonCosTheta & ProtonMomentum slices
												// -1 < ProtonCosTheta < 0
												{"SerialECal_ProtonCosThetaProtonMomentumPlot_0","SerialECal_ProtonCosThetaProtonMomentumPlot_1","SerialECal_ProtonCosThetaProtonMomentumPlot_2"},
												// 0 < ProtonCosTheta < 0.5
												{"SerialECal_ProtonCosThetaProtonMomentumPlot_3","SerialECal_ProtonCosThetaProtonMomentumPlot_4","SerialECal_ProtonCosThetaProtonMomentumPlot_5"},
												// 0.5 < ProtonCosTheta < 0.75
												{"SerialECal_ProtonCosThetaProtonMomentumPlot_6","SerialECal_ProtonCosThetaProtonMomentumPlot_7","SerialECal_ProtonCosThetaProtonMomentumPlot_8"},
												// 0.75 < ProtonCosTheta < 1
												{"SerialECal_ProtonCosThetaProtonMomentumPlot_9","SerialECal_ProtonCosThetaProtonMomentumPlot_10","SerialECal_ProtonCosThetaProtonMomentumPlot_11"},																																			
											}
	);	

	const int NCanvases = ThreeDimCanvases.size();
	cout << "Number of 3D Canvases = " << NCanvases << endl;

	//----------------------------------------//

	// TLine at 1

	TLine* line = new TLine(XaxisRangeMin,1,XaxisRangeMax,1);
	line->SetLineColor(kBlack);
	line->SetLineStyle(kDashed);
	line->SetLineWidth(2);	

	//--------------------------------------------//	

	// Loop over the 3D canvases

	for (int icanvas = 0; icanvas < NCanvases; icanvas++) {

		// Draw the main canvas

		TString PdfName = ThreeDimCanvases[icanvas][0][0]; // grab the 1st entry
		PdfName.ReplaceAll("Plot_0","");
		PdfName.ReplaceAll("Serial","");		

		TString CanvasName = "3D_PRD_" + PdfName;
		TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
		PlotCanvas->cd();
		PlotCanvas->SetBottomMargin(0.);
		PlotCanvas->SetTopMargin(0.);
		PlotCanvas->SetLeftMargin(0.1);
		PlotCanvas->SetRightMargin(0.);				
		PlotCanvas->Draw();	

		// Draw the pads for each slice

		int NPads = ThreeDimCanvases[icanvas].size();
		double XWidth = (XMaxCanvas - XMinCanvas) / NPads;			

		for (int ipad = 0; ipad < NPads; ipad++) {

			// Secondary slices

			int NPlots = ThreeDimCanvases[icanvas][ipad].size();

			TPad* pad[NPads][NPlots];

			double YWidth = (YMaxCanvas - YMinCanvas) / NPlots;				

			for (int iplot = 0; iplot < NPlots; iplot++) {

				PlotCanvas->cd();
				TString PadName = "pad_" + ToString(ipad) + "_" + ToString(iplot);
				pad[ipad][iplot] = new TPad(PadName,PadName,XMinCanvas + ipad * XWidth,YMinCanvas + iplot * YWidth,XMinCanvas + (ipad + 1) * XWidth,YMinCanvas + (iplot + 1) * YWidth);
				pad[ipad][iplot]->SetBottomMargin(0.);
				pad[ipad][iplot]->SetTopMargin(0.);
				pad[ipad][iplot]->SetLeftMargin(0.);
				pad[ipad][iplot]->SetRightMargin(0.);
				pad[ipad][iplot]->SetFillColor(kWhite);							
				pad[ipad][iplot]->Draw();					

				TH1D* Data = (TH1D*)(fXSec->Get("FullUnc_" + ThreeDimCanvases[icanvas][ipad][iplot]) );

				Data->GetXaxis()->CenterTitle();
				Data->GetXaxis()->SetLabelFont(FontStyle);
				Data->GetXaxis()->SetTitleFont(FontStyle);
				Data->GetXaxis()->SetLabelSize(0.);
				Data->GetXaxis()->SetTitleSize(0.);
				Data->GetXaxis()->SetNdivisions(6);					

				Data->GetYaxis()->CenterTitle();
				Data->GetYaxis()->SetLabelFont(FontStyle);
				Data->GetYaxis()->SetTitleFont(FontStyle);
				Data->GetYaxis()->SetLabelSize(TextSize);
				Data->GetYaxis()->SetTitleSize(TextSize);	
				Data->GetYaxis()->SetNdivisions(6);		

				TH2D* Cov = (TH2D*)fXSec->Get("UnfCov_" + ThreeDimCanvases[icanvas][ipad][iplot]);															

				//--------------------------------------------//

				TH1D* MCPlot[NMC];
				TH1D* DataClone[NMC];

				//--------------------------------------------//

				// TLegend & Chi2 calculation

				TLegend* leg = new TLegend(0.05,0.05,0.85,0.25);
				leg->SetBorderSize(0);
				leg->SetFillStyle(0);
				leg->SetNColumns(2);
				leg->SetTextFont(FontStyle);
				leg->SetTextSize(0.06);													

				double Chi2[NMC];				
				int Ndof[NMC];
				double pval[NMC];				

				//--------------------------------------------//												

				for (int imc = 0; imc < NMC; imc++) {									

					DataClone[imc] = (TH1D*)(Data->Clone());
					MCPlot[imc] = (TH1D*)(fXSec->Get(MCSample[imc] + "_" + ThreeDimCanvases[icanvas][ipad][iplot]));
					DataClone[imc]->Divide(MCPlot[imc]);

					DataClone[imc]->SetLineColor(MCColors[imc]);
					DataClone[imc]->SetMarkerColor(MCColors[imc]);						
					DataClone[imc]->SetMarkerStyle(MarkerStyles[iplot]);

					pad[ipad][iplot]->cd();
					DataClone[imc]->GetYaxis()->SetRangeUser(YaxisRangeMin,YaxisRangeMax);
					DataClone[imc]->Draw("e same");

					CalcChiSquared(Data,MCPlot[imc],Cov,Chi2[imc],Ndof[imc],pval[imc]);
					TString Chi2NdofAlt = "(" + to_string_with_precision(Chi2[imc],1) + "/" + TString(std::to_string(Ndof[imc])) +")";
					TLegendEntry* lGenie = leg->AddEntry(DataClone[imc],Label[imc] + Chi2NdofAlt,"l");
					lGenie->SetTextColor(MCColors[imc]); 										

					// ---------------------------------------- //		

					// Plotting the axes			

					if (imc == NMC-1) {

						// ---------------------------------------- //

						// Plotting the x axis						

						if (iplot == NPlots-1) {

							PlotCanvas->cd();
							TGaxis *xaxis = new TGaxis(XMinCanvas + ipad * XWidth,XMinCanvas,XMinCanvas + (ipad + 1) * XWidth,XMinCanvas,XaxisRangeMin,XaxisRangeMax,4,"");
							xaxis->SetLabelFont(FontStyle);
							xaxis->Draw();								

						}	

						// ---------------------------------------- //	

						// Plotting the y axis	

						if (ipad == 0) {

							PlotCanvas->cd();
							TGaxis *yaxis = new TGaxis(XMinCanvas,YMinCanvas + iplot * YWidth,XMinCanvas,YMinCanvas + (iplot+1) * YWidth,YaxisRangeMin,YaxisRangeMax,4,"");
							yaxis->SetLabelFont(FontStyle);
							yaxis->Draw();								

						}

						// ---------------------------------------- //												

					}	

					// ---------------------------------------- //																		

				} // End of the loop over the different generators/MC

				//--------------------------------------------//

				// Draw the line at 1
				// and then draw the ratio plots again so that they can be on top

				pad[ipad][iplot]->cd();
				line->Draw("same");

				for (int imc = 0; imc < NMC; imc++) {

					DataClone[imc]->Draw("e same");

				}

				// Legend

				leg->Draw();


				TLatex SliceLabel;
				SliceLabel.SetTextFont(FontStyle);
				SliceLabel.SetTextSize(0.06);
				SliceLabel.DrawLatexNDC(0.06,0.9,LatexLabel[ MapUncorCor[ ThreeDimCanvases[icanvas][ipad][iplot] ] ]);

				//--------------------------------------------//								

			} // End of the loop over the secondary slices

		} // End of the loop over the slices

		//----------------------------------------//

		// TLatex pads for the axes labels

		TLatex Ytex;
		Ytex.SetTextFont(FontStyle);
		Ytex.SetTextSize(0.05);
		PlotCanvas->cd();
		Ytex.SetTextAngle(90);			
		Ytex.DrawLatexNDC(0.04,0.37,"Data/Simulation");

		TLatex Xtex;
		Xtex.SetTextFont(FontStyle);
		Xtex.SetTextSize(0.05);
		PlotCanvas->cd();			
		Xtex.DrawLatexNDC(0.47,0.01,"E^{Cal} [GeV/c]");		

		//----------------------------------------//
		
		PlotCanvas->SaveAs("pdf/"+PdfName+".pdf");
		//delete Canvas;					

		//----------------------------------------//					

	} // End of the loop over the 3D canvases

	//----------------------------------------//

} // End of the program 