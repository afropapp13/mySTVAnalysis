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

#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision

#include "ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "ubana/AnalysisCode/Secondary_Code/GlobalSettings.cpp"
#include "ubana/AnalysisCode/Secondary_Code/myFunctions.cpp"

#include "ubana/myClasses/Util.h"

//---------------------------------------------//

void MakeBeautiful(TH1D* h){


}

//---------------------------------------------//

void Correlations_XSecSignalBkgUnc() {

	//---------------------------------------------//

	GlobalSettings();
	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	TGaxis::SetMaxDigits(3);
	gStyle->SetPaintTextFormat("4.4f");	

	std::cout << std::setprecision(3) << std::fixed;	

	TString ExactFileLocation = PathToFiles+CutExtension;

	//---------------------------------------------//

	vector<TString> PlotNames;
	PlotNames.push_back("DeltaPTPlot");
	PlotNames.push_back("MuonCosThetaSingleBinPlot");	 

	const int nplots = PlotNames.size();

	//---------------------------------------------//

	// XSec variations

	std::vector<TString> Vars;

	//Vars.push_back("AxFFCCQEshape_UBGenie");
	//Vars.push_back("DecayAngMEC_UBGenie");
	//Vars.push_back("NormCCCOH_UBGenie");
	//Vars.push_back("NormNCCOH_UBGenie");
	Vars.push_back("RPA_CCQE_UBGenie");
	//Vars.push_back("ThetaDelta2NRad_UBGenie");
	//Vars.push_back("Theta_Delta2Npi_UBGenie");
	//Vars.push_back("VecFFCCQEshape_UBGenie");
	Vars.push_back("XSecShape_CCMEC_UBGenie");
	//Vars.push_back("All_UBGenie");

	int nvars = Vars.size();							

	//---------------------------------------------//	

	// Loop over the plots

	for (int iplot = 0; iplot < nplots; iplot++) {						

		//---------------------------------------------//									

		// Loop over the xsec variations

		for (int ivar = 0; ivar < nvars; ivar++) {

			//---------------------------------------------//				

			// Total covariance matrix

			TString ivar_cov_xsecfile_name = MigrationMatrixPath+"IndividualWienerSVD_" + Vars[ivar] + "_CovarianceMatrices_Overlay9_Combined_"+UBCodeVersion+".root";
			cout << ivar_cov_xsecfile_name << endl;
			TFile* ivar_nomcov_xsecfile = TFile::Open(ivar_cov_xsecfile_name,"readonly");

			TH2D* ivar_nomfraccovmatrix = (TH2D*)(ivar_nomcov_xsecfile->Get(Vars[ivar] + "_FracCovariance_"+PlotNames[iplot]+"_Combined"));
			TH2D* signal_ivar_nomfraccovmatrix = (TH2D*)(ivar_nomcov_xsecfile->Get(Vars[ivar] + "_SignalFracCovariance_"+PlotNames[iplot]+"_Combined"));
			TH2D* bkg_ivar_nomfraccovmatrix = (TH2D*)(ivar_nomcov_xsecfile->Get(Vars[ivar] + "_BkgFracCovariance_"+PlotNames[iplot]+"_Combined"));							

			TString canvas_name_total = "total_fraccov_"+PlotNames[iplot]+"_"+Vars[ivar];
			TCanvas* plot_canvas_total = new TCanvas(canvas_name_total,canvas_name_total,205,34,1024,768);
			plot_canvas_total->cd();
			plot_canvas_total->SetRightMargin(0.15);					
			ivar_nomfraccovmatrix->SetMarkerColor(kWhite);
			ivar_nomfraccovmatrix->Draw("coltz");	
			ivar_nomfraccovmatrix->SaveAs(PlotPath+"Overlay9/"+canvas_name_total+".pdf");

			TString canvas_name_signal = "signal_fraccov_"+PlotNames[iplot]+"_"+Vars[ivar];
			TCanvas* plot_canvas_signal = new TCanvas(canvas_name_signal,canvas_name_signal,205,34,1024,768);
			plot_canvas_signal->cd();
			plot_canvas_signal->SetRightMargin(0.15);					
			signal_ivar_nomfraccovmatrix->SetMarkerColor(kWhite);
			signal_ivar_nomfraccovmatrix->Draw("coltz");	
			signal_ivar_nomfraccovmatrix->SaveAs(PlotPath+"Overlay9/"+canvas_name_signal+".pdf");

			TString canvas_name_bkg = "bkg_fraccov_"+PlotNames[iplot]+"_"+Vars[ivar];
			TCanvas* plot_canvas_bkg = new TCanvas(canvas_name_bkg,canvas_name_bkg,205,34,1024,768);
			plot_canvas_bkg->cd();
			plot_canvas_bkg->SetRightMargin(0.15);					
			bkg_ivar_nomfraccovmatrix->SetMarkerColor(kWhite);
			bkg_ivar_nomfraccovmatrix->Draw("coltz");	
			bkg_ivar_nomfraccovmatrix->SaveAs(PlotPath+"Overlay9/"+canvas_name_bkg+".pdf");											

			//---------------------------------------------//

		} // end of the loop over the xsec variations

	}

	//---------------------------------------------//

} // End of the program
