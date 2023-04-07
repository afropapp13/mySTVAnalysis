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

//---------------------------------------------//

TH2D* Multiply(TH2D* rotmatrix, TH2D* covmatrix) {

	TH2D* rotmatrix_clone = (TH2D*)(rotmatrix->Clone());

	int XBins = covmatrix->GetXaxis()->GetNbins();
	int YBins = covmatrix->GetYaxis()->GetNbins();

	if (XBins != YBins) { std::cout << "Not symmetric matrix" << std::endl; }

	TMatrixD rot(XBins,XBins);
	TMatrixD cov(XBins,XBins);

	H2M(rotmatrix, rot, kFALSE); // X axis: Reco, Y axis: True
	H2M(covmatrix, cov, kFALSE); // X axis: Reco, Y axis: True

	TMatrixD product = rot * cov;
	M2H(product, rotmatrix_clone);	

	return rotmatrix_clone;

}

//---------------------------------------------//

void WienerSVD_XSecVars() {

	//---------------------------------------------//

	GlobalSettings();
	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	TGaxis::SetMaxDigits(3);

	TString ExactFileLocation = PathToFiles+CutExtension;

	//---------------------------------------------//

	vector<TString> PlotNames;
	PlotNames.push_back("MuonCosThetaSingleBinPlot"); 
	const int nplots = PlotNames.size();

	//---------------------------------------------//

	// XSec variations

	std::vector<TString> Vars;

	Vars.push_back("AxFFCCQEshape_UBGenie");
	Vars.push_back("DecayAngMEC_UBGenie");
	Vars.push_back("NormCCCOH_UBGenie");
	Vars.push_back("NormNCCOH_UBGenie");
	Vars.push_back("RPA_CCQE_UBGenie");
	Vars.push_back("ThetaDelta2NRad_UBGenie");
	Vars.push_back("Theta_Delta2Npi_UBGenie");
	Vars.push_back("VecFFCCQEshape_UBGenie");
	Vars.push_back("XSecShape_CCMEC_UBGenie");
	Vars.push_back("All_UBGenie");

	int nvars = Vars.size();

	//---------------------------------------------//

	// Detailed XSec variations

	std::vector<TString> DetVars;

	DetVars.push_back("AGKYpT1pi_UBGenie");
	DetVars.push_back("AGKYxF1pi_UBGenie");
	DetVars.push_back("AhtBY_UBGenie");
	DetVars.push_back("BhtBY_UBGenie");
	DetVars.push_back("CV1uBY_UBGenie");
	DetVars.push_back("CV2uBY_UBGenie");
	DetVars.push_back("EtaNCEL_UBGenie");
	DetVars.push_back("FrAbs_N_UBGenie");
	DetVars.push_back("FrAbs_pi_UBGenie");
	DetVars.push_back("FrCEx_N_UBGenie");
	DetVars.push_back("FrCEx_pi_UBGenie");
	DetVars.push_back("FrInel_N_UBGenie");
	DetVars.push_back("FrInel_pi_UBGenie");
	DetVars.push_back("FrPiProd_N_UBGenie");
	DetVars.push_back("FrPiProd_pi_UBGenie");
	DetVars.push_back("FracDelta_CCMEC_UBGenie");
	DetVars.push_back("FracPN_CCMEC_UBGenie");
	DetVars.push_back("MFP_N_UBGenie");
	DetVars.push_back("MFP_pi_UBGenie");
	DetVars.push_back("MaCCQE_UBGenie");
	DetVars.push_back("MaCCRES_UBGenie");
	DetVars.push_back("MaNCEL_UBGenie");
	DetVars.push_back("MaNCRES_UBGenie");
	DetVars.push_back("MvCCRES_UBGenie");
	DetVars.push_back("MvNCRES_UBGenie");
	DetVars.push_back("NonRESBGvbarnCC1pi_UBGenie");
	DetVars.push_back("NonRESBGvbarnCC2pi_UBGenie");
	DetVars.push_back("NonRESBGvbarnNC1pi_UBGenie");
	DetVars.push_back("NonRESBGvbarnNC2pi_UBGenie");
	DetVars.push_back("NonRESBGvbarpCC1pi_UBGenie");
	DetVars.push_back("NonRESBGvbarpCC2pi_UBGenie");
	DetVars.push_back("NonRESBGvbarpNC1pi_UBGenie");
	DetVars.push_back("NonRESBGvbarpNC2pi_UBGenie");
	DetVars.push_back("NonRESBGvnCC1pi_UBGenie");
	DetVars.push_back("NonRESBGvnCC2pi_UBGenie");
	DetVars.push_back("NonRESBGvnNC1pi_UBGenie");
	DetVars.push_back("NonRESBGvnNC2pi_UBGenie");
	DetVars.push_back("NonRESBGvpCC1pi_UBGenie");
	DetVars.push_back("NonRESBGvpCC2pi_UBGenie");
	DetVars.push_back("NonRESBGvpNC1pi_UBGenie");
	DetVars.push_back("NonRESBGvpNC2pi_UBGenie");
	DetVars.push_back("NormCCMEC_UBGenie");
	DetVars.push_back("NormNCMEC_UBGenie");
	DetVars.push_back("RDecBR1eta_UBGenie");
	DetVars.push_back("RDecBR1gamma_UBGenie");					

	int ndetvars = DetVars.size();	

	//---------------------------------------------//

	// XSec file (includes covariance rotation matrix and its transpose)

	TString xsecfile_name = PathToExtractedXSec+"WienerSVD_ExtractedXSec_Overlay9_Combined_"+UBCodeVersion+".root";
	TFile* xsecfile = TFile::Open(xsecfile_name,"readonly");	

	// Two files with covariances

	// 1st: nominal covariances with full stats, including All_UBGenie

	TString nomcov_xsecfile_name = MigrationMatrixPath+"WienerSVD_XSec_CovarianceMatrices_Overlay9_Combined_"+UBCodeVersion+".root";
	TFile* nomcov_xsecfile = TFile::Open(nomcov_xsecfile_name,"readonly");		

	//---------------------------------------------//	

	// Loop over the plots

	for (int iplot = 0; iplot < nplots; iplot++) {	

		//---------------------------------------------// 

		// Rotation matrices to rotate into regularized phase space

		TH2D* rotmatrix = (TH2D*)(xsecfile->Get("CovRot"+PlotNames[iplot]));
		TH2D* transrotmatrix = (TH2D*)(xsecfile->Get("TranspCovRot"+PlotNames[iplot]));	

		// Total xsec syst		

		TH2D* nomcovmatrix = (TH2D*)(nomcov_xsecfile->Get("XSec_Covariance_"+PlotNames[iplot]+"_Combined"));		

		//---------------------------------------------//	

		double 	All_UBGenieValue = 0.;
		double 	TotalXSecFracUnc = 0.;								

		// Loop over the xsec variations

		for (int ivar = 0; ivar < nvars; ivar++) {

			//---------------------------------------------//				

			TString ivar_cov_xsecfile_name = MigrationMatrixPath+"IndividualWienerSVD_" + Vars[ivar] + "_CovarianceMatrices_Overlay9_Combined_"+UBCodeVersion+".root";
			TFile* ivar_nomcov_xsecfile = TFile::Open(ivar_cov_xsecfile_name,"readonly");

			TH2D* ivar_nomcovmatrix = (TH2D*)(ivar_nomcov_xsecfile->Get(Vars[ivar] + "_Covariance_"+PlotNames[iplot]+"_Combined"));
			TH2D* ivar_nomfraccovmatrix = (TH2D*)(ivar_nomcov_xsecfile->Get(Vars[ivar] + "_FracCovariance_"+PlotNames[iplot]+"_Combined"));			

			TH2D* ivar_unf_nomcovmatrix = Multiply( transrotmatrix , Multiply(ivar_nomcovmatrix,rotmatrix) );
			TH2D* ivar_unf_nomfraccovmatrix = Multiply( transrotmatrix , Multiply(ivar_nomfraccovmatrix,rotmatrix) );			

			cout << Vars[ivar] << " = " << TMath::Sqrt( ivar_nomcovmatrix->GetBinContent(1,1) );
			cout << ", fractional contrinution = " << TMath::Sqrt( ivar_nomfraccovmatrix->GetBinContent(1,1) ) * 100. << endl;

			if (Vars[ivar] == "All_UBGenie") { All_UBGenieValue = ivar_unf_nomcovmatrix->GetBinContent(1,1); }

			TotalXSecFracUnc += ivar_nomfraccovmatrix->GetBinContent(1,1);

			//---------------------------------------------//

		} // end of the loop over the xsec variations

		cout << "Total summed xsec fractional uncertainty = " << TMath::Sqrt(TotalXSecFracUnc) * 100. << "%" << endl;
		cout << "Total xsec fractional uncertainty = " << TMath::Sqrt(nomcovmatrix->GetBinContent(1,1)) * 100. << "%" << endl;		

		//---------------------------------------------//									

		double 	DetAll_UBGenieValue = 0.;	
/*
		// Loop over the detailed xsec variations

		for (int idetvar = 0; idetvar < ndetvars; idetvar++) {

				// 2nd: smaller production with detailed xsec variations

				TString redcov_xsecfile_name = MigrationMatrixPath+"WienerSVD_DetailedXSec_CovarianceMatrices_Overlay9_Run1_DecompXSecUnc_"+UBCodeVersion+".root";
				TFile* redcov_xsecfile = TFile::Open(redcov_xsecfile_name,"readonly");

			//---------------------------------------------//

			TH2D* detcovmatrix = (TH2D*)(nomcov_xsecfile->Get(DetVars[idetvar] + "_Covariance_"+PlotNames[iplot]+"_Run1_DecompXSecUnc"));				

			DetAll_UBGenieValue += detcovmatrix->GetBinContent(1,1);

			//---------------------------------------------//				

		} // end of the loop over the detailed xsec variations			

		//---------------------------------------------//	

		// Scaling factor to account for differences due to statistic fluctuations between the two xsec systematics files

		double sf = DetAll_UBGenieValue / DetAll_UBGenieValue;								

		// Loop over the detailed xsec variations

		for (int idetvar = 0; idetvar < ndetvars; idetvar++) {

			//---------------------------------------------//				

			TH2D* detcovmatrix = (TH2D*)(nomcov_xsecfile->Get(DetVars[idetvar] + "_Covariance_"+PlotNames[iplot]+"_Run1_DecompXSecUnc"));

			detcovmatrix->Scale(sf);

			TH2D* unf_detcovmatrix = Multiply( transrotmatrix , Multiply(detcovmatrix,rotmatrix) );

			cout << DetVars[idetvar] << " = " << TMath::Abs( unf_detcovmatrix->GetBinContent(1,1) ) << endl;

			//---------------------------------------------//

		} // end of the loop over the detailed xsec variations			
*/
	} // end of the loop over the plots

	//---------------------------------------------//

} // End of the program
