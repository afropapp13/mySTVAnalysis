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

void WienerSVD_XSecVars() {

	//---------------------------------------------//

	GlobalSettings();
	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	TGaxis::SetMaxDigits(3);

	std::cout << std::setprecision(3) << std::fixed;	

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

	// Unisims

	// DetVars.push_back("UnShortAxFFCCQEshape_UBGenie");
	// DetVars.push_back("UnShortDecayAngMEC_UBGenie");
	// DetVars.push_back("UnShortRPA_CCQE_UBGenie");
	// DetVars.push_back("UnShortTheta_Delta2Npi_UBGenie");
	// DetVars.push_back("UnShortVecFFCCQEshape_UBGenie");
	// DetVars.push_back("UnShortXSecShape_CCMEC_UBGenie");									

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

		TH1D* xsecfullunc = (TH1D*)(xsecfile->Get("XSecReco"+PlotNames[iplot]));

		// Total xsec syst		

		TH2D* nomcovmatrix = (TH2D*)(nomcov_xsecfile->Get("XSec_Covariance_"+PlotNames[iplot]+"_Combined"));
		TH2D* nomfraccovmatrix = (TH2D*)(nomcov_xsecfile->Get("XSec_FracCovariance_"+PlotNames[iplot]+"_Combined"));						

		//---------------------------------------------//

		double calcxsecunc = TMath::Sqrt(nomfraccovmatrix->GetBinContent(1,1)) * 100.;		

		double entry = xsecfullunc->GetBinContent(1);
		double error = xsecfullunc->GetBinError(1);		
		double ratio = error / entry * 100.;

		// Sacling factor due to efficiency corrections
		double sf = ratio / calcxsecunc;	

		//---------------------------------------------//	

		double 	All_UBGenieValue = 0.;
		double 	TotalUnfXSecFracUnc = 0.;
		double 	TotalEventsXSecFracUnc = 0.;										

		// Loop over the xsec variations

		for (int ivar = 0; ivar < nvars; ivar++) {

			//---------------------------------------------//				

			TString ivar_cov_xsecfile_name = MigrationMatrixPath+"IndividualWienerSVD_" + Vars[ivar] + "_CovarianceMatrices_Overlay9_Combined_"+UBCodeVersion+".root";
			TFile* ivar_nomcov_xsecfile = TFile::Open(ivar_cov_xsecfile_name,"readonly");
			TH2D* ivar_nomfraccovmatrix = (TH2D*)(ivar_nomcov_xsecfile->Get(Vars[ivar] + "_FracCovariance_"+PlotNames[iplot]+"_Combined"));					

			double localxsecevents = ivar_nomfraccovmatrix->GetBinContent(1,1);
			double localxsecunf = localxsecevents * TMath::Power(sf,2.);	

			cout << Vars[ivar] << ", fractional contribution = " << TMath::Sqrt( localxsecunf ) * 100. << " %"<< endl;			

			TotalEventsXSecFracUnc += localxsecevents;
			TotalUnfXSecFracUnc += localxsecunf;			

			if (Vars[ivar] == "All_UBGenie") { All_UBGenieValue = localxsecevents; }

			//---------------------------------------------//

		} // end of the loop over the xsec variations

		double summedxsecunc = TMath::Sqrt(TotalEventsXSecFracUnc) * 100.;

		cout << endl << "Sanity check" << endl;
		cout << "Total summed xsec fractional uncertainty = " << summedxsecunc << " %" << endl;
		cout << "Total xsec fractional uncertainty (events) = " << calcxsecunc << " %" << endl;	
		cout << "Total xsec fractional uncertainty (unf) = " << ratio << " %" << endl << endl;							

		//---------------------------------------------//	

		double DetAll_UBGenieValue = 0.;	

		for (int idetvar = 0; idetvar < ndetvars; idetvar++) {

			if ( !(string(DetVars[idetvar]).find("UnShort") != std::string::npos) ) {

				TString ivar_cov_xsecfile_name = MigrationMatrixPath+"IndividualWienerSVD_" + DetVars[idetvar] + "_CovarianceMatrices_Overlay9_Run1_DecompXSecUnc_"+UBCodeVersion+".root";
				TFile* ivar_nomcov_xsecfile = TFile::Open(ivar_cov_xsecfile_name,"readonly");
				TH2D* detfraccovmatrix = (TH2D*)(ivar_nomcov_xsecfile->Get(DetVars[idetvar] + "_FracCovariance_"+PlotNames[iplot]+"_Run1_DecompXSecUnc"));	

				double localxsecevents = detfraccovmatrix->GetBinContent(1,1);
				DetAll_UBGenieValue += localxsecevents;

			}

		}		

		double sf_sample = All_UBGenieValue	/ DetAll_UBGenieValue;

		//---------------------------------------------//						

		// Loop over the detailed xsec variations

		double sumxsecunc = 0.;

		for (int idetvar = 0; idetvar < ndetvars; idetvar++) {

			TString ivar_cov_xsecfile_name = MigrationMatrixPath+"IndividualWienerSVD_" + DetVars[idetvar] + "_CovarianceMatrices_Overlay9_Run1_DecompXSecUnc_"+UBCodeVersion+".root";
			TFile* ivar_nomcov_xsecfile = TFile::Open(ivar_cov_xsecfile_name,"readonly");			
			TH2D* detfraccovmatrix = (TH2D*)(ivar_nomcov_xsecfile->Get(DetVars[idetvar] + "_FracCovariance_"+PlotNames[iplot]+"_Run1_DecompXSecUnc"));							

			double localfracunc = detfraccovmatrix->GetBinContent(1,1) * TMath::Power(sf,2.) * TMath::Power(sf_sample,2.);
			if ( !(string(DetVars[idetvar]).find("UnShort") != std::string::npos) ) { sumxsecunc += localfracunc; }

			cout << DetVars[idetvar] << ", fractional contribution = " << TMath::Sqrt( detfraccovmatrix->GetBinContent(1,1) ) * 100. << " %"<< endl;						

		} // end of the loop over the detailed xsec variations		

		double detailedcalcunc = TMath::Sqrt(sumxsecunc) * 100.;
		cout << endl << "Total All_UBGenie detailed xsec unc = " << detailedcalcunc << " %" << endl;	

	} // end of the loop over the plots

	//---------------------------------------------//

} // End of the program
