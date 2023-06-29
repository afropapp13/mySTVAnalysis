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

void Correlations_XSecSignalBkgUnc() {

	//---------------------------------------------//

	GlobalSettings();
	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	TGaxis::SetMaxDigits(3);

	std::cout << std::setprecision(3) << std::fixed;	

	TString ExactFileLocation = PathToFiles+CutExtension;

	//---------------------------------------------//

	vector<TString> PlotNames;
	PlotNames.push_back("DeltaPTPlot"); 

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
	//Vars.push_back("XSecShape_CCMEC_UBGenie");
	//Vars.push_back("All_UBGenie");

	int nvars = Vars.size();							

	//---------------------------------------------//	

	// Loop over the plots

	for (int iplot = 0; iplot < nplots; iplot++) {						

		//---------------------------------------------//									

		// Loop over the xsec variations

		for (int ivar = 0; ivar < nvars; ivar++) {

			//---------------------------------------------//				

			TString ivar_cov_xsecfile_name = MigrationMatrixPath+"IndividualWienerSVD_" + Vars[ivar] + "_CovarianceMatrices_Overlay9_Combined_"+UBCodeVersion+".root";
			cout << ivar_cov_xsecfile_name << endl;
			TFile* ivar_nomcov_xsecfile = TFile::Open(ivar_cov_xsecfile_name,"readonly");
			TH2D* ivar_nomfraccovmatrix = (TH2D*)(ivar_nomcov_xsecfile->Get(Vars[ivar] + "_FracCovariance_"+PlotNames[iplot]+"_Combined"));		

			ivar_nomfraccovmatrix->Draw("coltz");			

			//---------------------------------------------//

		} // end of the loop over the xsec variations

	}

	//---------------------------------------------//

} // End of the program
