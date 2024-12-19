#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
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

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../../myClasses/Constants.h"
#include "../../myClasses/myFunctions.cpp"

using namespace std;
using namespace Constants;

static std::map<TString,TString> Mapping = {

		{ "ThetaVisPlot", "ThetaVisPlot" }

};

//----------------------------------------//

void thetavis_1d_print_xsecs() {

	//----------------------------------------//

	int DecimalAccuracy = 2;

	//----------------------------------------//

	vector<TString> PlotNames; vector<TString> Units;
	PlotNames.push_back("ThetaVisPlot"); Units.push_back("[$10^{-38}\\frac{\\mathrm{cm}^{2}}{\\mathrm{deg}\\,\\mathrm{Ar}}$]");

	const int NPlots = PlotNames.size();

	//----------------------------------------//

	// Open the file that contains all the xsecs

	TString XSecFileName = PathToExtractedXSec+"/WienerSVD_ExtractedXSec_Overlay9_Combined_"+UBCodeVersion+".root";
	TFile* fXSec = new TFile(XSecFileName,"readonly");

	//----------------------------------------//

	// Data release

	TString TxtName = "thetavis_1d_print_xsecs.tex";
	ofstream myTxtFile;
	myTxtFile.open(TxtName);

	//----------------------------------------//		

	// Loop over the plots

	for (int iplot = 0; iplot < NPlots; iplot ++) {

			//----------------------------------------//					
			
			TH1D* BeamOnFullUnc = (TH1D*)( fXSec->Get("Reco" + PlotNames[iplot]) );	
			int NBins = BeamOnFullUnc->GetXaxis()->GetNbins();

			TH2D* Ac = (TH2D*)fXSec->Get("Ac"+PlotNames[iplot]);	
			TH2D* Cov = (TH2D*)fXSec->Get("UnfCov"+PlotNames[iplot]);
			
			TString LatexLabelString = "$\\mathrm{"+LatexLabel[ Mapping[ PlotNames[iplot] ] ]+"}$";
			LatexLabelString.ReplaceAll("#","\\").ReplaceAll(" ","\\,");
	
			myTxtFile << "\\begin{table}[H]" << endl;
			myTxtFile << "\\raggedright" << endl;	
			myTxtFile << "\\begin{adjustbox}{width=\\textwidth}" << endl;						
			myTxtFile << "\\small" << endl;
			myTxtFile << "\\begin{tabular}{ |c|c|c|c|c| }" << endl;	
			myTxtFile << "\\hline" << endl;						
			myTxtFile << "\\multicolumn{5}{|c|}{Cross Section $\\theta_{\\mathrm{vis}}$, " << LatexLabelString << "} \\\\" << endl;
			myTxtFile << "\\hline" << endl;
			myTxtFile << "\\hline" << endl;			
			myTxtFile << "Bin \\# & Low edge [deg] & High edge [deg] & Cross Section " << Units[iplot] <<" & Uncertainty " << Units[iplot] << " \\\\" << endl;			
			myTxtFile << "\\hline" << endl;
			myTxtFile << "\\hline" << endl;	

			for (int ibin = 1; ibin <= NBins; ibin++) {

				double BinLow = BeamOnFullUnc->GetBinLowEdge(ibin);
				double BinWidth = BeamOnFullUnc->GetBinWidth(ibin);		
				double BinHigh = BinLow + BinWidth;	
				double BinValue = BeamOnFullUnc->GetBinContent(ibin);
				double BinError = BeamOnFullUnc->GetBinError(ibin);				

				myTxtFile << ibin << std::setprecision(4) << " & " << BinLow << " & " << BinHigh << std::setprecision(8) << " & " << BinValue << " & " <<  BinError << "\\\\" << endl;

			}
				
			myTxtFile << "\\hline" << endl;			
			myTxtFile << "\\end{tabular}" << endl;
			myTxtFile << "\\end{adjustbox}" << endl;		
			myTxtFile << "\\end{table}" << endl;				
			myTxtFile << endl << endl;

			//----------------------------------------//			

		} // End of the loop over the plots

		//----------------------------------------//					

} // End of the program 
