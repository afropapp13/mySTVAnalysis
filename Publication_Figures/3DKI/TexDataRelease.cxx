#include <TFile.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

static std::map<TString,TString> Mapping = {

		{ "DeltaPnPlot", "DeltaPnPlot" },
		{ "DeltaAlpha3DqPlot", "DeltaAlpha3DqPlot" },
		{ "DeltaPhi3DPlot", "DeltaPhi3DPlot" },		
		{ "DeltaPnPerpPlot", "DeltaPnPerpPlot" },
		{ "DeltaPnPerpxPlot", "DeltaPnPerpxPlot" },
		{ "DeltaPnPerpyPlot", "DeltaPnPerpyPlot" },					
		{ "DeltaPnParPlot", "DeltaPnParPlot" },					
		{ "SerialDeltaPn_DeltaAlpha3DqPlot_0", "DeltaPn_DeltaAlpha3Dq_0_00To45_00Plot" },
		{ "SerialDeltaPn_DeltaAlpha3DqPlot_1", "DeltaPn_DeltaAlpha3Dq_45_00To90_00Plot" },
		{ "SerialDeltaPn_DeltaAlpha3DqPlot_2", "DeltaPn_DeltaAlpha3Dq_90_00To135_00Plot" },
		{ "SerialDeltaPn_DeltaAlpha3DqPlot_3", "DeltaPn_DeltaAlpha3Dq_135_00To180_00Plot" },
		{ "SerialDeltaAlpha3Dq_DeltaPnPlot_0", "DeltaAlpha3Dq_DeltaPn_0_00To0_20Plot" },		
		{ "SerialDeltaAlpha3Dq_DeltaPnPlot_1", "DeltaAlpha3Dq_DeltaPn_0_20To0_40Plot" },
		{ "SerialDeltaAlpha3Dq_DeltaPnPlot_2", "DeltaAlpha3Dq_DeltaPn_0_40To1_00Plot" }
};

//----------------------------------------//

void TexDataRelease() {

	//----------------------------------------//

	vector<TString> PlotNames; vector<TString> Units; vector<TString> Var;
	PlotNames.push_back("DeltaPnPlot"); Units.push_back("GeV/\\textit{c}"); Var.push_back("$p_{n}$");
	PlotNames.push_back("DeltaAlpha3DqPlot"); Units.push_back("deg"); Var.push_back("$\\alpha_{3D}$");
	PlotNames.push_back("DeltaPhi3DPlot"); Units.push_back("deg"); Var.push_back("$\\phi_{3D}$");				
	PlotNames.push_back("DeltaPnPerpPlot"); Units.push_back("GeV/\\textit{c}"); Var.push_back("$p_{n\\perp}$");
	PlotNames.push_back("DeltaPnPerpxPlot"); Units.push_back("GeV/\\textit{c}"); Var.push_back("$p_{n\\perp,x}$");
	PlotNames.push_back("DeltaPnPerpyPlot"); Units.push_back("GeV/\\textit{c}"); Var.push_back("$p_{n\\perp,y}$");
	PlotNames.push_back("DeltaPnParPlot"); Units.push_back("GeV/\\textit{c}"); Var.push_back("$p_{n\\parallel}$");		
	PlotNames.push_back("SerialDeltaPn_DeltaAlpha3DqPlot_0"); Units.push_back("GeV/\\textit{c}"); Var.push_back("$p_{n}$");
	PlotNames.push_back("SerialDeltaPn_DeltaAlpha3DqPlot_3"); Units.push_back("GeV/\\textit{c}"); Var.push_back("$p_{n}$");
	PlotNames.push_back("SerialDeltaAlpha3Dq_DeltaPnPlot_0"); Units.push_back("deg"); Var.push_back("$\\alpha_{3D}$");
	PlotNames.push_back("SerialDeltaAlpha3Dq_DeltaPnPlot_2"); Units.push_back("deg"); Var.push_back("$\\alpha_{3D}$");		

	const int NPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << NPlots << endl;

	//----------------------------------------//

	TString Runs = "Combined";

	//----------------------------------------//

	// Open the file that contains all the xsecs

	TString XSecFileName = "../../myXSec/" + UBCodeVersion + "/GenXSec/All_XSecs_Combined_" + UBCodeVersion + ".root";
	TFile* fXSec = new TFile(XSecFileName,"readonly");

	//----------------------------------------//

	// Data release

	TString TxtName = "/home/afroditi/Dropbox/Apps/Overleaf/General Imbalance Variables for Measuring Nuclear Effects and Demonstration with MicroBooNE Data/XSec_DataRelease.tex";
	ofstream myTxtFile;
	myTxtFile.open(TxtName);		

	//----------------------------------------//		

	// Loop over the plots

	for (int iplot = 0; iplot < NPlots; iplot ++) {													

		//----------------------------------------//		

		// Data release

		TH1D* BeamOnFullUnc = (TH1D*)( fXSec->Get("FullUnc_" + PlotNames[iplot]) );	
		int NBins = BeamOnFullUnc->GetXaxis()->GetNbins();

		TString LatexLabelString = "$"+LatexLabel[ Mapping[ PlotNames[iplot] ] ]+"$";
		LatexLabelString.ReplaceAll("#","\\").ReplaceAll(" ","\\,");
		myTxtFile << "\\begin{table}[H]" << endl;
		myTxtFile << "\\raggedright" << endl;	
		myTxtFile << "\\begin{adjustbox}{width=\\textwidth}" << endl;						
		myTxtFile << "\\small" << endl;
		myTxtFile << "\\begin{tabular}{ |c|c|c|c|c| }" << endl;	
		myTxtFile << "\\hline" << endl;						
		myTxtFile << "\\multicolumn{5}{|c|}{Cross Section " + Var[iplot] + ", " << LatexLabelString << "} \\\\" << endl;
		myTxtFile << "\\hline" << endl;
		myTxtFile << "\\hline" << endl;			
		myTxtFile << "Bin \\# & Low edge [" << Units[iplot] << "] & High edge [" << Units[iplot] << "] & Cross Section [$10^{-38}\\frac{cm^{2}}{(" << Units[iplot] << ")\\,^{40}Ar}$]" << " & Uncertainty [$10^{-38}\\frac{cm^{2}}{(" << Units[iplot] << ")\\,^{40}Ar}$] \\\\" << endl;			
		myTxtFile << "\\hline" << endl;
		myTxtFile << "\\hline" << endl;			

		for (int ibin = 1; ibin <= NBins; ibin++) {

			double BinLow = BeamOnFullUnc->GetBinLowEdge(ibin);
			double BinWidth = BeamOnFullUnc->GetBinWidth(ibin);		
			double BinHigh = BinLow + BinWidth;	
			double BinValue = BeamOnFullUnc->GetBinContent(ibin);
			double BinError = BeamOnFullUnc->GetBinError(ibin);				

			myTxtFile << ibin << std::setprecision(4) << " & " << BinLow << " & " << BinHigh << std::setprecision(5) << " & " << BinValue << " & " <<  BinError << "\\\\" << endl;

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