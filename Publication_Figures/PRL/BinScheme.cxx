#include <TFile.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

//----------------------------------------//

TString PrintMultipleTimes(int times, TString string) {

	TString MergeTStrings = "";

    for (int itime = 0; itime < times; itime++) {

		MergeTStrings += string;

	} 

	return MergeTStrings;
 
}

//----------------------------------------//

void BinScheme() {

	//----------------------------------------//

	vector<TString> PlotNames;
	PlotNames.push_back("SerialDeltaPT_DeltaAlphaTPlot");
	PlotNames.push_back("SerialDeltaAlphaT_DeltaPTPlot");
	PlotNames.push_back("SerialDeltaPtx_DeltaPtyPlot");		

	const int NPlots = PlotNames.size();
	
	std::vector< std::vector<double> > binning;

	//----------------------------------------//

	// Open the file that contains all the xsecs

	TString XSecFileName = "../../myXSec/v08_00_00_52/GenXSec/All_XSecs_Combined_v08_00_00_52.root";
	TFile* fXSec = new TFile(XSecFileName,"readonly");

	//----------------------------------------//

	// Data release

	TString TxtName = "/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance_PRL_Rename/BinScheme.txt";
	ofstream myTxtFile;
	myTxtFile.open(TxtName);		

	//----------------------------------------//		

	// Loop over the plots

	for (int iplot = 0; iplot < NPlots; iplot ++) {

		//----------------------------------------//
	
		if (PlotNames[iplot] == "SerialDeltaPT_DeltaAlphaTPlot") { binning = TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices; }
		if (PlotNames[iplot] == "SerialDeltaAlphaT_DeltaPTPlot") { binning = TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices; }
		if (PlotNames[iplot] == "SerialDeltaPtx_DeltaPtyPlot") { binning = TwoDArrayNBinsDeltaPtxInDeltaPtySlices; }				

		//----------------------------------------//
		
		int slices = binning.size();
		vector<int> binsinslice;
		
		for (int islice = 0; islice < slices; islice++) {
					
			binsinslice.push_back(binning[islice].size()-1);
		
		}

		//----------------------------------------//

		TH1D* BeamOnFullUnc = (TH1D*)( fXSec->Get("FullUnc_" + PlotNames[iplot]) );	
		int NBins = BeamOnFullUnc->GetXaxis()->GetNbins();
					
		myTxtFile << PlotNames[iplot] << endl << endl;		
		myTxtFile << "Bin # & Slice & Low edge & High edge" << endl << endl;
		
		//cout << PlotNames[iplot] << endl << endl;		
		//cout << "Bin # & Slice & Low edge & High edge" << endl << endl;			
		
		int slice = 0;	
		int bincounter = 0;	

		for (int ibin = 1; ibin <= NBins; ibin++) {

			double BinLow = binning[slice][bincounter];		
			double BinHigh = binning[slice][bincounter+1];				

			myTxtFile << ibin << " & " << LatexLabel[ MapUncorCor[PlotNames[iplot] + "_" + TString(std::to_string(slice) ) ] ];
			myTxtFile << " & " << BinLow << " & " << BinHigh << endl;
			
			//cout << ibin << " & " << LatexLabel[ MapUncorCor[PlotNames[iplot] + "_" + TString(std::to_string(slice) ) ] ];
			//cout << " & " << BinLow << " & " << BinHigh << endl;			
				
			bincounter++;
			if (bincounter == binsinslice[slice]) { slice++; bincounter = 0;}

		}

		myTxtFile << endl << endl;
		//cout << endl << endl;		

		//----------------------------------------//

	} // End of the loop over the plots

} // End of the program 
