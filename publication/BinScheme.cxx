#include <TFile.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../../myClasses/Constants.h"

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
	PlotNames.push_back("SerialThetaVis_ECalPlot");
	PlotNames.push_back("SerialThetaVis_DeltaPnPlot");
	PlotNames.push_back("SerialThetaVis_PMissPlot");
	
	const int NPlots = PlotNames.size();
	
	std::vector< std::vector<double> > binning;

	//----------------------------------------//

	// Open the file that contains all the xsecs

	TString XSecFileName = PathToExtractedXSec+"/WienerSVD_ExtractedXSec_Overlay9_Combined_"+UBCodeVersion+".root";
	TFile* fXSec = new TFile(XSecFileName,"readonly");

	//----------------------------------------//

	// Data release

	TString TxtName = "./BinScheme.txt";
	ofstream myTxtFile;
	myTxtFile.open(TxtName);		

	//----------------------------------------//		

	// Loop over the plots

	for (int iplot = 0; iplot < NPlots; iplot ++) {

		//----------------------------------------//
	
		if (PlotNames[iplot] == "SerialThetaVis_ECalPlot") { binning = TwoDArrayNBinsThetaVisInECalSlices; }
		if (PlotNames[iplot] == "SerialThetaVis_DeltaPnPlot") { binning = TwoDArrayNBinsThetaVisInDeltaPnSlices; }
		if (PlotNames[iplot] == "SerialThetaVis_PMissPlot") { binning = TwoDArrayNBinsThetaVisInPMissSlices; }

		//----------------------------------------//
		
		int slices = binning.size();
		vector<int> binsinslice;
		
		for (int islice = 0; islice < slices; islice++) {
					
			binsinslice.push_back(binning[islice].size()-1);
		
		}

		//----------------------------------------//

		TH1D* BeamOnFullUnc = (TH1D*)( fXSec->Get("RecoFullUnc" + PlotNames[iplot]) );	
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
