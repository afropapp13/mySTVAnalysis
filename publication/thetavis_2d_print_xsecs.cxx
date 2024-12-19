#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../../myClasses/myFunctions.cpp"
#include "../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "../../myClasses/Util.h"

//----------------------------------------//

void thetavis_2d_print_xsecs() {
	
	//----------------------------------------//

	Tools tools;

	int DecimalAccuracy = 2;

	//----------------------------------------//

	vector<TString> PlotNames;
	vector<TString> Units;
	vector< vector<double> > SliceDiscriminators;
	vector< vector< vector<double> > > SliceBinning;

	//----------------------------------------//

	// Data release

	TString TxtName = "thetavis_2d_print_xsecs.tex";
	ofstream myTxtFile;
	myTxtFile.open(TxtName);

	//----------------------------------------//		

	// 2D analysis

	PlotNames.push_back("ThetaVis_ECalPlot"); Units.push_back("[$10^{-38}\\mathrm{\\frac{cm^{2}}{deg\\,GeV\\,Ar}}$]"); 
	PlotNames.push_back("ThetaVis_DeltaPnPlot");  Units.push_back("[$10^{-38}\\mathrm{\\frac{cm^{2}}{deg\\,(GeV/c)\\,Ar}}$]");
	PlotNames.push_back("ThetaVis_PMissPlot");  Units.push_back("[$10^{-38}\\mathrm{\\frac{cm^{2}}{deg\\,(GeV/c)\\,Ar}}$]");
	
	const int N1DPlots = PlotNames.size();

	//----------------------------------------//

	// vector of plots	
	vector<TH1D*> PlotsFullUncReco; PlotsFullUncReco.clear(); PlotsFullUncReco.resize(N1DPlots);

	// vector of vectors for vector slices
	vector<vector<TH1D*> > BeamOnFullUnc; BeamOnFullUnc.clear(); BeamOnFullUnc.resize(N1DPlots);

	TString FileSampleName = PathToExtractedXSec+"/WienerSVD_ExtractedXSec_Overlay9_Combined_"+UBCodeVersion+".root"; 
	TFile* FileSample = TFile::Open(FileSampleName,"readonly"); 
	
	for (int iplot = 0; iplot < N1DPlots; iplot ++) {

		//----------------------------------------//
		
		TH1D* histFullUncReco = (TH1D*)(FileSample->Get("RecoFullUncSerial"+PlotNames[iplot]));
		PlotsFullUncReco[iplot] = histFullUncReco;	
		
		TH2D* Ac = (TH2D*)FileSample->Get("AcSerial"+PlotNames[iplot]);

		TH2D* Cov = (TH2D*)FileSample->Get("UnfCovSerial"+PlotNames[iplot]);	
		Cov->Scale(1./TMath::Power(MultiDimScaleFactor["Serial" + PlotNames[iplot]],2.)); // includes scaling factor for multi dimensional analysis

		//----------------------------------------//

		// Setting up the relevant discriminators

		SliceDiscriminators.clear();
		SliceBinning.clear();

		if (PlotNames[iplot] == "ThetaVis_ECalPlot") {

			SliceDiscriminators.push_back(TwoDArrayNBinsECal); 
			SliceBinning.push_back(TwoDArrayNBinsThetaVisInECalSlices);

		}

		if (PlotNames[iplot] == "ThetaVis_DeltaPnPlot") {

			SliceDiscriminators.push_back(TwoDArrayNBinsDeltaPn); 
			SliceBinning.push_back(TwoDArrayNBinsThetaVisInDeltaPnSlices);

		}

		if (PlotNames[iplot] == "ThetaVis_PMissPlot") {

			SliceDiscriminators.push_back(TwoDArrayNBinsPMiss); 
			SliceBinning.push_back(TwoDArrayNBinsThetaVisInPMissSlices);

		}

		//----------------------------------------//

		// The covariance matrix needs to be scaled by the 2D bin width

		TH2D* CovClone = (TH2D*)(Cov->Clone()); 
		int n = Cov->GetXaxis()->GetNbins();

		for (int ix = 1; ix <= n; ix++) {

			for (int iy = 1; iy <= n; iy++) {

				double WidthX = Cov->GetXaxis()->GetBinWidth(ix);
				double WidthY = Cov->GetYaxis()->GetBinWidth(iy);

				double TwoDWidth = WidthX * WidthY;

				double BinContent = Cov->GetBinContent(ix,iy);
				// Division by bin width already included
				// Still need to include scaling due to slice range
				// That is done in Tools::Get2DHistoBins
				double NewBinContent = BinContent/TwoDWidth;
				CovClone->SetBinContent(ix,iy,NewBinContent);

			}					

		}	

		//------------------------------------//

		// Number of N-dimensional slices

		int NSlices = 1;
		vector<double> SerialVectorRanges;
		vector<int> SerialVectorBins;
		vector<int> SerialVectorLowBin;
		vector<int> SerialVectorHighBin;						

		int BinCounter = 1;

		// vector< vector<double> > SliceDiscriminators;
		// 1st index = how many disciminators
		// 2nd index = values / ranges of discriminators

		// How many discriminators do we have ? e.g. cos theta mu, Pmu et al
		for (int islice = 0; islice < (int)(SliceDiscriminators.size()); islice++) { 
				
			// For a given discriminator, how many slices do we have ? SliceDiscrimSize - 1
			int SliceDiscrimSize = SliceDiscriminators.at(islice).size()-1;
			NSlices *= SliceDiscrimSize; 

			for (int iSliceDiscrimSize = 0; iSliceDiscrimSize < SliceDiscrimSize; iSliceDiscrimSize++) {

				// Accessing the vector<double> with the bin ranges
				int SliceDiscrimValue = SliceBinning.at(0).at(iSliceDiscrimSize).size();

				// Storing the number of bins for a specific slice					
				SerialVectorBins.push_back(SliceDiscrimValue-1);
				for (int iBin = 0; iBin < SliceDiscrimValue; iBin++) {

					double BinValue = SliceBinning.at(0).at(iSliceDiscrimSize).at(iBin);
					// First bin number for a given slice
					if (iBin == 0) { SerialVectorLowBin.push_back(BinCounter); }
					// Last bin number for a given slice
					if (iBin == SliceDiscrimValue-2) { SerialVectorHighBin.push_back(BinCounter); }	

					// Storing the binning for a specific slice
					SerialVectorRanges.push_back(BinValue);

					// Increase the global bin counter
					// But not for the last bin
					if (iBin != SliceDiscrimValue-1) { BinCounter++; }

				} // End of the loop over the bins of a given slice
				//cout << "End of the loop over the bins of a given slice" << endl;

			} // End of the loop over the slices for a given discriminator

		} // End of the loop over the number of discriminators	

		//------------------------------------//
		       
		int StartIndex = 0;
		int BinStartIndex = 0;			

		//------------------------------------//

		// Loop over the N-dimensional slices

		BeamOnFullUnc[iplot].resize( NSlices );
		
		for (int NDimSlice = 0; NDimSlice < NSlices; NDimSlice++) {	

			//------------------------------------//

			TString NameCopy = PlotNames[iplot];

			NameCopy.ReplaceAll("unf_","");
			NameCopy.ReplaceAll("TrueUnf_","");
			NameCopy.ReplaceAll("True_","");
			NameCopy.ReplaceAll("True","");
			NameCopy.ReplaceAll("NoSmearAlt","");			

			NameCopy.ReplaceAll("_Run1","");
			NameCopy.ReplaceAll("_Run2","");
			NameCopy.ReplaceAll("_Run3","");
			NameCopy.ReplaceAll("_Run4","");
			NameCopy.ReplaceAll("_Run4a","");	
			NameCopy.ReplaceAll("_Run5","");
			NameCopy.ReplaceAll("_Combined","");	

			NameCopy = "Serial" + NameCopy + "_" + TString(std::to_string(NDimSlice));	

			//------------------------------------//		
				
			// Get the number of bins and the bin ranges for the specific slice	

			int SliceNBins = SerialVectorBins.at(NDimSlice);
			std::vector<double> SerialSliceBinning;		

			for (int iBin = 0; iBin < SliceNBins+1; iBin++) { 

				double value = SerialVectorRanges.at(StartIndex+iBin);
				SerialSliceBinning.push_back(value);
			
			} // End of the number of bins and the bin ranges declaration	
			
			//------------------------------------//

			// Slice out the xsec plot
			BeamOnFullUnc[iplot][NDimSlice] = tools.GetHistoBins(PlotsFullUncReco[iplot],SerialVectorLowBin.at(NDimSlice),\
                                                          SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"FullUnc");

			// Start printing out the xsecs

			TString LatexLabelString = "$\\mathrm{"+LatexLabel[ MapUncorCor[ NameCopy ] ]+"}$";
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

			int nbins = BeamOnFullUnc[iplot][NDimSlice]->GetXaxis()->GetNbins();

			for (int ibin = 1; ibin <= nbins; ibin++) {

				double BinLow = BeamOnFullUnc[iplot][NDimSlice]->GetBinLowEdge(ibin);
				double BinWidth = BeamOnFullUnc[iplot][NDimSlice]->GetBinWidth(ibin);		
				double BinHigh = BinLow + BinWidth;	
				double BinValue = BeamOnFullUnc[iplot][NDimSlice]->GetBinContent(ibin);
				double BinError = BeamOnFullUnc[iplot][NDimSlice]->GetBinError(ibin);				

				myTxtFile << ibin << std::setprecision(4) << " & " << BinLow << " & " << BinHigh << std::setprecision(8) << " & " << BinValue << " & " <<  BinError << "\\\\" << endl;

			}
				
			myTxtFile << "\\hline" << endl;			
			myTxtFile << "\\end{tabular}" << endl;
			myTxtFile << "\\end{adjustbox}" << endl;		
			myTxtFile << "\\end{table}" << endl;				
			myTxtFile << endl << endl;

			//------------------------------------//
				
			// Update the starting index to move to the next slice

			StartIndex += (SliceNBins+1);
			BinStartIndex += SliceNBins;
				
			//------------------------------------//
				
		} // End of the loop over the discriminators slices
			
		//----------------------------------------//

	} // End of the loop over the plots

} // End of the program 
