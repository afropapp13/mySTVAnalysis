#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TEfficiency.h>
#include <TMath.h>
#include <TLatex.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

using namespace std;

// --------------------------------------------------------------------------------------------------------------------------------------------

void PrintFinalUnc() {

	std::cout << std::fixed << std::setprecision(2);

	// --------------------------------------------------------------------------------------------------------------------------------------------

	// Numbers reported in %
	// Intentionally leaving out the POT (2%) & NTarget (1%) contrinution
	// Given that they are overall normalization uncertainties

	// --------------------------------------------------------------------------------------------------------------------------------------------

	// Run 1

	double Run1Stat = 14.97;
	double Run1POT = 0.; // 2%
	double Run1NTarget = 0.; // 1%
	double Run1LY = 5.87;
	double Run1TPC = 16.96;
	double Run1G4 = 0.;
	double Run1Genie = 0.;
	double Run1Flux = 0.;
	double Run1Dirt = 0.91;

	double Run1StatMC = 28.86;
	double Run1POTSmEff = 0.; // 2%
	double Run1NTargetSmEff = 0.; // 1%
	double Run1LYSmEff = 0.;
	double Run1TPCSmEff = 0.;
	double Run1G4SmEff = 0.;
	double Run1GenieSmEff = 0.;
	double Run1FluxSmEff = 0.;

	// --------------------------------------------------------------------------------------------------------------------------------------------

	// Run 3

	double Run3Stat = 23.90;
	double Run3POT = 0.; // 2%
	double Run3NTarget = 0.; // 1%
	double Run3LY = 1.08;
	double Run3TPC = 1.36;
	double Run3G4 = 0.;
	double Run3Genie = 0.;
	double Run3Flux = 0.;
	double Run3Dirt = 0.66;

	double Run3StatMC = 0.;
	double Run3POTSmEff = 0.; // 2%
	double Run3NTargetSmEff = 0.; // 1%
	double Run3LYSmEff = 0.;
	double Run3TPCSmEff = 0.;
	double Run3G4SmEff = 0.;
	double Run3GenieSmEff = 0.;
	double Run3FluxSmEff = 0.;

	// --------------------------------------------------------------------------------------------------------------------------------------------

	// Run 1: Merge the relevant uncertainties by adding them in quadrature
	// Print the table to be added in the AN

	double Run1TotalStat = TMath::Sqrt( TMath::Power(Run1Stat,2.) + TMath::Power(Run1StatMC,2.) );
	double Run1TotalPOT = TMath::Sqrt( TMath::Power(Run1POT,2.) + TMath::Power(Run1POTSmEff,2.) );
	double Run1TotalNTarget = TMath::Sqrt( TMath::Power(Run1NTarget,2.) + TMath::Power(Run1NTargetSmEff,2.) );
	double Run1TotalDet = TMath::Sqrt( TMath::Power(Run1LY,2.) + TMath::Power(Run1TPC,2.) + TMath::Power(Run1LYSmEff,2.) + TMath::Power(Run1TPCSmEff,2.) );
	double Run1TotalG4 = TMath::Sqrt( TMath::Power(Run1G4,2.) + TMath::Power(Run1G4SmEff,2.));
	double Run1TotalGenie = TMath::Sqrt( TMath::Power(Run1Genie,2.) + TMath::Power(Run1GenieSmEff,2.));
	double Run1TotalFlux = TMath::Sqrt( TMath::Power(Run1Flux,2.) + TMath::Power(Run1FluxSmEff,2.));
	double Run1TotalDirt = TMath::Sqrt( TMath::Power(Run1Dirt,2.));

	double Run1Total = TMath::Sqrt( TMath::Power(Run1TotalStat,2.)+TMath::Power(Run1TotalDet,2.)+TMath::Power(Run1TotalG4,2.)+TMath::Power(Run1TotalGenie,2.)+TMath::Power(Run1TotalFlux,2.)+TMath::Power(Run1TotalDirt,2.));

	cout << endl;
	cout << "Stat & " << Run1TotalStat << "\\\\" << endl;
	cout << "Detector & " << Run1TotalDet << "\\\\" << endl;
	cout << "G4 & " << Run1TotalG4 << "\\\\" << endl;
	cout << "XSec & " << Run1TotalGenie << "\\\\" << endl;
	cout << "Flux & " << Run1TotalFlux << "\\\\" << endl;
	cout << "Dirt & " << Run1TotalDirt << "\\\\" << endl;

	cout << "\\hline \\hline" << endl; 
	cout << "Total & " << Run1Total << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << endl << endl;

	// --------------------------------------------------------------------------------------------------------------------------------------------

	// Run 3: Merge the relevant uncertainties by adding them in quadrature
	// Print the table to be added in the AN

	double Run3TotalStat = TMath::Sqrt( TMath::Power(Run3Stat,2.) + TMath::Power(Run3StatMC,2.) );
	double Run3TotalPOT = TMath::Sqrt( TMath::Power(Run3POT,2.) + TMath::Power(Run3POTSmEff,2.) );
	double Run3TotalNTarget = TMath::Sqrt( TMath::Power(Run3NTarget,2.) + TMath::Power(Run3NTargetSmEff,2.) );
	double Run3TotalDet = TMath::Sqrt( TMath::Power(Run3LY,2.) + TMath::Power(Run3TPC,2.) + TMath::Power(Run3LYSmEff,2.) + TMath::Power(Run3TPCSmEff,2.) );
	double Run3TotalG4 = TMath::Sqrt( TMath::Power(Run3G4,2.) + TMath::Power(Run3G4SmEff,2.));
	double Run3TotalGenie = TMath::Sqrt( TMath::Power(Run3Genie,2.) + TMath::Power(Run3GenieSmEff,2.));
	double Run3TotalFlux = TMath::Sqrt( TMath::Power(Run3Flux,2.) + TMath::Power(Run3FluxSmEff,2.));
	double Run3TotalDirt = TMath::Sqrt( TMath::Power(Run3Dirt,2.));

	double Run3Total = TMath::Sqrt( TMath::Power(Run3TotalStat,2.)+TMath::Power(Run3TotalDet,2.)+TMath::Power(Run3TotalG4,2.)+TMath::Power(Run3TotalGenie,2.)+TMath::Power(Run3TotalFlux,2.)+TMath::Power(Run3TotalDirt,2.));

	cout << endl;
	cout << "Stat & " << Run3TotalStat << "\\\\" << endl;
	cout << "Detector & " << Run3TotalDet << "\\\\" << endl;
	cout << "G4 & " << Run3TotalG4 << "\\\\" << endl;
	cout << "XSec & " << Run3TotalGenie << "\\\\" << endl;
	cout << "Flux & " << Run3TotalFlux << "\\\\" << endl;
	cout << "Dirt & " << Run3TotalDirt << "\\\\" << endl;

	cout << "\\hline \\hline" << endl; 
	cout << "Total & " << Run3Total << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << endl << endl;

	// --------------------------------------------------------------------------------------------------------------------------------------------

} // End of the program 
