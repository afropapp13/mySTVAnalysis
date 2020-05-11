#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "TString.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace std;

namespace Constants {

	// Argon 

	static const double A = 40.;
	static const double Z = 18.;

	static const int FontStyle = 132;

	// UBCodeVersion

	static TString UBCodeVersion = "v08_00_00_33";

	// Labels

	static const TString PlotXAxis[] = {
	  "#frac{d#sigma}{d#deltaP_{T}} [10^{-38} cm^{2}]"
	 ,"#frac{d#sigma}{d#delta#alpha_{T}} []10^{-38} cm^{2}]"
	 ,"#frac{d#sigma}{d#delta#phi_{T}} [10^{-38} cm^{2}]"
	 ,"#frac{d#sigma}{dP_{#mu}} [10^{-38} cm^{2}]"
	 ,"#frac{d#sigma}{dcos#theta_{#mu}} [10^{-38} cm^{2}]"
	 ,"#frac{d#sigma}{d#phi_{#mu}} [10^{-38} cm^{2}]"
	 ,"#frac{d#sigma}{dP_{p}} [10^{-38} cm^{2}]"
	 ,"#frac{d#sigma}{dcos#theta_{p}} [10^{-38} cm^{2}]"
	 ,"#frac{d#sigma}{d#phi_{p}} [10^{-38} cm^{2}]"
	 ,"#frac{d#sigma}{dE^{Cal}} [10^{-38} cm^{2}]"
	 ,"#frac{d#sigma}{dE^{QE}} [10^{-38} cm^{2}]"
	 ,"#frac{d#sigma}{dQ^{2}} [10^{-38} cm^{2}]"
	};

	// Genie Constants

	static const double FluxIntegratedXSection = 26.5736; // e-38 cm^2
	static const int NGenieEvents = 1E5;

	// Global Constants

	static const double Units = 1E38; // so that the extracted cross-section is 10^{-38} cm^{2}
	static const double Flux = 3.601E10; //cm^(-2)
	static const double NTargets = 1.1782E30;

	// --------------------------------------------------------------------------------------------------------------------------------------------------------

	// POT Normalization

	// v33 Run 1 

	static const double tor860_wcut_Run1 = 4.131e+19;
	static const double E1DCNT_wcut_Run1 = 9147384.0;
	static const double EXT_Run1 = 34150796.0;

	// --------------------------------------------------------------------------------------------------------------------------------------------------------

	// Binning & XLabels for Selection Cuts & Kinematic Variables

	static const int NBinsECal = 8; static const double ArrayNBinsECal[NBinsECal+1] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8}; 
	static TString LabelXAxisECal = ";E^{Cal} [GeV]"; static TString LabelXAxisTrueECal = ";True E^{Cal} [GeV]";

	static const int NBinsEQE = 8; static const double ArrayNBinsEQE[NBinsEQE+1] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8}; 
	static TString LabelXAxisEQE = ";E^{QE} [GeV]"; static TString LabelXAxisTrueEQE = ";True E^{QE} [GeV]";

	static const int NBinsQ2 = 8; static const double ArrayNBinsQ2[NBinsQ2+1] = { 0.0,0.1,0.2,0.3,0.4,0.5,0.65,0.8,1.}; 
	static TString LabelXAxisQ2 = ";Q^{2} [GeV^{2}/c^{2}]"; static TString LabelXAxisTrueQ2 = ";True Q^{2} [GeV^{2}/c^{2}]";

	static const int NBinsMuonMomentum = 8; static const double ArrayNBinsMuonMomentum[NBinsMuonMomentum+1] = { 0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7}; 
	static TString LabelXAxisMuonMomentum = ";P_{#mu} [GeV/c]"; static TString LabelXAxisTrueMuonMomentum = ";True P_{#mu} [GeV/c]";

	static const int NBinsProtonMomentum = 10; static const double ArrayNBinsProtonMomentum[NBinsProtonMomentum+1] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2}; 
	static TString LabelXAxisProtonMomentum = ";P_{p} [GeV/c]"; static TString LabelXAxisTrueProtonMomentum = ";True P_{p} [GeV/c]";

	static const int NBinsMuonPhi = 7; static const double ArrayNBinsMuonPhi[NBinsMuonPhi+1] = { -180.,-128.6,-77.1,-25.7,25.7,77.1,128.6,180.};
	static TString LabelXAxisMuonPhi = ";#phi_{#mu} [deg]"; static TString LabelXAxisTrueMuonPhi = ";True #phi_{#mu} [deg]";

	static const int NBinsProtonPhi = 7; static const double ArrayNBinsProtonPhi[NBinsProtonPhi+1] = { -180.,-128.6,-77.1,-25.7,25.7,77.1,128.6,180.};
	static TString LabelXAxisProtonPhi = ";#phi_{p} [deg]"; static TString LabelXAxisTrueProtonPhi = ";True #phi_{p} [deg]";

	static const int NBinsMuonCosTheta = 9; static const double ArrayNBinsMuonCosTheta[NBinsMuonCosTheta+1] = { -1.,-0.8,-0.6,-0.4,-0.2,-0.0,0.2,0.5,0.8,1.}; 
	static TString LabelXAxisMuonCosTheta = ";cos#theta_{#mu}"; static TString LabelXAxisTrueMuonCosTheta = ";True cos#theta_{#mu}";

	static const int NBinsProtonCosTheta = 9; static const double ArrayNBinsProtonCosTheta[NBinsProtonCosTheta+1] = { -1.,-0.8, -0.6,-0.3,0.,0.2,0.4,0.6,0.8,1. }; 
	static TString LabelXAxisProtonCosTheta = ";cos#theta_{p}"; static TString LabelXAxisTrueProtonCosTheta = ";True cos#theta_{p}";

	static const int NBinsChi2 = 25; TString RecoLabelXAxisChi2 = ";#chi^{2}_{p}";
	static const double MinChi2 = 0., MaxChi2 = 500.;

	static const int NBinsThreePlaneChi2 = 25; TString RecoLabelXAxisThreePlaneChi2 = ";3-plane #chi^{2}_{p}";
	static const double MinThreePlaneChi2 = 0., MaxThreePlaneChi2 = 500.;

	static const int NBinsThreePlaneChi2LogLikelihood = 40; TString RecoLabelXAxisThreePlaneChi2LogLikelihood = ";3-Plane LogLikelihood";
	static const double MinThreePlaneChi2LogLikelihood = -5., MaxThreePlaneChi2LogLikelihood = 5.;
	TString RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood = ";3-Plane LogLikelihood #mu-Candidate";
	TString RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood = ";3-Plane LogLikelihood p-Candidate";

	static const int NBinsDeltaTheta = 18; TString RecoLabelXAxisDeltaTheta = ";#Delta#theta";
	static const double MinDeltaTheta = 0., MaxDeltaTheta = 180.;

	static const int NBinsDeltaPhi = 15; TString RecoLabelXAxisDeltaPhi = ";#Delta#phi";
	static const double MinDeltaPhi = 0., MaxDeltaPhi = 360.;

	static const int NBinsdYZ = 20; TString RecoLabelXAxisdYZ = ";d_{YZ} (cm)";
	static const double MindYZ = 0., MaxdYZ = 500.;

	static const int NBinsNPE = 20; TString RecoLabelXAxisNPE = ";# PE";
	static const double MinNPE = 0., MaxNPE = 1000.;

	static const int NBinsDistance = 15; TString RecoLabelXAxisDistance = ";#mu - p distance (cm)";
	static const double MinDistance = 0., MaxDistance = 5.;

	static const int NBinsNuScore = 25; TString RecoLabelXAxisNuScore = ";#nu score";
	static const double MinNuScore = 0.05, MaxNuScore = 1.;

	static const int NBinsFlashScore = 25; TString RecoLabelXAxisFlashScore = ";Flash score";
	static const double MinFlashScore = 0., MaxFlashScore = 50.;

	// --------------------------------------------------------------------------------------------------------------------------------------------------------

	// Binning & XLabels for STV

	static const int NBinsDeltaPT = 9; static const double ArrayNBinsDeltaPT[NBinsDeltaPT+1] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9}; 
	static TString LabelXAxisDeltaPT = ";#deltap_{T} [GeV/c]"; static TString LabelXAxisTrueDeltaPT = ";True #deltap_{T} [GeV/c]";

	static const int NBinsDeltaAlphaT = 7; static const double ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT+1] = { 0.,25.,50.,75.,100.,125.,150.,180. }; 
	static TString LabelXAxisDeltaAlphaT = ";#delta#alpha_{T} [deg]"; static TString LabelXAxisTrueDeltaAlphaT = ";True #delta#alpha_{T} [deg]";

	static const int NBinsDeltaPhiT = 10; static const double ArrayNBinsDeltaPhiT[NBinsDeltaPhiT+1] = { 0.,12.5,25.,37.5,50.,62.5,75.,90., 105., 120., 135. }; 
	static TString LabelXAxisDeltaPhiT = ";#delta#phi_{T} [deg]"; static TString LabelXAxisTrueDeltaPhiT = ";True #delta#phi_{T} [deg]";

	// --------------------------------------------------------------------------------------------------------------------------------------------------------

	// Labels for 2D Plots

	static TString LabelXAxisDeltaPT2D = LabelXAxisTrueDeltaPT+";Reco #deltap_{T} [GeV/c]";
	static TString LabelXAxisDeltaAlphaT2D = LabelXAxisTrueDeltaAlphaT+";Reco #delta#alpha_{T} [deg]";
	static TString LabelXAxisDeltaPhiT2D = LabelXAxisTrueDeltaPhiT+";Reco #delta#phi_{T} [deg]";

	static TString LabelXAxisMuonMomentum2D = LabelXAxisTrueMuonMomentum+";Reco P_{#mu} [GeV/c]";
	static TString LabelXAxisMuonPhi2D = LabelXAxisTrueMuonPhi+";Reco #phi_{#mu} [deg]";
	static TString LabelXAxisMuonCosTheta2D = LabelXAxisTrueMuonCosTheta+";Reco cos#theta_{#mu}";

	static TString LabelXAxisProtonMomentum2D = LabelXAxisTrueProtonMomentum+";Reco P_{p} [GeV/c]";
	static TString LabelXAxisProtonPhi2D = LabelXAxisTrueProtonPhi+";Reco #phi_{p} [deg]";
	static TString LabelXAxisProtonCosTheta2D = LabelXAxisTrueProtonCosTheta+";Reco cos#theta_{p}";

	static TString LabelXAxisECal2D = ";True E^{Cal} [GeV];Reco E^{Cal} [GeV]";
	static TString LabelXAxisEQE2D = ";True E^{QE} [GeV];Reco E^{QE} [GeV]";
	static TString LabelXAxisQ22D = ";True Q^{2} [GeV^{2}/c^{2}];Reco Q^{2} [GeV^{2}/c^{2}]";

	// --------------------------------------------------------------------------------------------------------------------------------------------------------

	static int OverlayColor = kAzure-4;
	static int BeamOnColor = kBlack;
	static int GenieColor = 610;

	// --------------------------------------------------------------------------------------------------------------------------------------------------------

	// Constants, Cuts & Thresholds

	static const double NPE = 10;

	static const double MuonChi2Cut = 100.;
	static const double ProtonChi2Cut = 80.;

	static const double MuonThreePlaneChi2LogLikelihoodCut = -1.;
	static const double ProtonThreePlaneChi2LogLikelihoodCut = -0.5;

	static const int MuonPdg = 13, ProtonPdg = 2212, AbsChargedPionPdg = 211, NeutralPionPdg = 111;

	static const double MuonMass = 106, ProtonMass = 938; // MeV
	static const double MuonMass_GeV = 0.106, ProtonMass_GeV = 0.938; // GeV

	static const double CosmicPID = -999.;
	static const int CosmicPdg = -99;

	static const double BE = 0.04; // GeV

	static const double DeltaRCut = 3; // cm
	static const double ThresholdPurity = 0.1; // 10%
	static const double UBSpaceReso = 0.3; // cm // 3mm spacing between the wires

	static const double DeltaThetaCut = 180.; // deg
	static const double DeltaPhiCentralValue = 180.; // deg
	static const double DeltaPhiOpeningAngle = 35.; // deg
	static const double MaxTransMissMomentum = 0.35; // deg

	static const double MaxMuPDistance = 2.; // cm

	static const double MaxPFParticleTrackDistance = 20.; // cm

	static const double YZBeamFlashVertexMaxDistance = 500.; // cm
	static const double BeamFlashPEThreshold = 100; // PEs

	static const double ChargedPionMomentumThres = 0.07;
	static const double NeutralPionMomentumThres = 0.0;

	static const double PurityThreshold = 0.1;

	static const double MinimumNuScore = 0.6;

	// --------------------------------------------------------------------------------------------------------------------------------------------------------

}
#endif
