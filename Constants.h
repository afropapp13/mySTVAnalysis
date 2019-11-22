#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "TString.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace std;

namespace Constants {

	static const int FontStyle = 132;

	// Run Period

	static const TString WhichRun = "Run1";
//	static const TString WhichRun = "Run3";

	// Samples

	static TString UBCodeVersion = "v08_00_00_22";

	// Labels

	static const TString PlotXAxis[] = {
	  "#frac{d#sigma}{d#deltaP_{T}} (10^{-38} cm^{2})" 
	 ,"#frac{d#sigma}{d#delta#alpha_{T}} (10^{-38} cm^{2})"
	 ,"#frac{d#sigma}{d#delta#phi_{T}} (10^{-38} cm^{2})"
	 ,"#frac{d#sigma}{dP_{#mu}} (10^{-38} cm^{2})"
	 ,"#frac{d#sigma}{dcos(#theta_{#mu})} (10^{-38} cm^{2})"
	 ,"#frac{d#sigma}{d#phi_{#mu}} (10^{-38} cm^{2})"
	 ,"#frac{d#sigma}{dP_{p}} (10^{-38} cm^{2})"
	 ,"#frac{d#sigma}{dcos(#theta_{p})} (10^{-38} cm^{2})"
	 ,"#frac{d#sigma}{d#phi_{p}} (10^{-38} cm^{2})"
	 ,"#frac{d#sigma}{dE^{Cal}} (10^{-38} cm^{2})"
	 ,"#frac{d#sigma}{dQ^{2}} (10^{-38} cm^{2})"
	};

	// Genie Constants

	static const double FluxIntegratedXSection = 26.5736; // e-38 cm^2
	static const int NGenieEvents = 1E6;

	// Global Constants

	static const double Units = 1E38; // so that the extracted cross-section is 10^{-38} cm^{2}
	static const double Flux = 3.601E10; //cm^(-2)
	static const double NTargets = 1.1782E30;

	// --------------------------------------------------------------------------------------------------------------------------------------------------------

	// POT Normalization

	// v20 & v22 Run 1 

	static const double tor860_wcut = 4.405e+19;
	static const double E1DCNT_wcut = 9771236.0;

	static const double EXT = 14675888.0;
	static const double DirtPOT = 3.19485e+20;
	static const double OverlayPOT = 8.92987e+19; // CV

	// Systematics

	static const double OverlayPOT_SCE = 1.11365e+20; // SCE
	static const double OverlayPOT_DLdown = 1.03462e+20; // DLdown

	// --------------------------------------------------------------------------------------------------------------------------------------------------------

	// Binning & XLabels for Selection Cuts & Kinematic Variables

	static const int NBinsECal = 7; static const double ArrayNBinsECal[NBinsECal+1] = { 0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6}; 
	static TString LabelXAxisECal = ";E^{Cal} (GeV)";

	static const int NBinsQ2 = 7; static const double ArrayNBinsQ2[NBinsQ2+1] = { 0.0,0.1,0.2,0.3,0.4,0.5,0.65,0.8}; 
	static TString LabelXAxisQ2 = ";Q^{2} (GeV^{2}/c^{2})";

	static const int NBinsMuonMomentum = 7; static const double ArrayNBinsMuonMomentum[NBinsMuonMomentum+1] = { 0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5 }; 
	static TString LabelXAxisMuonMomentum = ";P_{#mu} (GeV/c)"; static TString LabelXAxisTrueMuonMomentum = ";True P_{#mu} (GeV/c)";

	static const int NBinsProtonMomentum = 7; static const double ArrayNBinsProtonMomentum[NBinsProtonMomentum+1] = { 0.3,0.4,0.5,0.6,0.7,0.8,0.9,1. }; 
	static TString LabelXAxisProtonMomentum = ";P_{p} (GeV/c)"; static TString LabelXAxisTrueProtonMomentum = ";True P_{p} (GeV/c)";

	static const int NBinsMuonPhi = 7; static const double ArrayNBinsMuonPhi[NBinsMuonPhi+1] = { -180.,-128.6,-77.1,-25.7,25.7,77.1,128.6,180. }; 
	static TString LabelXAxisMuonPhi = ";#phi_{#mu} (deg)"; static TString LabelXAxisTrueMuonPhi = ";True #phi_{#mu} (deg)";

	static const int NBinsProtonPhi = 7; static const double ArrayNBinsProtonPhi[NBinsProtonPhi+1] = { -180.,-128.6,-77.1,-25.7,25.7,77.1,128.6,180. }; 
	static TString LabelXAxisProtonPhi = ";#phi_{p} (deg)"; static TString LabelXAxisTrueProtonPhi = ";True #phi_{p} (deg)";

	static const int NBinsMuonCosTheta = 7; static const double ArrayNBinsMuonCosTheta[NBinsMuonCosTheta+1] = { -0.6,-0.4,-0.2,-0.0,0.2,0.5,0.8,1. }; 
	static TString LabelXAxisMuonCosTheta = ";cos(#theta_{#mu})"; static TString LabelXAxisTrueMuonCosTheta = ";True cos(#theta_{#mu})";

	static const int NBinsProtonCosTheta = 7; static const double ArrayNBinsProtonCosTheta[NBinsProtonCosTheta+1] = { -0.6,-0.3,0.,0.2,0.4,0.6,0.8,1. }; 
	static TString LabelXAxisProtonCosTheta = ";cos(#theta_{p})"; static TString LabelXAxisTrueProtonCosTheta = ";True cos(#theta_{p})";

	static const int NBinsChi2 = 25; TString RecoLabelXAxisChi2 = ";#chi^{2}_{p}";
	static const double MinChi2 = 0., MaxChi2 = 500.;

	static const int NBinsThreePlaneChi2 = 25; TString RecoLabelXAxisThreePlaneChi2 = ";3-plane #chi^{2}_{p}";
	static const double MinThreePlaneChi2 = 0., MaxThreePlaneChi2 = 500.;

	static const int NBinsThreePlaneChi2LogLikelihood = 20; TString RecoLabelXAxisThreePlaneChi2LogLikelihood = ";3-Plane LogLikelihood";
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
	static const double MinNuScore = 0., MaxNuScore = 1.;

	static const int NBinsFlashScore = 25; TString RecoLabelXAxisFlashScore = ";Flash score";
	static const double MinFlashScore = 0., MaxFlashScore = 50.;

	// --------------------------------------------------------------------------------------------------------------------------------------------------------

	// Binning & XLabels for STV

	static const int NBinsDeltaPT = 9; static const double ArrayNBinsDeltaPT[NBinsDeltaPT+1] = { 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 }; 
	static TString LabelXAxisDeltaPT = ";#deltaP_{T} (GeV/c)";

	static const int NBinsDeltaAlphaT = 7; static const double ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT+1] = { 0.,25.,50.,75.,100.,125.,150.,180. }; 
	static TString LabelXAxisDeltaAlphaT = ";#delta#alpha_{T} (deg)";

	static const int NBinsDeltaPhiT = 10; static const double ArrayNBinsDeltaPhiT[NBinsDeltaPhiT+1] = { 0.,12.5,25.,37.5,50.,62.5,75.,90., 105., 120., 135. }; 
	static TString LabelXAxisDeltaPhiT = ";#delta#phi_{T} (deg)";

	static TString LabelXAxisDeltaPT2D = ";True #deltaP_{T} (GeV/c);Reco #deltaP_{T} (GeV/c)";
	static TString LabelXAxisDeltaAlphaT2D = ";True #delta#alpha_{T} (deg);Reco #delta#alpha_{T} (deg)";
	static TString LabelXAxisDeltaPhiT2D = ";True #delta#phi_{T} (deg);Reco #delta#phi_{T} (deg)";

	static TString LabelXAxisMuonMomentum2D = ";True P_{#mu} (GeV/c);Reco P_{#mu} (GeV/c)";
	static TString LabelXAxisMuonPhi2D = ";True #phi_{#mu} (deg);Reco #phi_{#mu} (deg)";
	static TString LabelXAxisMuonCosTheta2D = ";True cos(#theta_{#mu});Reco cos(#theta_{#mu})";

	static TString LabelXAxisProtonMomentum2D = ";True P_{p} (GeV/c);Reco P_{p} (GeV/c)";
	static TString LabelXAxisProtonPhi2D = ";True #phi_{p} (deg);Reco #phi_{p} (deg)";
	static TString LabelXAxisProtonCosTheta2D = ";True cos(#theta_{p});Reco cos(#theta_{p})";

	static TString LabelXAxisECal2D = ";True E^{Cal} (GeV);Reco E^{Cal} (GeV)";
	static TString LabelXAxisQ22D = ";True Q^{2} (GeV^{2}/c^{2});Reco Q^{2} (GeV^{2}/c^{2})";

	// --------------------------------------------------------------------------------------------------------------------------------------------------------

	// Constants, Cuts & Thresholds

	static const double NPE = 100;

	static const double MuonChi2Cut = 100.;
	static const double ProtonChi2Cut = 80.;

	static const double MuonThreePlaneChi2LogLikelihoodCut = -1.;
	static const double ProtonThreePlaneChi2LogLikelihoodCut = -2.;

	static const int MuonPdg = 13, ProtonPdg = 2212, AbsChargedPionPdg = 211, NeutralPionPdg = 111;

	static const double MuonMass = 106, ProtonMass = 938; // MeV
	static const double MuonMass_GeV = 0.106, ProtonMass_GeV = 0.938; // GeV

	static const double CosmicPID = -999.;
	static const int CosmicPdg = -99;

	static const double BE = 0.04; // GeV

	static const double DeltaRCut = 3; // cm
	static const double ThresholdPurity = 0.1; // 10%
	static const double UBSpaceReso = 0.3; // cm // 3mm spacing between the wires

	static const double DeltaThetaCentralValue = 90.; // deg
	static const double DeltaThetaOpeningAngle = 70.; // deg
	static const double DeltaPhiCentralValue = 180.; // deg
	static const double DeltaPhiOpeningAngle = 35.; // deg // default 35
	static const double MaxTransMissMomentum = 0.35; // deg

	static const double MaxMuPDistance = 2.; // cm

	static const double MaxPFParticleTrackDistance = 20.; // cm

	static const double YZBeamFlashVertexMaxDistance = 500.; // cm
	static const double BeamFlashPEThreshold = 100; // PEs

	static const double ChargedPionMomentumThres = 0.07;
	static const double NeutralPionMomentumThres = 0.0;

	static const double PurityThreshold = 0.1;

	static const double MinimumNuScore = 0.4;

}
#endif
