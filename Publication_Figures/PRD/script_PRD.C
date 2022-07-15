#include <iostream>
#include <vector>

void script_PRD() {

	gROOT->ProcessLine(".x script_PRD_DeltaPT_InDeltaAlphaT_Gene.C");
	gROOT->ProcessLine(".x script_PRD_DeltaPT_InDeltaAlphaT_Genie.C");	

	gROOT->ProcessLine(".x script_PRD_DeltaPT_InCosThetaMu_Gene.C");
	gROOT->ProcessLine(".x script_PRD_DeltaPT_InCosThetaMu_Genie.C");	
	gROOT->ProcessLine(".x script_PRD_DeltaPT_InCosThetaMuSlices_Inte.C");

	gROOT->ProcessLine(".x script_PRD_DeltaPT_InCosThetaP_Gene.C");
	gROOT->ProcessLine(".x script_PRD_DeltaPT_InCosThetaP_Genie.C");	

	gROOT->ProcessLine(".x script_PRD_DeltaAlphaT_InDeltaPT_Gene.C");
	gROOT->ProcessLine(".x script_PRD_DeltaAlphaT_InDeltaPT_Genie.C");	
//	gROOT->ProcessLine(".x script_PRD_DeltaAlphaT_InDeltaPT_FSI.C");



}
