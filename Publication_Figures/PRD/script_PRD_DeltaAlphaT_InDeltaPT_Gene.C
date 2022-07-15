#include <iostream>
#include <vector>

void script_PRD_DeltaAlphaT_InDeltaPT_Gene() {

	gROOT->ProcessLine(".L PRD_DeltaAlphaTInDeltaPTSlices_Gene.cxx++");
	gROOT->ProcessLine("PRD_DeltaAlphaTInDeltaPTSlices_Gene()");

	gROOT->ProcessLine(".q");		

}
