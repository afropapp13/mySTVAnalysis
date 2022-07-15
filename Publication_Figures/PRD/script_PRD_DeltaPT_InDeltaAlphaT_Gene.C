#include <iostream>
#include <vector>

void script_PRD_DeltaPT_InDeltaAlphaT_Gene() {

	gROOT->ProcessLine(".L PRD_DeltaPTInDeltaAlphaTSlices_Gene.cxx++");
	gROOT->ProcessLine("PRD_DeltaPTInDeltaAlphaTSlices_Gene()");

	gROOT->ProcessLine(".q");		

}
