#include <iostream>
#include <vector>

void script_PRD_DeltaAlphaT_InCosThetaMu_Gene() {

	gROOT->ProcessLine(".L PRD_DeltaAlphaTInCosThetaMuSlices_Gene.cxx++");
	gROOT->ProcessLine("PRD_DeltaAlphaTInCosThetaMuSlices_Gene()");

	gROOT->ProcessLine(".q");		

}
