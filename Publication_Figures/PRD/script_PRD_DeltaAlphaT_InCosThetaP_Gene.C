#include <iostream>
#include <vector>

void script_PRD_DeltaAlphaT_InCosThetaP_Gene() {

	gROOT->ProcessLine(".L PRD_DeltaAlphaTInCosThetaPSlices_Gene.cxx++");
	gROOT->ProcessLine("PRD_DeltaAlphaTInCosThetaPSlices_Gene()");

	gROOT->ProcessLine(".q");		

}
