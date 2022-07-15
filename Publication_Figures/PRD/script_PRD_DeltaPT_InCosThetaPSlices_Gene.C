#include <iostream>
#include <vector>

void script_PRD_DeltaPT_InCosThetaPSlices_Gene() {

	gROOT->ProcessLine(".L PRD_DeltaPTInCosThetaPSlices_Gene.cxx++");
	gROOT->ProcessLine("PRD_DeltaPTInCosThetaPSlices_Gene()");	

	gROOT->ProcessLine(".q");	

}
