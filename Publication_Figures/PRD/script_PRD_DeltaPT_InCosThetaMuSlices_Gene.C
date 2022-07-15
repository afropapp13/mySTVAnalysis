#include <iostream>
#include <vector>

void script_PRD_DeltaPT_InCosThetaMuSlices_Gene() {

	gROOT->ProcessLine(".L PRD_DeltaPTInCosThetaMuSlices_Gene.cxx++");
	gROOT->ProcessLine("PRD_DeltaPTInCosThetaMuSlices_Gene()");	

	gROOT->ProcessLine(".q");	

}
