#include <iostream>
#include <vector>

void script_PRD_DeltaPhiT_InDeltaPT_Gene() {

	gROOT->ProcessLine(".L PRD_DeltaPhiTInDeltaPTSlices_Gene.cxx++");
	gROOT->ProcessLine("PRD_DeltaPhiTInDeltaPTSlices_Gene()");

	gROOT->ProcessLine(".q");		

}
