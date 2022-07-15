#include <iostream>
#include <vector>

void script_PRD_DeltaPtx_InDeltaPty_Gene() {

	gROOT->ProcessLine(".L PRD_DeltaPtxInDeltaPtySlices_Gene.cxx++");
	gROOT->ProcessLine("PRD_DeltaPtxInDeltaPtySlices_Gene()");

	gROOT->ProcessLine(".q");		

}
