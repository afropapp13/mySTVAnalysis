#include <iostream>
#include <vector>

void script_PRD_DeltaPtx_InDeltaPty_Genie() {

	gROOT->ProcessLine(".L PRD_DeltaPtxInDeltaPtySlices_Genie.cxx++");
	gROOT->ProcessLine("PRD_DeltaPtxInDeltaPtySlices_Genie()");

	gROOT->ProcessLine(".q");		

}
