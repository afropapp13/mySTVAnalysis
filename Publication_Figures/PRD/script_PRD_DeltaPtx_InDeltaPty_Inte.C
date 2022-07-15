#include <iostream>
#include <vector>

void script_PRD_DeltaPtx_InDeltaPty_Inte() {

	gROOT->ProcessLine(".L PRD_DeltaPtxInDeltaPty_InteBreakdown.cxx++");
	gROOT->ProcessLine("PRD_DeltaPtxInDeltaPty_InteBreakdown()");

	gROOT->ProcessLine(".q");		

}
