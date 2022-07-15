#include <iostream>
#include <vector>

void script_PRD_DeltaPhiT_InDeltaPTSlices_Inte() {

	gROOT->ProcessLine(".L PRD_DeltaPhiTInDeltaPT_InteBreakdown.cxx++");
	gROOT->ProcessLine("PRD_DeltaPhiTInDeltaPT_InteBreakdown()");	

	gROOT->ProcessLine(".q");	

}
