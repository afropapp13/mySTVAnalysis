#include <iostream>
#include <vector>

void script_PRD_DeltaPhiT_InDeltaPT_Genie() {

	gROOT->ProcessLine(".L PRD_DeltaPhiTInDeltaPTSlices_Genie.cxx++");
	gROOT->ProcessLine("PRD_DeltaPhiTInDeltaPTSlices_Genie()");

	gROOT->ProcessLine(".q");		

}
