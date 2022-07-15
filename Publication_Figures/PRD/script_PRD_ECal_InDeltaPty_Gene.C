#include <iostream>
#include <vector>

void script_PRD_ECal_InDeltaPty_Gene() {

	gROOT->ProcessLine(".L PRD_ECalInDeltaPtySlices_Gene.cxx++");
	gROOT->ProcessLine("PRD_ECalInDeltaPtySlices_Gene()");

	gROOT->ProcessLine(".q");		

}
