#include <iostream>
#include <vector>

void script_PRD_DeltaPT_InCosThetaMuSlices_Genie() {

	gROOT->ProcessLine(".L PRD_DeltaPTInCosThetaMuSlices_Genie.cxx++");
	gROOT->ProcessLine("PRD_DeltaPTInCosThetaMuSlices_Genie()");	

	gROOT->ProcessLine(".q");	

}
