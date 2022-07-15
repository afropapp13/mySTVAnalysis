#include <iostream>
#include <vector>

void script_PRD_DeltaPT_InCosThetaPSlices_Genie() {

	gROOT->ProcessLine(".L PRD_DeltaPTInCosThetaPSlices_Genie.cxx++");
	gROOT->ProcessLine("PRD_DeltaPTInCosThetaPSlices_Genie()");	

	gROOT->ProcessLine(".q");	

}
