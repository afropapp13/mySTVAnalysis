#include <iostream>
#include <vector>

void script_PRD_DeltaPT_InCosThetaMu_Inte() {

	gROOT->ProcessLine(".L PRD_DeltaPTInCosThetaMu_InteBreakdown.cxx++");
	gROOT->ProcessLine("PRD_DeltaPTInCosThetaMu_InteBreakdown()");	

	gROOT->ProcessLine(".q");	

}
