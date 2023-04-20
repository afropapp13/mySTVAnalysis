#include <iostream>
#include <vector>

void script_PRD_DeltaPT_InCosThetaMu_Inte_GiBUU() {

	gROOT->ProcessLine(".L PRD_DeltaPTInCosThetaMu_InteBreakdown_GiBUU.cxx++");
	gROOT->ProcessLine("PRD_DeltaPTInCosThetaMu_InteBreakdown_GiBUU()");		

	gROOT->ProcessLine(".q");	

}
