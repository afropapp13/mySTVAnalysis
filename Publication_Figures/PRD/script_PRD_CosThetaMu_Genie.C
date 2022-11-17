#include <iostream>
#include <vector>

void script_PRD_CosThetaMu_Genie() {

	gROOT->ProcessLine(".L PRD_CosThetaMu_Genie.cxx++");
	gROOT->ProcessLine("PRD_CosThetaMu_Genie()");

	gROOT->ProcessLine(".q");		

}
