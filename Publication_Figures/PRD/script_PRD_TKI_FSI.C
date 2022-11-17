#include <iostream>
#include <vector>

void script_PRD_TKI_FSI() {

	gROOT->ProcessLine(".L PRD_TKI_FSI.cxx++");
	gROOT->ProcessLine("PRD_TKI_FSI()");	

	gROOT->ProcessLine(".q");	

}
