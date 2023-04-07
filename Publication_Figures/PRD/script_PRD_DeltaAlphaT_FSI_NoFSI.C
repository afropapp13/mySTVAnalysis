#include <iostream>
#include <vector>

void script_PRD_DeltaAlphaT_FSI_NoFSI() {

	gROOT->ProcessLine(".L PRD_DeltaAlphaT_FSI_NoFSI.cxx++");
	gROOT->ProcessLine("PRD_DeltaAlphaT_FSI_NoFSI()");	

	gROOT->ProcessLine(".q");	

}
