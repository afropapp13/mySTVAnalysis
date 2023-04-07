#include <iostream>
#include <vector>

void script_PRD_DeltaPT_FSI_NoFSI() {

	gROOT->ProcessLine(".L PRD_DeltaPT_FSI_NoFSI.cxx++");
	gROOT->ProcessLine("PRD_DeltaPT_FSI_NoFSI()");	

	gROOT->ProcessLine(".q");	

}
