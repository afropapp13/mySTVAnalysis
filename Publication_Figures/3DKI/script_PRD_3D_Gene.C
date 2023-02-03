#include <iostream>
#include <vector>

void script_PRD_3D_Gene() {

	gROOT->ProcessLine(".L PRD_3D_Gene.cxx++");
	gROOT->ProcessLine("PRD_3D_Gene()");

//	gROOT->ProcessLine(".q");		

}
