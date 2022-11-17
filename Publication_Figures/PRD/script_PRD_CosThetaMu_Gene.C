#include <iostream>
#include <vector>

void script_PRD_CosThetaMu_Gene() {

	gROOT->ProcessLine(".L PRD_CosThetaMu_Gene.cxx++");
	gROOT->ProcessLine("PRD_CosThetaMu_Gene()");

	gROOT->ProcessLine(".q");		

}
