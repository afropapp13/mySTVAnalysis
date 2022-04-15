#include <iostream>
#include <vector>

void script_PRL_Fig2() {

	gROOT->ProcessLine(".L PRL_Fig2_DeltaAlphaTInDeltaPTSlices.cxx++");
	gROOT->ProcessLine("PRL_Fig2_DeltaAlphaTInDeltaPTSlices()");	

}
