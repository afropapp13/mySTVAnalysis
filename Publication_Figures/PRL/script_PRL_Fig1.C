#include <iostream>
#include <vector>

void script_PRL_Fig1() {

	gROOT->ProcessLine(".L PRL_Fig1_DeltaPTInDeltaAlphaTSlices.cxx++");
	gROOT->ProcessLine("PRL_Fig1_DeltaPTInDeltaAlphaTSlices()");	

}
