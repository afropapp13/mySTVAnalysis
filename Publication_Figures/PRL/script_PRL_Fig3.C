#include <iostream>
#include <vector>

void script_PRL_Fig3() {

	gROOT->ProcessLine(".L PRL_Fig3_DeltaPtxInDeltaPtySlices.cxx++");
	gROOT->ProcessLine("PRL_Fig3_DeltaPtxInDeltaPtySlices()");	

}
