#include <iostream>
#include <vector>

void script_SuppMat_PRL_Fig3() {

	gROOT->ProcessLine(".L PRL_SuppMat_Fig3_DeltaPtxInDeltaPtySlices.cxx++");
	gROOT->ProcessLine("PRL_SuppMat_Fig3_DeltaPtxInDeltaPtySlices()");	

}
