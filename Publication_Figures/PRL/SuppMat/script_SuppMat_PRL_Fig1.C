#include <iostream>
#include <vector>

void script_SuppMat_PRL_Fig1() {

	gROOT->ProcessLine(".L PRL_SuppMat_Fig1_DeltaPTInDeltaAlphaTSlices.cxx++");
	gROOT->ProcessLine("PRL_SuppMat_Fig1_DeltaPTInDeltaAlphaTSlices()");	

}
