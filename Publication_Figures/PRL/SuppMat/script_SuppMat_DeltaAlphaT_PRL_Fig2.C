#include <iostream>
#include <vector>

void script_SuppMat_DeltaAlphaT_PRL_Fig2() {

	gROOT->ProcessLine(".L PRL_SuppMat_GiBUU_Fig2_DeltaAlphaTInDeltaPTSlices.cxx++");
	gROOT->ProcessLine("PRL_SuppMat_GiBUU_Fig2_DeltaAlphaTInDeltaPTSlices()");	

}
