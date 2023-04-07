#include <iostream>
#include <vector>

void script_PRL_Noah_NuclModels() {

	gROOT->ProcessLine(".L PRL_Noah_NuclModel_DeltaAlphaTInDeltaPTSlices.cxx++");
	gROOT->ProcessLine("PRL_Noah_NuclModel_DeltaAlphaTInDeltaPTSlices()");	

}
