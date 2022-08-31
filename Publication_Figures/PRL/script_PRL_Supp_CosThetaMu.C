#include <iostream>
#include <vector>

void script_PRL_Supp_CosThetaMu() {

	gROOT->ProcessLine(".L PRL_SuppMat_CosThetaMu.cxx++");
	gROOT->ProcessLine("PRL_SuppMat_CosThetaMu()");	

}
