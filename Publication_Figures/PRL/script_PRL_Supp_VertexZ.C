#include <iostream>
#include <vector>

void script_PRL_Supp_VertexZ() {

	gROOT->ProcessLine(".L PRL_Supp_VertexZ.cxx++");
	gROOT->ProcessLine("PRL_Supp_VertexZ()");	

}
