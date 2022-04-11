#include <iostream>
#include <vector>

void script_OverlayRecoTruth() {

	gROOT->ProcessLine(".L ../../myClasses/Util.C++");
	gROOT->ProcessLine(".L OverlayRecoTrueGen.cpp++");
	gROOT->ProcessLine("OverlayRecoTrueGen()");	

}
