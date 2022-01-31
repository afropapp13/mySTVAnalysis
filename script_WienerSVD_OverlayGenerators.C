{

	gROOT->ProcessLine(".L ../myClasses/Util.C");
	gROOT->ProcessLine(".L WienerSVD_OverlayGenerators.cpp");
	//GENIE versions
	gROOT->ProcessLine("WienerSVD_OverlayGenerators()");
	//AltGen
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,true)");
	//GENIE FSI tweaks
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,true)");
	//GENIE Flag tweaks
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,true)");
	//GENIE Closure Test
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,false,true)");			
	//GENIE Nuclear model tweaks
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,false,false,true)");
	//NuWro Comparisons out-of-the-box vs larsoft
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,false,false,false,true)");
	//Just data vs Nominal MC
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,false,false,false,false,true)");		

}
