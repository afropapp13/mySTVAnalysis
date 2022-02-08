{

	gROOT->ProcessLine(".L ../myClasses/Util.C");
	gROOT->ProcessLine(".L ../myClasses/Tools.cxx");	
	gROOT->ProcessLine(".L TwoDimWienerSVD_OverlayGenerators.cpp");
	//GENIE versions
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators()");
	//AltGen
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,true)");
	//GENIE FSI tweaks
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,false,true)");
	//GENIE Flag tweaks
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,false,false,true)");
	//GENIE Closure Test
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,false,false,false,true)");			
	//GENIE Nuclear model tweaks
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,false,false,false,false,true)");
	//NuWro Comparisons out-of-the-box vs larsoft
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,false,false,false,false,false,true)");
	//Just data vs Nominal MC
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,false,false,false,false,false,false,true)");	
	//All
	//Just data vs Nominal MC
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,false,false,false,false,false,false,false,false,true)");		

}
