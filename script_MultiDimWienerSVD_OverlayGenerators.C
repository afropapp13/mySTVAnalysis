{

	gROOT->ProcessLine(".L ../myClasses/Util.C");
	gROOT->ProcessLine(".L ../myClasses/Tools.cxx");	
	gROOT->ProcessLine(".L MultiDimWienerSVD_OverlayGenerators.cpp");
	//GENIE versions
//	gROOT->ProcessLine("MultiDimWienerSVD_OverlayGenerators()");
	//AltGen
	gROOT->ProcessLine("MultiDimWienerSVD_OverlayGenerators(false,true)");
	//GENIE FSI tweaks
//	gROOT->ProcessLine("MultiDimWienerSVD_OverlayGenerators(false,false,true)");
	//GENIE Flag tweaks
//	gROOT->ProcessLine("MultiDimWienerSVD_OverlayGenerators(false,false,false,true)");
	//GENIE Closure Test
//	gROOT->ProcessLine("MultiDimWienerSVD_OverlayGenerators(false,false,false,false,true)");			
	//GENIE Nuclear model tweaks
//	gROOT->ProcessLine("MultiDimWienerSVD_OverlayGenerators(false,false,false,false,false,true)");
	//NuWro Comparisons out-of-the-box vs larsoft
//	gROOT->ProcessLine("MultiDimWienerSVD_OverlayGenerators(false,false,false,false,false,false,true)");
	//Just data vs Nominal MC
//	gROOT->ProcessLine("MultiDimWienerSVD_OverlayGenerators(false,false,false,false,false,false,false,true)");		

}
