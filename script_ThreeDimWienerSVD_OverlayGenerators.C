{

	gROOT->ProcessLine(".L ../myClasses/Util.C");
	gROOT->ProcessLine(".L ../myClasses/Tools.cxx");	
	gROOT->ProcessLine(".L ThreeDimWienerSVD_OverlayGenerators.cpp");
	//GENIE versions
//	gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators()");
	//AltGen
	gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(false,true)");
	//GENIE FSI tweaks
//	gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(false,false,true)");
	//GENIE Flag tweaks
//	gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(false,false,false,true)");
	//GENIE Closure Test
//	gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(false,false,false,false,true)");			
	//GENIE Nuclear model tweaks
//	gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(false,false,false,false,false,true)");
	//NuWro Comparisons out-of-the-box vs larsoft
//	gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(false,false,false,false,false,false,true)");
	//Just data vs Nominal MC
//	gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(false,false,false,false,false,false,false,true)");		

}
