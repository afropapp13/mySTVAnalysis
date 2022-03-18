{

	//------------------------------//

	std::vector<TString> PlotNames; std::vector<int> NumberFirstDiscrSlices;  std::vector<int> NumberSecondDiscrSlices; 

	gROOT->ProcessLine(".L ../myClasses/Util.C");
	gROOT->ProcessLine(".L ../myClasses/Tools.cxx");	
	gROOT->ProcessLine(".L ThreeDimWienerSVD_OverlayGenerators.cpp");	

	// PlotName, First Discriminator (e.g. DeltaPT) , First Discriminator  (e.g. DeltaAlphaT)

	PlotNames.push_back("SerialECal_DeltaPTDeltaAlphaTPlot"); NumberFirstDiscrSlices.push_back(3); NumberSecondDiscrSlices.push_back(4);
	PlotNames.push_back("SerialECal_DeltaPtxDeltaPtyPlot"); NumberFirstDiscrSlices.push_back(3); NumberSecondDiscrSlices.push_back(3);
	PlotNames.push_back("SerialECal_MuonCosThetaMuonMomentumPlot"); NumberFirstDiscrSlices.push_back(4); NumberSecondDiscrSlices.push_back(3);
	PlotNames.push_back("SerialECal_ProtonCosThetaProtonMomentumPlot"); NumberFirstDiscrSlices.push_back(4); NumberSecondDiscrSlices.push_back(3);		

	for (int i = 0; i < (int)(PlotNames.size()); i++) {


		for (int j = 0; j < NumberFirstDiscrSlices.at(i); j++) {		

			//GENIE versions
			gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(\""+PlotNames[i]+"\","+TString(std::to_string(j))+","+TString( std::to_string(j * NumberSecondDiscrSlices[i]) ) + ")");
			//AltGen
			gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(\""+PlotNames[i]+"\","+TString(std::to_string(j))+","+TString( std::to_string(j * NumberSecondDiscrSlices[i]) ) +",false,true)");
			//GENIE FSI tweaks
			gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(\""+PlotNames[i]+"\","+TString(std::to_string(j))+","+TString( std::to_string(j * NumberSecondDiscrSlices[i]) ) +",false,false,true)");
			//GENIE Flag tweaks
			gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(\""+PlotNames[i]+"\","+TString(std::to_string(j))+","+TString( std::to_string(j * NumberSecondDiscrSlices[i]) ) +",false,false,false,true)");
			//GENIE Closure Test
			gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(\""+PlotNames[i]+"\","+TString(std::to_string(j))+","+TString( std::to_string(j * NumberSecondDiscrSlices[i]) ) +",false,false,false,false,true)");
			//GENIE Nuclear model tweaks
			gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(\""+PlotNames[i]+"\","+TString(std::to_string(j))+","+TString( std::to_string(j * NumberSecondDiscrSlices[i]) ) +",false,false,false,false,false,true)");
			//NuWro Comparisons out-of-the-box vs larsoft
			gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(\""+PlotNames[i]+"\","+TString(std::to_string(j))+","+TString( std::to_string(j * NumberSecondDiscrSlices[i]) ) +",false,false,false,false,false,false,true)");
			//Just data vs Nominal MC
			gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(\""+PlotNames[i]+"\","+TString(std::to_string(j))+","+TString( std::to_string(j * NumberSecondDiscrSlices[i]) ) +",false,false,false,false,false,false,false,true)");
			//All
			gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(\""+PlotNames[i]+"\","+TString(std::to_string(j))+","+TString( std::to_string(j * NumberSecondDiscrSlices[i]) ) +",false,false,false,false,false,false,false,false,false,true)");
			// w/ & w/o FSI + same for GiBUU
			gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(\""+PlotNames[i]+"\","+TString(std::to_string(j))+","+TString( std::to_string(j * NumberSecondDiscrSlices[i]) ) +",false,false,false,false,false,false,false,false,false,false,true)");							
			// G21 FSI
			gROOT->ProcessLine("ThreeDimWienerSVD_OverlayGenerators(\""+PlotNames[i]+"\","+TString(std::to_string(j))+","+TString( std::to_string(j * NumberSecondDiscrSlices[i]) ) +",false,false,false,false,false,false,false,false,false,false,false,true)");

		}

	}

}
