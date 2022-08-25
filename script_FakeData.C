void script_FakeData() {

	gROOT->ProcessLine(".L ../../myClasses/Util.C++");
	gROOT->ProcessLine(".L ../../myClasses/WienerSVD.C++");
	gROOT->ProcessLine(".L FakeData_WienerSVD_XSection_Extraction.cpp++");

	// Regular Overlay MC as the base MC, alternative BeamOn Fake Data samples

	gROOT->ProcessLine("FakeData_WienerSVD_XSection_Extraction(\"Overlay9\",\"Overlay9NuWro\")");
	gROOT->ProcessLine("FakeData_WienerSVD_XSection_Extraction(\"Overlay9\",\"NoTuneOverlay9\")");
	gROOT->ProcessLine("FakeData_WienerSVD_XSection_Extraction(\"Overlay9\",\"GENIEv2Overlay9\")");	
	gROOT->ProcessLine("FakeData_WienerSVD_XSection_Extraction(\"Overlay9\",\"TwiceMECOverlay9\")");

	// mixed study: alternative MC (GENIEv2,NoTuneOverlay9,TwiceMECOverlay9) & fake data (NuWro)
	// need corresponding response & total covariance matrix

	gROOT->ProcessLine("FakeData_WienerSVD_XSection_Extraction(\"GENIEv2Overlay9\",\"Overlay9NuWro\")");
	gROOT->ProcessLine("FakeData_WienerSVD_XSection_Extraction(\"NoTuneOverlay9\",\"Overlay9NuWro\")");
	gROOT->ProcessLine("FakeData_WienerSVD_XSection_Extraction(\"TwiceMECOverlay9\",\"Overlay9NuWro\")");

}
