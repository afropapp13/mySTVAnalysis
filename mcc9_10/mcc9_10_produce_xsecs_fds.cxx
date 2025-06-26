void mcc9_10_produce_xsecs_fds() {

	gROOT->ProcessLine(".L ../../../generators/Util.C++");
	gROOT->ProcessLine(".L ../../../generators/WienerSVD.C++");
	gROOT->ProcessLine(".L mcc9_10_fds_extract_xsec.cxx++");

	// Regular Overlay MC as the base MC, alternative BeamOn Fake Data samples

	gROOT->ProcessLine("mcc9_10_fds_extract_xsec(\"mcc9_10_Overlay9\",\"mcc9_10_RSOverlay9\")");

}
