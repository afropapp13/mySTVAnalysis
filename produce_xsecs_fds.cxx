void produce_xsecs_fds() {

	gROOT->ProcessLine(".L ../myClasses/Util.C++");
	gROOT->ProcessLine(".L ../myClasses/WienerSVD.C++");
	gROOT->ProcessLine(".L fds_extract_xsec.cxx++");

	// Regular Overlay MC as the base MC, alternative BeamOn Fake Data samples

	gROOT->ProcessLine("fds_extract_xsec(\"Overlay9\",\"Overlay9NuWro\")");
	gROOT->ProcessLine("fds_extract_xsec(\"Overlay9\",\"NoTuneOverlay9\")");
	gROOT->ProcessLine("fds_extract_xsec(\"Overlay9\",\"TwiceMECOverlay9\")");

}
