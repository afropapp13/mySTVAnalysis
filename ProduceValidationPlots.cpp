{

//gROOT->ProcessLine(".x Create2DPlots.cpp"); // Not needed, keep it for historical reasons

gROOT->ProcessLine(".L Create1DPlotsTotal.cpp++"); gROOT->ProcessLine("Create1DPlotsTotal()");

//gROOT->ProcessLine(".L DataDistributions.cpp++"); gROOT->ProcessLine("DataDistributions()");

//gROOT->ProcessLine(".x MigrationMatrices.cpp");

//gROOT->ProcessLine(".L XSection_Extraction.cpp++"); gROOT->ProcessLine("XSection_Extraction()");

//gROOT->ProcessLine(".L Systematics.cpp++"); gROOT->ProcessLine("Systematics()");

//gROOT->ProcessLine(".q");
};
