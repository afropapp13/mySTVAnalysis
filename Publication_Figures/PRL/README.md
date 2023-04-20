root -b script_PRL_Fig1.C
root -b script_PRL_Fig2.C
root -b script_PRL_Fig3.C

root -b script_PRL_Supp_CosThetaMu.C
root -b script_PRL_Supp_VertexZ.C


./PRL_Copy_SuppMat_Plots.sh
./PRL_Copy_PreSelection_Plots.sh
./PRL_Copy_FakeDataStudies.sh
./PRL_InteractionBreakDown.sh

root -b DataRelease.cxx
root -b BinScheme.cxx
