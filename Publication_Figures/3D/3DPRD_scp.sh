export OutPutDir=/home/afroditi/Dropbox/Apps/Overleaf/General\ Imbalance\ Variables\ for\ Measuring\ Nuclear\ Effects\ and\ Demonstration\ with\ MicroBooNE\ Data/figures/gen/Afro/
export InPutDir=/uboone/data/users/apapadop/FlatTTreePlots/

scp apapadop@uboonegpvm07.fnal.gov:"${InPutDir}"*.pdf "${OutPutDir}"
