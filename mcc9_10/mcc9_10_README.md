root -b mcc9_10_run_cv.cxx
/*root -b mcc9_10_run_det.cxx*/ // we don't have det vars yet
root -b mcc9_10_run_g4.cxx
root -b mcc9_10_run_xsec.cxx
root -b mcc9_10_run_flux.cxx
root -b mcc9_10_run_mc_stat.cxx
/*root -b mcc9_10_run_nuwro.cxx*/ // we don't have the nuwro alternative samples yet

#################################################################################################################################

root -b mcc9_10_systematics.cxx
root -b mcc9_10_merge_covariances.cxx
root -b mcc9_10_produce_xsecs.cxx

root -b mcc9_10_event_rate_systematics.cxx
root -b mcc9_10_event_rate_merge_covariances.cxx

#################################################################################################################################

# Fake data studies

# We need the stat & xsec uncertainties only

root -b
.L mcc9_10_fds_stat_covariance_matrices.cxx
mcc9_10_fds_stat_covariance_matrices("Stat","mcc9_10_Overlay9","mcc9_10_RSOverlay9","mcc9_10_ExtBNB9","mcc9_10_OverlayDirt9")

# Merge the alternative MC covariances
root -b
.L mcc9_10_merge_covariances.cxx++
mcc9_10_merge_covariances("mcc9_10_Overlay9","mcc9_10_RSOverlay9","mcc9_10_RSOverlay9")

#################################################################################################################################

# Fake data studies with Wiener SVD
root -b mcc9_10_produce_xsecs_fds.cxx

# xsec uncertainties
root -b mcc9_10_plot_xsec_unc.cxx

#################################################################################################################################

cd ../../myEvents/mcc9_10

root -b mcc9_10_topological_breakdown.cxx
root -b mcc9_10_interaction_breakdown.cxx