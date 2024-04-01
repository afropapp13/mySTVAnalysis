root -b run_cv.cxx
root -b run_det.cxx
root -b run_g4.cxx
root -b run_xsec.cxx
root -b run_flux.cxx
root -b run_mc_stat.cxx
root -b run_nuwro.cxx

#################################################################################################################################

root -b systematics.cxx
root -b merge_covariances.cxx
root -b produce_xsecs.cxx

#################################################################################################################################

# Fake data studies

# We need the stat & xsec uncertainties only

root -b
.L fds_stat_covariance_matrices.cxx
fds_stat_covariance_matrices("Stat","Overlay9","Overlay9NuWro","ExtBNB9","OverlayDirt9")
fds_stat_covariance_matrices("Stat","Overlay9","NoTuneOverlay9","ExtBNB9","OverlayDirt9")
fds_stat_covariance_matrices("Stat","Overlay9","TwiceMECOverlay9","ExtBNB9","OverlayDirt9")

# Merge the alternative MC covariances
root -b
.L merge_covariances.cxx++
merge_covariances("Overlay9","Overlay9NuWro","Overlay9NuWro")
merge_covariances("Overlay9","NoTuneOverlay9","NoTuneOverlay9")
merge_covariances("Overlay9","TwiceMECOverlay9","TwiceMECOverlay9") 

#################################################################################################################################

# Fake data studies with Wiener SVD
root -b produce_xsecs_fds.cxx

# xsec uncertainties
root -b plot_xsec_unc.cxx

#################################################################################################################################

cd ../myEvents

root -b topological_breakdown.cxx
root -b interaction_breakdown.cxx

