# myEvents

# Extracting the nominal xsec predictions with statistical uncertainties

root -l script.C

# Extracting the detector variations xsecs

root -l script_Detector_Systematics.C

# Extracting the GEANT4 variations xsecs

root -l script_G4_Systematics.C

# Extracting the Genie variations xsecs

root -l script_Genie_Systematics.C

# Extracting the Flux variations xsecs

root -l script_Flux_Systematics.C


#################################################################################################################################


# Plotting the detector variation xsecs and storing the relevant systematics with respect to the CV sample

root -l Detector_Systematics.cpp

# Plotting the nominal xsecs, adding the 2% POT uncertainty and storing the relevant systematics with respect to the CV sample

root -l POT_Systematics.cpp

# Storing the GEANT4 uncertainties

root -l G4_Systematics.cpp

# Storing the Genie uncertainties

root -l Genie_Systematics.cpp

# Storing the Flux uncertainties

root -l Flux_Systematics.cpp


#################################################################################################################################


# Plotting the final results / total uncertainties / putting everything together

root -l Systematics.cpp


