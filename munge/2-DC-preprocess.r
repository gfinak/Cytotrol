library(ProjectTemplate)
load.project()

panel <- "DC"
path_cytotrol <- "/loc/no-backup/ramey/Cytotrol"
path_FCS_files <- "/loc/no-backup/ramey/Cytotrol/DC FCS files/"
path_workspace <- "/loc/no-backup/ramey/Cytotrol/XML//CA_CytoTrol DC.xml"

# These are the markers that we will keep after the data have been preprocessed.
# NOTE: NHLBI does not include a Live marker
markers_of_interest <- c("FSC-A", "SSC-A", "CD56", "CD123", "CD11c", "CD16",
                         "Lineage", "CD14", "HLADR")

# Determines the centers from the FCS paths
centers <- sapply(strsplit(dir(path_FCS_files), split = " "), tail, n = 1)

# Parses each center's workspace and then extract its flowSet.
message("Parsing Lyoplate workspaces...")
fs_list <- sapply(centers, function(center) {
  message("Center: ", center)

  fcs_path <- file.path(path_FCS_files, paste(panel, center))
  ws <- openWorkspace(path_workspace)

  # Because there is a typo in the workspace for Kings, we manually set the
  # center in this case.
  if (center == "Kings") {
    gating_set <- parseWorkspace(ws, name = "KIngs", path = fcs_path, isNcdf = FALSE)
  } else {
    gating_set <- parseWorkspace(ws, name = center, path = fcs_path, isNcdf = FALSE)
  }
  closeWorkspace(ws)
  flow_set <- flowData(gating_set)

  pData(flow_set)$Center <- center
  varM <- varMetadata(phenoData(flow_set))
  varM[-1,] <- rownames(varM)[-1]
  varMetadata(phenoData(flow_set)) <- varM
  
  flow_set
}, simplify = FALSE)


# Swaps the channels and markers for the current 'flowSet' object. This ensures
# that we can 'rbind2' the 'GatingSetList' below because the stain names do not
# match otherwise.
message ("Swapping flowSet channels and markers")
fs_list <- sapply(centers, function(center) {
  message("Center: ", center)

  fsApply(fs_list[[center]], preprocess_flowframe,
          markers_keep = markers_of_interest)
}, simplify = FALSE)

# Merges the list of flowSet objects into a single flowSet object. This code is
# verbose but it circumvents an issue introduced recently in flowIncubator.
flow_set <- fs_list[[1]]
for (i in seq.int(2, length(fs_list))) {
  flow_set <- rbind2(flow_set, fs_list[[i]])
}

# Renames BIIR to Baylor to match Lyoplate data
pData(flow_set)$Center[pData(flow_set)$Center == "BIIR"] <- "Baylor"

# Creates GatingSet from the flowSet
gs_DC <- GatingSet(flow_set)

# Archives the results
save_gs(gs_DC, path = file.path(path_cytotrol, "gating-sets/gs-DC"))

