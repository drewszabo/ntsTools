#' Export MS2 spectra from patRoon to GNPS format
#' This function exports MS2 spectra from patRoon v3.0 to a format compatible with GNPS. It generates a feature table and an MGF file containing the MS2 spectra, which can be used for molecular networking in GNPS.
#' @param fGroups The feature groups object from patRoon (default is fGroups).
#' @param mslists The MS lists object from patRoon (default is mslists).
#' @param polarity The polarity of the MS data, either "positive" or "negative" (default is c("positive", "negative")).
#' @param path The file path where the output files will be saved (default is the current working directory).
#' @return A list containing the feature table and the MGF file path for GNPS.
#' @importFrom xcms featureDefinitions featureValues
#' @importFrom patRoon as.data.table
#' @export
generateGNPSInfo <- function(fGroups = fGroups, mslists = mslists, polarity = c("positive", "negative"), path) {
  # Generate feature table
  featuresDef <- featureDefinitions(fGroups@xdata)
  featuresIntensities <- featureValues(fGroups@xdata, value = "into")
  peakID <- as.character(colnames(fGroups@groups))
  
  dataTable <- merge(featuresDef, featuresIntensities, by = 0, all = TRUE)
  dataTable <- dataTable[, !(colnames(dataTable) %in% c("peakidx"))]
  dataTable$peak_ID <- peakID
  
  write.table(dataTable, paste0(path, "output_gnps.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  resultsmslists <- patRoon::as.data.table(mslists)
  groups <- unique(resultsmslists$group)
  
  fileConn <- file(paste0(path, "output_gnps.mgf"), "w")
  writeLines(paste0("COM=Exported by ", Sys.getenv("USERNAME"), " on ", format(Sys.time(), "%d%m%y %H%M")), fileConn)
  
  for (peak_ID in groups) {
    feature_list <- dataTable[dataTable$peak_ID == peak_ID, ]
    peaks <- resultsmslists[group == ..peak_ID]
    MS1peaks <- peaks[type == "MS", .(mz, intensity)]
    MS2peaks <- peaks[type == "MSMS", .(mz, intensity)]
    
    ret <- round(feature_list$rtmed, 2)
    mz <- round(feature_list$mzmed, 4)
    numPeaksMS1 <- nrow(MS1peaks)
    numPeaksMS2 <- nrow(MS2peaks)
    
    # MS1 Block
    writeLines(c("BEGIN IONS",
                 paste0("FEATURE_ID=", feature_list$Row.names),
                 paste0("PEAK_ID=", peak_ID),
                 "MSLEVEL=1",
                 paste0("RTINSECONDS=", ret),
                 paste0("PEPMASS=", mz),
                 paste0("CHARGE=", ifelse(polarity == "positive", "1+", "1-")),
                 paste0("SCANS=", as.integer(sub("^FT", "", feature_list$Row.names))),
                 paste0("Num peaks=", numPeaksMS1),
                 apply(MS1peaks, 1, function(x) paste(round(x[1], 6), round(x[2], 0), sep = " ")),
                 "END IONS", ""), fileConn)
    
    # MS2 Block
    writeLines(c("BEGIN IONS",
                 paste0("FEATURE_ID=", feature_list$Row.names),
                 paste0("PEAK_ID=", peak_ID),
                 "MSLEVEL=2",
                 paste0("RTINSECONDS=", ret),
                 paste0("PEPMASS=", mz),
                 paste0("CHARGE=", ifelse(polarity == "positive", "1+", "1-")),
                 paste0("SCANS=", as.integer(sub("^FT", "", feature_list$Row.names))),
                 paste0("Num peaks=", numPeaksMS2),
                 apply(MS2peaks, 1, function(x) paste(round(x[1], 6), round(x[2], 0), sep = " ")),
                 "END IONS", ""), fileConn)
  }
  
  close(fileConn)
}

#' Generate MGF file for SIRIUS from patRoon data
#' This function generates an MGF file for SIRIUS from patRoon data. It creates a feature table and an MGF file containing the MS2 spectra, which can be used for compound identification in SIRIUS.
#' @param fGroups The feature groups object from patRoon (default is fGroups).
#' @param mslists The MS lists object from patRoon (default is mslists).
#' @param polarity The polarity of the MS data, either "positive" or "negative" (default is c("positive", "negative")).
#' @param path The file path where the output MGF file will be saved (default is the current working directory).
#' @return The file path to the generated MGF file for SIRIUS.
#' @importFrom xcms featureDefinitions featureValues
#' @importFrom patRoon as.data.table
#' @import data.table
#' @export
generateSIRIUSmgf <- function(fGroups = fGroups, mslists = mslists, polarity = c("positive", "negative"), path) {
  library(data.table)
  
  # Match polarity argument safely
  polarity <- match.arg(polarity)
  charge_str <- if (polarity == "positive") "1+" else "1-"
  
  # 1. Efficiently generate and format the feature data table
  featuresDef <- as.data.table(xcms::featureDefinitions(fGroups@xdata))
  featuresIntensities <- xcms::featureValues(fGroups@xdata, value = "into")
  peakID <- as.character(colnames(fGroups@groups))
  
  dataTable <- as.data.table(featuresIntensities)
  dataTable <- cbind(featuresDef, dataTable)
  dataTable[, peak_ID := peakID]
  
  # Keep only the columns we actually need for the loop and set keys for instant lookups
  feature_lookup <- dataTable[, .(peak_ID, rtmed = round(rtmed, 2), mzmed = round(mzmed, 4))]
  setkey(feature_lookup, peak_ID)
  
  # 2. Process MS lists using data.table splitting
  resultsmslists <- patRoon::as.data.table(mslists)
  setDT(resultsmslists) # Ensure it's treated as a data.table
  
  # Pre-round columns to vectorize computation before splitting
  resultsmslists[, mz_str := sprintf("%.6f", mz)]
  resultsmslists[, int_str := sprintf("%.0f", intensity)]
  resultsmslists[, peak_line := paste(mz_str, int_str)]
  
  # Split the massive list into a named list of data.tables by group (very fast)
  ms_split <- split(resultsmslists, by = "group", keep.by = FALSE)
  groups_list <- names(ms_split)
  
  # 3. Open file connection
  out_file <- file.path(path, "output_sirius.mgf")
  fileConn <- file(out_file, "w")
  
  # Write Header
  writeLines(paste0("COM=Exported by ", Sys.getenv("USERNAME"), " on ", format(Sys.time(), "%d%m%y %H%M")), fileConn)
  
  # 4. Process efficiently group by group
  for (this_group in groups_list) {
    # Instant O(1) key-based lookup for features
    f_row <- feature_lookup[.(this_group)]
    if (nrow(f_row) == 0 || is.na(f_row$mzmed)) next
    
    peaks_subset <- ms_split[[this_group]]
    
    # Separate MS1 and MS2 peak text lines using vector indexing
    ms1_lines <- peaks_subset[type == "MS", peak_line]
    ms2_lines <- peaks_subset[type == "MSMS", peak_line]
    
    # Construct MS1 Block
    ms1_block <- c(
      "BEGIN IONS",
      paste0("FEATURE_ID=", this_group),
      paste0("PEPMASS=", f_row$mzmed),
      "MSLEVEL=1",
      paste0("RTINSECONDS=", f_row$rtmed),
      paste0("CHARGE=", charge_str),
      ms1_lines,
      "END IONS", ""
    )
    
    # Construct MS2 Block
    ms2_block <- c(
      "BEGIN IONS",
      paste0("FEATURE_ID=", this_group),
      paste0("PEPMASS=", f_row$mzmed),
      "MSLEVEL=2",
      paste0("RTINSECONDS=", f_row$rtmed),
      paste0("CHARGE=", charge_str),
      ms2_lines,
      "END IONS", ""
    )
    
    # Single write operation per group
    writeLines(c(ms1_block, ms2_block), fileConn)
  }
  
  close(fileConn)
  return(out_file)
}