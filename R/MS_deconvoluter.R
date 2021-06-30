MS_deconvoluter <- function(MassSpec_file, MS_level = 1) {
  xfile <- openMSfile(MassSpec_file)
  peakTable <- header(xfile) # this gets table of details for each spectra
  spectraList <- spectra(xfile) # this gets the spectra values
  x_MS <- which(peakTable$peaksCount > 0 & peakTable$msLevel == MS_level) # peaks from soft ionization channel ## some files may not have data in the column re-calibration period.
  spectraList <- spectraList[x_MS]
  peakTable <- peakTable[x_MS, ]
  RetentionTime <- as.data.frame(peakTable[, 7], drop = FALSE)/60  # Retention times in minute
  RetentionTime <- as.matrix(RetentionTime)
  outputer <- list(spectraList, RetentionTime)
  return(outputer)
}
