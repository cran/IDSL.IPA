IPA_MSdeconvoluter <- function(HRMS_path, MSfile, MS_level = 1) {
  p2l <- IDSL.MXP::peak2list(HRMS_path, MSfile)
  peakTable <- p2l[["peakTable"]]
  spectraList <- p2l[["spectraList"]]
  x_MS <- which(peakTable$peaksCount > 0 & peakTable$msLevel == MS_level) ## some files may not have data in the column re-calibration period.
  spectraList <- spectraList[x_MS]
  peakTable <- peakTable[x_MS, ]
  RetentionTime <- peakTable$retentionTime  # Retention times in minute
  MS_polarity <- ifelse(peakTable$polarity[x_MS[1]] == 1, "+", "-")
  outputer <- list(spectraList, RetentionTime, MS_polarity)
  return(outputer)
}
