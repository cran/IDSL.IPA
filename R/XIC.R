XIC <- function(spectraList.xic, scan_number_start, mz_target, mass_accuracy_xic) {
  chrom_ScN_Int <- do.call(rbind, lapply(1:length(spectraList.xic), function(t) {
    PEAKS <- spectraList.xic[[t]]
    PEAKSt <- c(t, 0, 0)
    x <- which(abs(PEAKS[ ,1] - mz_target) <= mass_accuracy_xic)
    if (length(x) > 0) {
      if (length(x) > 1) {
        x2 <- which.min(abs(PEAKS[x, 1] - mz_target))
        x <- x[x2[1]]
      }
      PEAKSt <- c(t, PEAKS[x, 1], PEAKS[x, 2]) # c(Scan number, m/z, Intensity)
    }
    PEAKSt
  }))
  chrom_ScN_Int[, 1] <- chrom_ScN_Int[, 1] + scan_number_start - 1
  return(chrom_ScN_Int)
}
