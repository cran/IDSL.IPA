carbon_isotopes_explorer <- function (spectraList, int_threshold, mass_accuracy_13c, max_R13C) {
  NumScans <- length(spectraList)
  spec_scan <- do.call(rbind, lapply(1:NumScans, function(t) {
    spec_scan_j <- c()
    Spec <- spectraList[[t]]
    if (length(Spec) > 0) {
      mz <- Spec[, 1]
      int <- Spec[, 2]
      for (j in 1:length(mz)) {
        if (int[j] >= int_threshold) {        # Intensity threshold in each scan
          x13C <- which(abs(mz[j] - (mz - 1.00335484)) <= mass_accuracy_13c)
          L_x13C <- length(x13C)
          if (L_x13C > 0) {
            if (L_x13C > 1) {
              x13C_min <- which.min(abs(mz[j] - (mz[x13C] - 1.00335484)))
              x13C <- x13C[x13C_min[1]]
            }
            if (int[x13C]/int[j] <= max_R13C/100) {
              spec_scan_x <- c(mz[j], int[j], t, mz[x13C], int[x13C])
              spec_scan_j <- rbind(spec_scan_j, spec_scan_x)
              mz[x13C] <- 0
            }
          }
        }
      }
    }
    spec_scan_j
  }))
  rownames(spec_scan) <- c()
  spec_scan <- spec_scan[order(spec_scan[, 2], decreasing = TRUE), ]   # Sort spec_scan rows by their intensity
  return(spec_scan)
}
