carbon_isotopes_explorer <- function(spectraList, int_threshold, mass_accuracy_13c, max_R13C) {
  rR13C <- max_R13C/100
  NumScans <- length(spectraList)
  spec_scan <- do.call(rbind, lapply(1:NumScans, function(t) {
    ##
    Spec <- spectraList[[t]]
    L_Spec <- nrow(Spec)
    if (L_Spec > 0) {
      spec_scan_j <- matrix(rep(0, L_Spec*5), ncol = 5)
      counter_j <- 0
      ##
      mz <- Spec[, 1]
      int <- Spec[, 2]
      for (j in 1:L_Spec) {
        if ((int[j] >= int_threshold) & (mz[j] != 0)) {        # Intensity threshold in each scan
          x13C <- which(abs(mz[j] - (mz - 1.00335484)) <= mass_accuracy_13c)
          L_x13C <- length(x13C)
          if (L_x13C > 0) {
            if (L_x13C > 1) {
              x13C_min <- which.min(abs(mz[j] - (mz[x13C] - 1.00335484)))
              x13C <- x13C[x13C_min]
              ##
              if (length(x13C) >= 2) {
                xR13C <- which(int[x13C]/int[j] <= rR13C)
                x13C <- x13C[xR13C[1]]
              }
            }
            if (int[x13C]/int[j] <= rR13C) {
              counter_j <- counter_j + 1
              spec_scan_j[counter_j, ] <- c(mz[j], int[j], t, mz[x13C], int[x13C])
              mz[x13C] <- 0
            }
          }
        }
      }
      ##
      if (counter_j > 0) {
        spec_scan_j <- spec_scan_j[1:counter_j, ]
      }
    }
  }))
  rownames(spec_scan) <- c()
  spec_scan <- spec_scan[order(spec_scan[, 2], decreasing = TRUE), ]   # Sort spec_scan rows by their intensity
  return(spec_scan)
}
