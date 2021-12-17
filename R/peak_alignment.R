peak_alignment <- function(input_path_pl, file_names_pl, RT_pl, mz_error, rt_tol, n_quantile, number_processing_cores) {
  L_PL <- length(file_names_pl)
  imzRTXcol_main <- do.call(rbind, lapply(1:L_PL, function(i) {
    peaklist <- loadRdata(paste0(input_path_pl, "/", file_names_pl[i]))
    cbind(rep(i, nrow(peaklist)), peaklist[, 8], RT_pl[[i]], peaklist[, Xcol = 4], 1:nrow(peaklist))
  }))
  imzRTXcol_main <- imzRTXcol_main[order(imzRTXcol_main[, 2], decreasing = TRUE), ]
  ##
  MZ_Q <- quantile(imzRTXcol_main[, 2], probs = c(1:n_quantile)/n_quantile)
  MZ_Q_boundaries <- cbind(c(0, MZ_Q[1 :(n_quantile - 1)]), MZ_Q)
  MZ_Q_boundaries[, 1] <- MZ_Q_boundaries[, 1] - mz_error*1.5
  MZ_Q_boundaries[, 2] <- MZ_Q_boundaries[, 2] + mz_error*1.5
  #
  FeatureTable_Main_call <- function (q) {
    x_Q <- which(imzRTXcol_main[, 2] >= MZ_Q_boundaries[q, 1] &
                   imzRTXcol_main[, 2] <= MZ_Q_boundaries[q, 2])
    imzRTXcol <- imzRTXcol_main[x_Q, ]
    imzRTXcol <- imzRTXcol[order(imzRTXcol[, 4], decreasing = TRUE), ]
    N_imzRTXcol <- length(x_Q)
    FeatureTable <- matrix(rep(0, (L_PL + 2)*N_imzRTXcol), nrow = N_imzRTXcol)
    counter <- 1
    i <- 1
    repeat {
      FeatureTable[counter, 1] <- imzRTXcol[i, 2]
      FeatureTable[counter, 2] <- imzRTXcol[i, 3]
      FeatureTable[counter, (imzRTXcol[i, 1] + 2)] <- imzRTXcol[i, 5]
      x <- which(abs(imzRTXcol[i, 2] - imzRTXcol[, 2]) <= mz_error &
                   abs(imzRTXcol[i, 3] - imzRTXcol[, 3]) <= rt_tol &
                   imzRTXcol[i, 1] != imzRTXcol[, 1] & imzRTXcol[, 1] != 0)
      if (length(x) > 0) {
        iSamples <- imzRTXcol[x, 1]
        if (length(iSamples) != length(unique(iSamples))) {
          xi <- table(iSamples)
          xii <- which(xi > 1)
          Rx <- as.numeric(names(xii))  # Repeated peaks in the same sample
          A <- cbind(imzRTXcol[x, ], x)
          xiii <- sapply(1:length(Rx), function(j) {
            xj <- which(A[, 1] == Rx[j])
            xjj <- which.min(abs(A[xj, 3] - imzRTXcol[i, 3]))
            A[xj[-xjj], 6]
          })
          x <- setdiff(x, xiii)
          iSamples <- imzRTXcol[x, 1]
        }
        for (j in 1:length(iSamples)) {
          FeatureTable[counter, iSamples[j] + 2] <- imzRTXcol[x[j], 5]
        }
        imzRTXcol[x, 1] <- 0
        imzRTXcol[i, 1] <- 0
      }
      if (i < N_imzRTXcol) {
        i <- i + which(imzRTXcol[(i + 1):N_imzRTXcol, 1] != 0)[1]
        if (is.na(i)) {
          break
        }
      } else {
        break
      }
      counter <- counter + 1
    }
    FeatureTable[1:counter, ]
  }
  ##
  osType <- Sys.info()[['sysname']]
  if (osType == "Linux") {
    FeatureTable_Main <- do.call(rbind, mclapply(1:n_quantile, function (q) {
      FeatureTable_Main_call(q)
    }, mc.cores = number_processing_cores))
  }
  ##
  if(osType == "Windows") {
    cl <- makeCluster(number_processing_cores)
    registerDoSNOW(cl)
    FeatureTable_Main <- foreach(q = 1:n_quantile, .combine="rbind", .verbose = FALSE) %dopar% {
      FeatureTable_Main_call(q)
    }
    stopCluster(cl)
  }
  ## To resolve redundant peaks in the peak matrix table
  L_PL2 <- L_PL + 2
  L_PL3 <- L_PL + 3
  x_s <- sapply(1:dim(FeatureTable_Main)[1], function(i) {
    length(which(FeatureTable_Main[i, 3:L_PL2] > 0))
  })
  FeatureTable_Main <- cbind(x_s, FeatureTable_Main)
  FeatureTable_Main <- FeatureTable_Main[order(FeatureTable_Main[, 1], decreasing = TRUE), ]
  ##
  progressBARboundaries <- txtProgressBar(min = 1, max = dim(FeatureTable_Main)[1], initial = 1, style = 3)
  for (i in 1:dim(FeatureTable_Main)[1]) {
    setTxtProgressBar(progressBARboundaries, i)
    x_c <- which(abs(FeatureTable_Main[i, 2] - FeatureTable_Main[, 2]) <= mz_error &
                   abs(FeatureTable_Main[i, 3] - FeatureTable_Main[, 3]) <= rt_tol &
                   FeatureTable_Main[i, 1] != 0)
    L_x_c <- length(x_c)
    if (L_x_c > 1) {
      x_diff <- setdiff(x_c, i)
      if (FeatureTable_Main[i, 1] < L_PL) {
        table_c <- do.call(rbind, lapply(1:L_x_c, function(j) {
          FeatureTable_Main[x_c[j], 1:L_PL3]
        }))
        x_table_main0 <- which(table_c[1, ] == 0)
        for (j in x_table_main0) {
          x_non0 <- which(table_c[, j] > 0)
          if (length(x_non0) > 0) {
            if (length(x_non0) > 1) {
              x_min <- which.min(abs(table_c[1, 3] - table_c[x_non0, 3]))
              x_non0 <- x_non0[x_min[1]]
            }
            table_c[1, j] <- table_c[x_non0, j]
          }
        }
        FeatureTable_Main[i, 4:L_PL3] <- table_c[1, 4:L_PL3]
      }
      FeatureTable_Main[x_diff, ] <- 0
    }
  }
  close(progressBARboundaries)
  x_non0 <- which(FeatureTable_Main[, 1] != 0)
  FeatureTable_Main <- FeatureTable_Main[x_non0, ]
  FeatureTable_Main <- FeatureTable_Main[, -1]
  FeatureTable_Main <- FeatureTable_Main[order(FeatureTable_Main[, 1], decreasing = FALSE), ]
  return(FeatureTable_Main)
}
