IPA_Workflow <- function(spreadsheet) {
  PARAM <- IPA_xlsxAnalyzer(spreadsheet)
  if (length(PARAM) > 0) {
    if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0001'), 2]) == "yes") {
      IPA_PeakAnalyzer(PARAM)
    }
    if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0002'), 2]) == "yes") {
      IPA_PeakAlignment(PARAM)
    } else {
      x0004 <- PARAM[which(PARAM[, 1] == 'PARAM0004'), 2]
      x0045 <- PARAM[which(PARAM[, 1] == 'PARAM0045'), 2]
      if (!is.na(x0004) & !is.na(x0045)) {
        if (tolower(x0004) == "yes" & tolower(x0045) == "yes") {
          IPA_PeakAlignment(PARAM)
        }
      }
    }
    if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0003'), 2]) == "yes") {
      IPA_GapFiller(PARAM)
    }
    if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0004'), 2]) == "yes") {
      if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0046'), 2]) == "yes") {
        IPA_CompoundsAnnotation(PARAM)
      }
      if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0047'), 2]) == "yes") {
        IPA_PeaklistAnnotation(PARAM)
      }
    }
    print("Completed computation successfully!")
  }
}
