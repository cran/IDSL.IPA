\name{gaussianity_measurement}
\alias{gaussianity_measurement}
\title{
gaussianity measurement
}
\description{
This module measures gaussianity of chromatographic peak using Pearson correlation coefficients (\eqn{\rho}) at top 80 percent of peak.
}
\usage{
gaussianity_measurement(RT, Int, BL, gauge = 0.8)
}
\arguments{
  \item{RT}{
a vector of retention times of the chromatographic peak.
}
  \item{Int}{
a vector of intensities of the chromatographic peak.
}
  \item{BL}{
a vector of baseline of the chromatographic peak.
}
  \item{gauge}{
represents the gauge height of peak for Gaussianity measurement.
}
}
  \value{
Gaussianity of the chromatographic peak.
}
\examples{
data("peak_spline")
RT <- peak_spline[, 1]
Int <- peak_spline[, 2]
BL <- peak_spline[, 3]
gaussianity_measurement(RT, Int, BL, gauge = 0.8)
}
