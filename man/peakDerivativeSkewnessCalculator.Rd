\name{peakDerivativeSkewnessCalculator}
\alias{peakDerivativeSkewnessCalculator}
\title{
Peak Derivative Skewness Calculator 
}
\description{
This function calculates skewness of a chromatographic peak using first order degree of numerical differentiation.
}
\usage{
peakDerivativeSkewnessCalculator(rt, int)
}
\arguments{
  \item{rt}{
a vector representing retention times of the chromatographic peak.
}
  \item{int}{
a vector representing intensities of the chromatographic peak.
}
}
\value{
Skewness of a chromatographic peak. 1 is for very symmetric peak. Minimum is 0 from this function.
}
\examples{
data(peak_spline)
rt <- peak_spline[, 1]
int <- peak_spline[, 2]
peakDerivativeSkewnessCalculator(rt, int)
}
