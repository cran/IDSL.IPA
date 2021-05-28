\name{snr_xcms}
\alias{snr_xcms}
\title{
SNR xcms
}
\description{
This function calculates S/N values using a method suggested in the xcms paper (Tautenhahn, 2008).
}
\usage{
snr_xcms(int)
}
\arguments{
  \item{int}{
a vector of intensities corresponding to the vector of retention times for the
chromatographic peak.
}
}
\value{
S/N value
}
\examples{
data(peak_spline)
int <- peak_spline[, 2]
snr_xcms(int)
}
