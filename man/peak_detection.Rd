\name{peak_detection}
\alias{peak_detection}
\title{
peak detection
}
\description{
This function detects separated chromatographical peaks on the chromatogram.
}
\usage{
peak_detection(int)
}
\arguments{
  \item{int}{
a vector of intensities of the chromatogram.
}
}
\value{
A matrix of 2 columns. Each row indicates peak boundary indices on the 'int' vector.
}
\examples{
data(chromatogram_builder)
int <- chromatogram_builder$SmoothedChromatogram
peak_detection(int)
}
