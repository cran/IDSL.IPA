\name{MS_deconvoluter}
\alias{MS_deconvoluter}
\title{
MS deconvoluter
}
\description{
This function deconvolute mass spectrometry files into a list of mass spectrals and a vector of retention times.
}
\usage{
MS_deconvoluter(MassSpec_file, MS_level = 1)
}
\arguments{
  \item{MassSpec_file}{
mass spectromtery file.
}
  \item{MS_level}{
MS level to extract information.
}
}
\value{
  \item{spectraList}{
 a list of mass spectra.
}
  \item{RetentionTime}{
 a vector of retention times for scan numbers.
}
  \item{MS_polarity}{
mass spectromtery ionization mode (+/-)
}
}
