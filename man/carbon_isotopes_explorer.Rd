\name{carbon_isotopes_explorer}
\alias{carbon_isotopes_explorer}
\title{
Carbon Isotopes Explorer
}
\description{
This function isolates 12C/13C isotopologue pairs in high-resolution mass
spectral datasets
}
\usage{
carbon_isotopes_explorer(spectraList, int_threshold, mass_accuracy_13c,
max_R13C)
}
\arguments{
  \item{spectraList}{
a list of mass spectra in each chromatogram scan.
}
  \item{int_threshold}{
a value to represent intensity threshold at each chromatogram scan.
}
  \item{mass_accuracy_13c}{
a mass error to detect 13C isotoplogues.
}
  \item{max_R13C}{
a maximum allowed value of R13C for 12C/13C isotopologue pairs in each
chromatogram scan.
}
}
\value{
A matrix consists of 5 columns. The column contents are the m/z of 12C
isotoplogues, intensity of 12C isotoplogues, scan number (t), m/z of 13C
isotoplogues, and intensity of 13C isotoplogues, respectively.
}
