% Generated by roxygen2 (4.1.1.9000): do not edit by hand
% Please edit documentation in R/OD_Analysis.R
\name{makeData.OD}
\alias{makeData.OD}
\title{Create simulated Occupancy data using the Occupancy/Detection model}
\usage{
makeData.OD(Design, Covariate.data = NULL, Occupancy.formula = ~1,
  Occupancy.params = c(0), Detection.formula = ~1,
  Detection.params = c(0))
}
\arguments{
\item{Design}{}

\item{Covariate.data}{}

\item{Occupancy}{model}

\item{Occupancy}{params}

\item{Detection}{model}

\item{Detection}{params}
}
\description{
Creates simulated data from the
Occupancy/Availability model. The primary use of
this function is for teaching the package.
}

