% Generated by roxygen2 (4.1.1.9000): do not edit by hand
% Please edit documentation in R/RN_Analysis.R
\name{Occ.RN}
\alias{Occ.RN}
\title{Royle-Nichols Occupancy Model}
\usage{
Occ.RN(Y, Covariate.data = NULL, Abundance.formula = ~1,
  Detection.formula = ~1, N.limit = 10, n.chains = 4, n.adapt = 1000,
  n.iter = 1000, num.cores = 1, fit = NA, inits = NULL)
}
\arguments{
\item{Y}{A data frame with 5 columns - Species, Site, Observer, Effort, NumDetections}

\item{Covariate.data}{Data frame with the same number of rows as Y}

\item{Abundance.formula}{Model formula specifying the relationship between Abundance and the covariates}

\item{Detection.formula}{Model formula specifying the relationship between Detection probability and the covariates}

\item{N.limit}{Maximum conceivable number of individuals at a site. Making this too large slows down the model fitting, though.}

\item{n.chains}{How many MCMC chains to fit.}

\item{n.adapt}{How long should the burn-in phase be.}

\item{n.iter}{How many samples to take after the burn-in.}

\item{num.cores}{How many computer cores to use.}

\item{inits}{Initial values for the model chains.}
}
\value{
An object of type RN_OccupancyModel
}
\description{
Blah blah some description of the model
}

