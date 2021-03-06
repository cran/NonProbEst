% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NonProbEst.R
\name{mean_estimation}
\alias{mean_estimation}
\title{Estimates the population means}
\usage{
mean_estimation(sample, weights, estimated_vars, N = NULL)
}
\arguments{
\item{sample}{A data frame containing the sample with the variables for which the means are to be calculated.}

\item{weights}{A vector of pre-calculated weights.}

\item{estimated_vars}{String vector specifying the variables in the sample to be estimated.}

\item{N}{An integer specifying the population size (optional).}
}
\value{
A vector with the corresponding estimations.
}
\description{
Estimates the means for the specified variables measured in a sample given some pre-calculated weights.
}
\examples{
covariates = c("education_primaria", "education_secundaria")
data_propensities = propensities(sampleNP, sampleP, covariates)
psa_weights = sc_weights(data_propensities$convenience)
mean_estimation(sampleNP, psa_weights, c("vote_pens"))
}
