% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NonProbEst.R
\name{prop_estimation}
\alias{prop_estimation}
\title{Estimates the population proportion}
\usage{
prop_estimation(sample, weights, estimated_vars, class, N = NULL)
}
\arguments{
\item{sample}{A data frame containing the sample with the variables for which the means are to be calculated.}

\item{weights}{A vector of pre-calculated weights.}

\item{estimated_vars}{String vector specifying the variables in the sample to be estimated.}

\item{class}{String vector specifying which class (value) proportion is to be estimated in each variable. The \emph{i}-th element of this vector corresponds to the class of which proportion is desired to estimate of the \emph{i}-th variable of the vector specified in \code{estimated_vars}.}

\item{N}{An integer specifying the population size (optional).}
}
\value{
A vector with the corresponding estimations.
}
\description{
Estimates the proportion of a given class or classes for the specified variables measured in a sample given some pre-calculated weights.
}
\examples{
covariates = c("education_primaria", "education_secundaria")
data_propensities = propensities(sampleNP, sampleP, covariates)
psa_weights = sc_weights(data_propensities$convenience)

#The function will estimate the proportion of individuals
#with the 0 value in vote_pens and the 1 value in vote_pir
prop_estimation(sampleNP, psa_weights, c("vote_pens", "vote_pir"), c(0, 1))
}
