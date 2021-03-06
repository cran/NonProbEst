% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NonProbEst.R
\name{valliant_weights}
\alias{valliant_weights}
\title{Calculates Valliant weights}
\usage{
valliant_weights(propensities)
}
\arguments{
\item{propensities}{A vector with the propensities associated to the elements of the convenience sample.}
}
\value{
A vector with the corresponding weights.
}
\description{
Computes weights from propensity estimates using the 1/pi_i formula introduced in Valliant (2019).
}
\details{
The function takes the vector of propensities \eqn{\pi(x)} and calculates the weights to be applied in the Hajek estimator using the formula that can be found in Valliant (2019). For an individual \emph{i}, weight is calculated as follows:
\deqn{w_i = 1/\pi_i (x)}
}
\examples{
covariates = c("education_primaria", "education_secundaria")
data_propensities = propensities(sampleNP, sampleP, covariates)
valliant_weights(data_propensities$convenience)
}
\references{
Valliant, R. (2019). \emph{Comparing Alternatives for Estimation from Nonprobability Samples}. Journal of Survey Statistics and Methodology, smz003, \url{https://doi.org/10.1093/jssam/smz003}
}
