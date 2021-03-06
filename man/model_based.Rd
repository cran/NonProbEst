% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NonProbEst.R
\name{model_based}
\alias{model_based}
\title{Calculates a model based estimation}
\usage{
model_based(
  sample_data,
  full_data,
  covariates,
  estimated_var,
  estimate_mean = FALSE,
  positive_label = NULL,
  algorithm = "glm",
  proc = NULL,
  ...
)
}
\arguments{
\item{sample_data}{Data frame containing the sample.}

\item{full_data}{Data frame containing all the individuals contained in the population.}

\item{covariates}{String vector specifying the common variables to use for training.}

\item{estimated_var}{String specifying the variable to estimate.}

\item{estimate_mean}{Boolean specifying whether the mean estimation should be returned. Otherwise, the total estimation is returned by default.}

\item{positive_label}{String specifying the label to be considered positive if the estimated variable is categorical. Leave it as the default NULL otherwise.}

\item{algorithm}{A string specifying which classification or regression model to use (same as caret's method).}

\item{proc}{A string or vector of strings specifying if any of the data preprocessing techniques available in \link[caret]{train} function from `caret` package should be applied to data prior to the propensity estimation. By default, its value is NULL and no preprocessing is applied.}

\item{...}{Further parameters to be passed to the \link[caret]{train} function.}
}
\value{
The population total estimation (or mean if specified by the `estimate_mean` parameter).
}
\description{
It uses the model based estimator. The idea in order to estimate the population total is to add the sample responses and the predicted responses for the individuals not contained in the sample. See for example Valliant et al. (2000).
}
\details{
Training of the models is done via the `caret` package. The algorithm specified in \code{algorithm} must match one of the names in the list of algorithms supported by `caret`.
}
\examples{
#Simple example
covariates = c("education_primaria", "education_secundaria",
   "age", "sex", "language")
if (is.numeric(sampleNP$vote_gen))
   sampleNP$vote_gen = factor(sampleNP$vote_gen, c(0, 1), c('F', 'T'))
model_based(sampleNP, population, covariates,
   "vote_gen", positive_label = 'T', algorithm = 'glmnet')
}
\references{
Valliant, R., Dorfman, A. H., & Royall, R. M. (2000) \emph{Finite population sampling and inference: a prediction approach.} Wiley, New York.
}
