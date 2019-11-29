#' Calculates each sample propensities
#'
#' Given a convenience sample and a reference sample, computes estimates on the propensity to participate in the convenience sample based on classification models to be selected by the user.
#'
#' Training of the propensity estimation models is done via the `caret` package. The algorithm specified in \code{algorithm} must match one of the names in the list of algorithms supported by `caret`. Case weights are used to balance classes (for models that accept them).
#' The smoothing formula for propensities avoids mathematical irregularities in the calculation of sample weight when an estimated propensity is 0 or 1. Further details can be found in Buskirk and Kolenikov (2015).
#'
#' @param convenience_sample Data frame containing the non-probabilistic sample.
#' @param reference_sample Data frame containing the probabilistic sample.
#' @param covariates String vector specifying the common variables to use for training.
#' @param algorithm A string specifying which classification or regression model to use (same as caret's method).
#' @param smooth A logical value; if TRUE, propensity estimates pi_i are smoothed applying the formula (1000*pi_i + 0.5)/1001
#' @param proc A string or vector of strings specifying if any of the data preprocessing techniques available in \link[caret]{train} function from `caret` package should be applied to data prior to the propensity estimation. By default, its value is NULL and no preprocessing is applied.
#' @param trControl A trainControl specifying the computational nuances of the \link[caret]{train} function.
#' @param ... Further parameters to be passed to the \link[caret]{train} function.
#' @return A list containing `convenience` propensities and `reference` propensities.
#' @references Buskirk, T. D., & Kolenikov, S. (2015). \emph{Finding respondents in the forest: A comparison of logistic regression and random forest models for response propensity weighting and stratification.} Survey Methods: Insights from the Field, 17.
#' @examples
#' #Simple example with default parameters
#' covariates = c("education_primaria", "education_secundaria", "education_terciaria")
#' propensities(sampleNP, sampleP, covariates)
propensities = function(convenience_sample, reference_sample, covariates, algorithm = "glm", smooth = FALSE, proc = NULL, trControl = trainControl(classProbs = TRUE), ...) {
	n_convenience = nrow(convenience_sample)
	n_reference = nrow(reference_sample)

	data = rbind(convenience_sample[, covariates, drop = FALSE], reference_sample[, covariates, drop = FALSE])
	labels = append(rep(1, n_convenience), rep(0, n_reference))
	model_weights = append(rep(1, n_convenience), rep(n_convenience / n_reference, n_reference))
	
	trControl$classProbs = TRUE
	model = train(data, factor(labels, levels = c(1, 0), labels = c("Positive", "Negative")), algorithm, weights = model_weights,
		preProcess = proc, trControl = trControl, ...)
	probabilities = predict(model, data, type = "prob")$Positive

	if (smooth)
		probabilities = (1000 * probabilities + 0.5) / 1001

	list(
		convenience = probabilities[1:n_convenience],
		reference = probabilities[(n_convenience + 1):length(probabilities)]
	)
}

#' Predicts unknown responses
#'
#' It uses the matching method introduced by Rivers (2007). The idea is to model the relationship between y_k and x_k using the convenience sample in order to predict y_k for the reference sample. You can then predict the total using the `total_estimation` method.
#'
#' Training of the models is done via the `caret` package. The algorithm specified in \code{algorithm} must match one of the names in the list of algorithms supported by `caret`. If the estimated variable is categorical, probabilities are returned.
#'
#' @param convenience_sample Data frame containing the non-probabilistic sample.
#' @param reference_sample Data frame containing the probabilistic sample.
#' @param covariates String vector specifying the common variables to use for training.
#' @param estimated_var String specifying the variable to estimate.
#' @param positive_label String specifying the label to be considered positive if the estimated variable is categorical. Leave it as the default NULL otherwise.
#' @param algorithm A string specifying which classification or regression model to use (same as caret's method).
#' @param proc A string or vector of strings specifying if any of the data preprocessing techniques available in \link[caret]{train} function from `caret` package should be applied to data prior to the propensity estimation. By default, its value is NULL and no preprocessing is applied.
#' @param ... Further parameters to be passed to the \link[caret]{train} function.
#' @return A vector containing the estimated responses for the reference sample.
#' @references Rivers, D. (2007). \emph{Sampling for Web Surveys}. Presented in Joint Statistical Meetings, Salt Lake City, UT.
#' @examples
#' #Simple example with default parameters
#' N = 50000
#' covariates = c("education_primaria", "education_secundaria", "education_terciaria")
#' estimated_votes = data.frame(
#'    vote_gen = matching(sampleNP, sampleP, covariates, "vote_gen")
#' )
#' total_estimation(estimated_votes, N / nrow(estimated_votes), c("vote_gen"), N)
matching = function(convenience_sample, reference_sample, covariates, estimated_var, positive_label = NULL, algorithm = "glm", proc = NULL, ...) {
	data = convenience_sample[, covariates, drop = FALSE]
	values = convenience_sample[, estimated_var]
	test = reference_sample[, covariates, drop = FALSE]
	model = train(data, values, algorithm, preProcess = proc, trControl = trainControl(classProbs =  TRUE), ...)

	if (is.null(positive_label))
		return(predict(model, test))
	else
		return(predict(model, test, type = "prob")[, positive_label])
}

#' Calculates a model based estimation
#'
#' It uses the model based estimator. The idea in order to estimate the population total is to add the sample responses and the predicted responses for the individuals not contained in the sample. See for example Valliant et al. (2000).
#'
#' Training of the models is done via the `caret` package. The algorithm specified in \code{algorithm} must match one of the names in the list of algorithms supported by `caret`.
#'
#' @param sample_data Data frame containing the sample.
#' @param full_data Data frame containing all the individuals contained in the population.
#' @param covariates String vector specifying the common variables to use for training.
#' @param estimated_var String specifying the variable to estimate.
#' @param estimate_mean Boolean specifying whether the mean estimation should be returned. Otherwise, the total estimation is returned by default.
#' @param positive_label String specifying the label to be considered positive if the estimated variable is categorical. Leave it as the default NULL otherwise.
#' @param algorithm A string specifying which classification or regression model to use (same as caret's method).
#' @param proc A string or vector of strings specifying if any of the data preprocessing techniques available in \link[caret]{train} function from `caret` package should be applied to data prior to the propensity estimation. By default, its value is NULL and no preprocessing is applied.
#' @param ... Further parameters to be passed to the \link[caret]{train} function.
#' @return The population total estimation (or mean if specified by the `estimate_mean` parameter).
#' @references Valliant, R., Dorfman, A. H., & Royall, R. M. (2000) \emph{Finite population sampling and inference: a prediction approach.} Wiley, New York.
#' @examples
#' #Simple example with default parameters
#' covariates = c("education_primaria", "education_secundaria",
#'    "education_terciaria", "age", "sex", "language")
#' model_based(sampleNP, population, covariates, "vote_gen")
model_based = function(sample_data, full_data, covariates, estimated_var, estimate_mean = FALSE, positive_label = NULL, algorithm = "glm", proc = NULL, ...) {
	known_values = sample_data[, estimated_var]
	if (!is.null(positive_label))
		known_values = known_values == positive_label
	all_data = rbind(sample_data[, covariates, drop = FALSE], full_data[, covariates, drop = FALSE])
	all_predicted_values = matching(sample_data, all_data, covariates, estimated_var, positive_label, algorithm = algorithm, proc = proc, ...)
	known_predicted_values = all_predicted_values[1:length(known_values)]
	predicted_values = all_predicted_values[(length(known_values) + 1):length(all_predicted_values)]
	total = sum(known_values, predicted_values, -known_predicted_values)
	
	if (estimate_mean)
		return(total / nrow(full_data))
	else
		return(total)
}

#' Calculates a model assisted estimation
#'
#' It uses the model assisted estimator introduced by Särndal et al. (1992).
#'
#' Training of the models is done via the `caret` package. The algorithm specified in \code{algorithm} must match one of the names in the list of algorithms supported by `caret`.
#'
#' @param sample_data Data frame containing the sample.
#' @param weights Vector containing the sample weights.
#' @param full_data Data frame containing all the individuals contained in the population.
#' @param covariates String vector specifying the common variables to use for training.
#' @param estimated_var String specifying the variable to estimate.
#' @param estimate_mean Boolean specifying whether the mean estimation should be returned. Otherwise, the total estimation is returned by default.
#' @param positive_label String specifying the label to be considered positive if the estimated variable is categorical. Leave it as the default NULL otherwise.
#' @param algorithm A string specifying which classification or regression model to use (same as caret's method).
#' @param proc A string or vector of strings specifying if any of the data preprocessing techniques available in \link[caret]{train} function from `caret` package should be applied to data prior to the propensity estimation. By default, its value is NULL and no preprocessing is applied.
#' @param ... Further parameters to be passed to the \link[caret]{train} function.
#' @return The population total estimation (or mean if specified by the `estimate_mean` parameter).
#' @references Särndal, C. E., Swensson, B., & Wretman, J. (1992). \emph{Model assisted survey sampling.} Springer, New York.
#' @examples
#' #Simple example with default parameters
#' covariates = c("education_primaria", "education_secundaria",
#'    "education_terciaria", "age", "sex", "language")
#' model_assisted(sampleNP, nrow(population) / nrow(sampleNP),
#'    population, covariates, "vote_gen")
model_assisted = function(sample_data, weights, full_data, covariates, estimated_var, estimate_mean = FALSE, positive_label = NULL, algorithm = "glm", proc = NULL, ...) {
	known_values = sample_data[, estimated_var]
	if (!is.null(positive_label))
		known_values = known_values == positive_label
	all_data = rbind(sample_data[, covariates, drop = FALSE], full_data[, covariates, drop = FALSE])
	predicted_values = matching(sample_data, all_data, covariates, estimated_var, positive_label, algorithm = algorithm, proc = proc, ...)
	known_predicted_values = predicted_values[1:length(known_values)]
	total = sum(predicted_values, -known_predicted_values) + sum((known_values - known_predicted_values) * weights)
	
	if (estimate_mean)
		return(total / nrow(full_data))
	else
		return(total)
}

#' Calculates Valliant weights
#'
#' Computes weights from propensity estimates using the 1/pi_i formula introduced in Valliant (2019).
#'
#' The function takes the vector of propensities \eqn{\pi(x)} and calculates the weights to be applied in the Hajek estimator using the formula that can be found in Valliant (2019). For an individual \emph{i}, weight is calculated as follows:
#' \deqn{w_i = 1/\pi_i (x)}
#'
#' @param propensities A vector with the propensities associated to the elements of the convenience sample.
#' @return A vector with the corresponding weights.
#' @references Valliant, R. (2019). \emph{Comparing Alternatives for Estimation from Nonprobability Samples}. Journal of Survey Statistics and Methodology, smz003, \url{https://doi.org/10.1093/jssam/smz003}
#' @examples
#' covariates = c("education_primaria", "education_secundaria", "education_terciaria")
#' data_propensities = propensities(sampleNP, sampleP, covariates)
#' valliant_weights(data_propensities$convenience)
valliant_weights = function(propensities) {
	1 / propensities
}

#' Calculates Schonlau and Couper weights
#'
#' Computes weights from propensity estimates using the (1 - pi_i)/pi_i formula introduced in Schonlau and Couper (2017).
#'
#' The function takes the vector of propensities \eqn{\pi(x)} and calculates the weights to be applied in the Hajek estimator using the formula that can be found in Schonlau and Couper (2017). For an individual \emph{i}, weight is calculated as follows:
#' \deqn{w_i = \frac{1 - \pi_i (x)}{\pi_i (x)}}
#'
#' @param propensities A vector with the propensities associated to the elements of the convenience sample.
#' @return A vector with the corresponding weights.
#' @references Schonlau, M., & Couper, M. P. (2017). \emph{Options for conducting web surveys.} Statistical Science, 32(2), 279-292.
#' @examples
#' covariates = c("education_primaria", "education_secundaria", "education_terciaria")
#' data_propensities = propensities(sampleNP, sampleP, covariates)
#' sc_weights(data_propensities$convenience)
sc_weights = function(propensities) {
	(1 - propensities) / propensities
}

#' Calculates Valliant and Dever weights
#'
#' Computes weights from propensity estimates using the propensity stratification 1/p_i averaging formula introduced in Valliant and Dever (2011).
#'
#' The function takes the vector of propensities \eqn{\pi(x)} and calculates the weights to be applied in the Horvitz-Thompson estimator using the formula that can be found in Valliant and Dever (2019). The vector of propensities is divided in \emph{g} strata (ideally five according to Cochran, 1968) aiming to have individuals with similar propensities in each strata. After the stratification, weight is calculated as follows for an individual \emph{i}:
#' \deqn{w_i = \frac{n(g_i)}{ \sum_{k \in g_i} \pi_k (x)}}
#' where \eqn{g_i} represents the strata to which \emph{i} belongs, and \eqn{n(g_i)} is the number of individuals in the \eqn{g_i} strata.
#'
#' @param convenience_propensities A vector with the propensities associated with the convenience sample.
#' @param reference_propensities A vector with the propensities associated with the reference sample.
#' @param g The number of strata to use; by default, its value is 5.
#' @return A vector with the corresponding weights.
#' @references Valliant, R., & Dever, J. A. (2011). \emph{Estimating propensity adjustments for volunteer web surveys.} Sociological Methods & Research, 40(1), 105-137.
#' @references Cochran, W. G. (1968). \emph{The Effectiveness of Adjustment by Subclassification in Removing Bias in Observational Studies.} Biometrics, 24(2), 295-313
#' @examples
#' covariates = c("education_primaria", "education_secundaria", "education_terciaria")
#' data_propensities = propensities(sampleNP, sampleP, covariates)
#' vd_weights(data_propensities$convenience, data_propensities$reference)
vd_weights = function(convenience_propensities, reference_propensities, g = 5) {
	propensities = append(convenience_propensities, reference_propensities)
	
	cuts = cut(1:length(propensities), g, labels = FALSE)
	data_order = order(propensities)
	strata = vector()
	
	for (i in 1:g) {
		strata[data_order[cuts == i]] = i
	}

	stratum_weights = sapply(1:g, function(i) {
		1 / mean(propensities[strata == i])
	})

	stratum_weights[strata[1:length(convenience_propensities)]]
}

#' Calculates Lee weights
#'
#' Computes weights from propensity estimates using the propensity stratification design weights averaging formula introduced in Lee (2006) and Lee and Valliant (2009).
#'
#' The function takes the vector of propensities \eqn{\pi(x)} and calculates the weights to be applied in the Horvitz-Thompson estimator using the formula that can be found in Lee (2006) and Lee and Valliant (2009). The vector of propensities is divided in \emph{g} strata (ideally five according to Cochran, 1968) aiming to have individuals with similar propensities in each strata. After the stratification, weight is calculated as follows for an individual \emph{i}:
#' \deqn{w_i = \frac{n_r(g_i) / n_r}{n_v(g_i) / n_v}}
#' where \eqn{g_i} represents the strata to which \emph{i} belongs, \eqn{n_r (g_i)} and \eqn{n_v (g_i)} are the number of individuals in the \eqn{g_i} strata from the reference and the convenience sample respectively, and \eqn{n_r} and \eqn{n_v} are the sample sizes for the reference and the convenience sample respectively.
#'
#' @param convenience_propensities A vector with the propensities associated with the convenience sample.
#' @param reference_propensities A vector with the propensities associated with the reference sample.
#' @param g The number of strata to use; by default, its value is 5.
#' @return A vector with the corresponding weights.
#' @references Lee, S. (2006). \emph{Propensity score adjustment as a weighting scheme for volunteer panel web surveys.} Journal of official statistics, 22(2), 329.
#' @references Lee, S., & Valliant, R. (2009). \emph{Estimation for volunteer panel web surveys using propensity score adjustment and calibration adjustment.} Sociological Methods & Research, 37(3), 319-343.
#' @references Cochran, W. G. (1968). \emph{The Effectiveness of Adjustment by Subclassification in Removing Bias in Observational Studies.} Biometrics, 24(2), 295-313
#' @examples
#' covariates = c("education_primaria", "education_secundaria", "education_terciaria")
#' data_propensities = propensities(sampleNP, sampleP, covariates)
#' lee_weights(data_propensities$convenience, data_propensities$reference)
lee_weights = function(convenience_propensities, reference_propensities, g = 5) {
	propensities = append(convenience_propensities, reference_propensities)
	
	cuts = cut(1:length(propensities), g, labels = FALSE)
	data_order = order(propensities)
	strata = vector()
	
	for (i in 1:g) {
		strata[data_order[cuts == i]] = i
	}

	n_convenience = length(convenience_propensities)
	n_reference = length(reference_propensities)
	convenience_strata = strata[1:n_convenience]
	refence_strata = strata[(n_convenience + 1):length(strata)]

	stratum_weights = sapply(1:g, function(i) {
		(sum(refence_strata == i) / n_reference) / (sum(convenience_strata == i) / n_convenience)
	})

	stratum_weights[convenience_strata]
}

#' Weights of the calibration estimator
#'
#' Calculates the calibration weights from a disjunct matrix of covariates, a vector of population totals and a vector of initial weights.
#'
#' The function uses the `calib` function from the `sampling` package for the estimation of g-weights, which are multiplied by the initial weights to obtain the final calibration weights. The initial weights can be calculated previously from the propensities for any of the implemented methods (see functions \code{lee_weights}, \code{sc_weights}, \code{valliant_weights}, \code{vd_weights}). The population size is used to scale said initial weights so they are easier to calibrate.
#'
#' @param Xs Matrix of calibration variables.
#' @param totals A vector containing population totals for each column (class) of the calibration variables matrix.
#' @param initial_weights A vector containing the initial weights for each individual.
#' @param N Integer indicating the population size.
#' @param ... Further arguments to be passed to the `calib` function from the `sampling` package.
#' @return A vector with the corresponding weights.
#' @examples
#' n = nrow(sampleNP)
#' N = 50000
#' language_total = 45429
#' covariates = c("education_primaria","education_secundaria",
#'    "education_terciaria","age","sex")
#' pi = propensities(sampleNP, sampleP, covariates, algorithm = "glm", smooth = FALSE)
#' wi = sc_weights(pi$convenience)
#' calib_weights(sampleNP$language, language_total, wi, N, method = "raking")
calib_weights = function(Xs, totals, initial_weights, N, ...) {
	initial_weights = initial_weights * (N / sum(initial_weights))
	calib(Xs, initial_weights, totals, ...) * initial_weights
}

#' Estimates the population means
#'
#' Estimates the means for the specified variables measured in a sample given some pre-calculated weights.
#'
#' @param sample A data frame containing the sample with the variables for which the means are to be calculated.
#' @param weights A vector of pre-calculated weights.
#' @param estimated_vars String vector specifying the variables in the sample to be estimated.
#' @param N An integer specifying the population size (optional).
#' @return A vector with the corresponding estimations.
#' @examples
#' covariates = c("education_primaria", "education_secundaria", "education_terciaria")
#' data_propensities = propensities(sampleNP, sampleP, covariates)
#' psa_weights = sc_weights(data_propensities$convenience)
#' mean_estimation(sampleNP, psa_weights, c("vote_pens"))
mean_estimation = function(sample, weights, estimated_vars, N = NULL) {
	total = ifelse(is.null(N), sum(weights), N)

	sapply(estimated_vars, function(var_name) {
		sum(sample[, var_name] * weights) / total
	})
}

#' Estimates the population totals
#'
#' Estimates the population totals for the specified variables measured in a sample given some pre-calculated weights.
#'
#' @param sample A data frame containing the sample with the variables for which the estimated population totals are to be calculated.
#' @param weights A vector of pre-calculated weights.
#' @param estimated_vars String vector specifying the variables in the sample to be estimated.
#' @param N An integer specifying the population size.
#' @return A vector with the corresponding estimations.
#' @examples
#' covariates = c("education_primaria", "education_secundaria", "education_terciaria")
#' data_propensities = propensities(sampleNP, sampleP, covariates)
#' psa_weights = sc_weights(data_propensities$convenience)
#' total_estimation(sampleNP, psa_weights, c("vote_pens"), 50000)
total_estimation = function(sample, weights, estimated_vars, N) {
	sum_weights = sum(weights)

	sapply(estimated_vars, function(var_name) {
		sum(sample[, var_name] * weights) / sum_weights * N
	})
}

#' Estimates the population proportion
#'
#' Estimates the proportion of a given class or classes for the specified variables measured in a sample given some pre-calculated weights.
#'
#' @param sample A data frame containing the sample with the variables for which the means are to be calculated.
#' @param weights A vector of pre-calculated weights.
#' @param estimated_vars String vector specifying the variables in the sample to be estimated.
#' @param class String vector specifying which class (value) proportion is to be estimated in each variable. The \emph{i}-th element of this vector corresponds to the class of which proportion is desired to estimate of the \emph{i}-th variable of the vector specified in \code{estimated_vars}.
#' @param N An integer specifying the population size (optional).
#' @return A vector with the corresponding estimations.
#' @examples
#' covariates = c("education_primaria", "education_secundaria", "education_terciaria")
#' data_propensities = propensities(sampleNP, sampleP, covariates)
#' psa_weights = sc_weights(data_propensities$convenience)
#'
#' #The function will estimate the proportion of individuals
#' #with the 0 value in vote_pens and the 1 value in vote_gen
#' prop_estimation(sampleNP, psa_weights, c("vote_pens", "vote_gen"), c(0, 1))
prop_estimation = function(sample, weights, estimated_vars, class, N = NULL) {
  total = ifelse(is.null(N), sum(weights), N)

  sapply(estimated_vars, function(var_name) {
    sum(sample[, var_name] == class[estimated_vars %in% var_name] * weights) / total
  })
}

#' Calculates Jackknife variance with reweighting
#'
#' Calculates the variance of a given estimator by Leave-One-Out Jackknife (Quenouille, 1956) with reweighting in each iteration.
#'
#' The estimation of the variance requires a recalculation of the estimates in each iteration which might involve weighting adjustments, leading to an increase in computation time. It is expected that the estimated variance captures the weighting adjustments' variability and the estimator's variability.
#'
#' @param estimated_vars A string vector specifying the variables for which the estimators' variance are to be estimated.
#' @param convenience_sample Data frame containing the non-probabilistic sample.
#' @param reference_sample Data frame containing the probabilistic sample.
#' @param covariates String vector specifying the common variables to use for training.
#' @param N Integer indicating the population size. Optional.
#' @param algorithm A string specifying which classification or regression model to use (same as caret's method). By default, its value is "glm" (logistic regression).
#' @param smooth A logical value; if TRUE, propensity estimates pi_i are smoothed applying the formula (1000*pi_i + 0.5)/1001
#' @param proc A string or vector of strings specifying if any of the data preprocessing techniques available in \link[caret]{train} function from `caret` package should be applied to data prior to the propensity estimation. By default, its value is NULL and no preprocessing is applied.
#' @param trControl A trainControl specifying the computational nuances of the \link[caret]{train} function.
#' @param weighting.func A string specifying which function should be used to compute weights from propensity scores. Available functions are the following: \itemize{ \item \code{sc} calls \link{sc_weights}. \item \code{valliant} calls \link{valliant_weights}. \item \code{lee} calls \link{lee_weights}. \item \code{vd} calls \link{vd_weights}. }
#' @param g If \code{weighting.func = "lee"} or \code{weighting.func = "vd"}, this element specifies the number of strata to use; by default, its value is 5.
#' @param calib A logical value; if TRUE, PSA weights are used as initial weights for calibration. By default, its value is FALSE.
#' @param calib_vars A string or vector of strings specifying the variables to be used for calibration. By default, its value is NULL.
#' @param totals A vector containing population totals for each column (class) of the calibration variables matrix. Ignored if \code{calib} is set to FALSE.
#' @param args.calib A list containing further arguments to be passed to the \link{calib_weights} function.
#' @param ... Further parameters to be passed to the \link[caret]{train} function.
#' @return The resulting variance.
#' @references Quenouille, M. H. (1956). \emph{Notes on bias in estimation.} Biometrika, 43(3/4), 353-360.
#' @examples
#' \donttest{
#' #A simple example without calibration and default parameters
#' covariates = c("education_primaria", "education_secundaria", "education_terciaria")
#' jackknife_variance("vote_pens",sampleNP, sampleP, covariates)
#'
#' #An example with linear calibration and default parameters
#' covariates = c("education_primaria", "education_secundaria", "education_terciaria")
#' calib_vars = c("age", "sex")
#' totals = c(2544377, 24284)
#'
#' jackknife_variance("vote_pens",sampleNP, sampleP, covariates,
#' calib = T, calib_vars, totals, args.calib = list(method = "linear"))
#' }
jackknife_variance = function(estimated_vars, convenience_sample, reference_sample,
                              covariates, N = NULL, algorithm = "glm", smooth = FALSE, proc = NULL, trControl = trainControl(classProbs = TRUE), weighting.func = "sc", g = 5,
                              calib = FALSE, calib_vars = NULL, totals = NULL, args.calib = NULL, ...) {
	sample_size = nrow(convenience_sample)
	correction_factor = ifelse(is.null(N), 1, 1 - sample_size / N)

	estimations = sapply(1:sample_size, function(i) {
		sample <- convenience_sample[-i,]
		pi <- propensities(sample, reference_sample, covariates, algorithm, smooth, proc, trControl, ...)
		if(weighting.func == "sc") wi <- sc_weights(pi$convenience)
		if(weighting.func == "valliant") wi <- valliant_weights(pi$convenience)
		if(weighting.func == "lee") wi <- lee_weights(pi$convenience, pi$reference, g)
		if(weighting.func == "vd") wi <- vd_weights(pi$convenience, pi$reference, g)
		if(calib) wi <- do.call(calib_weights, append(list(Xs = sample[, calib_vars],
	                                           totals = totals, initial_weights = wi), args.calib))
		sum(sample[, estimated_vars] * wi) / sum(wi)
	})

	(sample_size - 1) / sample_size * sum((estimations - mean(estimations))^2) * correction_factor
}

#' Calculates Jackknife variance with reweighting
#'
#' Calculates the variance of a given estimator by Leave-One-Out Jackknife (Quenouille, 1956) with reweighting in each iteration.
#'
#' The estimation of the variance requires a recalculation of the estimates in each iteration which might involve weighting adjustments, leading to an increase in computation time. It is expected that the estimated variance captures the weighting adjustments' variability and the estimator's variability.
#'
#' @param sample Data frame containing the non-probabilistic sample.
#' @param estimator Function that, given a sample as a parameter, returns an estimation.
#' @param N Integer indicating the population size. Optional.
#' @return The resulting variance.
#' @references Quenouille, M. H. (1956). \emph{Notes on bias in estimation.} Biometrika, 43(3/4), 353-360.
#' @examples
#' \donttest{
#' covariates = c("education_primaria", "education_secundaria",
#'    "education_terciaria", "age", "sex", "language")
#' vote_gen_estimator = function(sample) {
#'    model_based(sample, population, covariates, "vote_gen")
#' }
#' generic_jackknife_variance(sampleNP, vote_gen_estimator)
#' }
generic_jackknife_variance = function(sample, estimator, N = NULL) {
	sample_size = nrow(sample)
	correction_factor = ifelse(is.null(N), 1, 1 - sample_size / N)
	
	estimations = sapply(1:sample_size, function(i) {
		estimator(sample[-i,])
	})

	(sample_size - 1) / sample_size * sum((estimations - mean(estimations))^2) * correction_factor
}

#' Calculates Jackknife variance without reweighting
#'
#' Calculates the variance of a given estimator by Leave-One-Out Jackknife (Quenouille, 1956) with the original adjusted weights.
#'
#' The variance estimation is performed by eliminating an individual at each iteration with its corresponding weight and estimating the mean of the corresponding subsample, which is further used in the Jackknife formula as the usual procedure. The calculation of variance estimates through this procedure might take less computation time but also might not take into account the variance of the weighting method.
#'
#' @param sample A data frame containing the sample.
#' @param weights A vector containing the pre-calculated weights.
#' @param estimated_vars A string vector specifying the variables for which the estimators' variance are to be estimated.
#' @param N Integer indicating the population size. Optional.
#' @return A vector containing the resulting variance for each variable.
#' @references Quenouille, M. H. (1956). \emph{Notes on bias in estimation.} Biometrika, 43(3/4), 353-360.
#' @examples
#' covariates = c("education_primaria", "education_secundaria", "education_terciaria")
#' data_propensities = propensities(sampleNP, sampleP, covariates)
#' psa_weights = sc_weights(data_propensities$convenience)
#' fast_jackknife_variance(sampleNP, psa_weights, c("vote_pens"), 50000)
fast_jackknife_variance = function(sample, weights, estimated_vars, N = NULL) {
	sample_size = nrow(sample)
	correction_factor = ifelse(is.null(N), 1, 1 - sample_size / N)

	sapply(estimated_vars, function(var_name) {
		estimations = sapply(1:sample_size, function(i) {
			sum(sample[-i, var_name] * weights[-i]) / sum(weights[-i])
		})

		(sample_size - 1) / sample_size * sum((estimations - mean(estimations))^2) * correction_factor
	})
}

#' Confidence interval
#'
#' Calculates the confidence interval for the estimator considered.
#'
#' @param estimation A numeric value specifying the point estimation.
#' @param std_dev A numeric value specifying the standard deviation of the point estimation.
#' @param confidence A numeric value between 0 and 1 specifying the confidence level, taken as 1 - alpha (1 - Type I error). By default, its value is 0.95.
#' @return A vector containing the lower and upper bounds.
#' @examples
#' covariates = c("education_primaria","education_secundaria",
#' "education_terciaria", "age", "sex")
#' pi = propensities(sampleNP, sampleP, covariates, algorithm = "glm", smooth = FALSE)
#' psa_weights = sc_weights(pi$convenience)
#' N = 50000
#' Y_est = total_estimation(sampleNP, psa_weights, estimated_vars = "vote_pens", N = N)
#' VY_est = fast_jackknife_variance(sampleNP, psa_weights,
#'    estimated_vars = "vote_pens") * N^2
#' confidence_interval(Y_est, sqrt(VY_est), confidence = 0.90)
confidence_interval = function(estimation, std_dev, confidence = 0.95) {
	deviation = std_dev * qnorm(1 - (1 - confidence) / 2)
	c(lower = estimation - deviation, upper = estimation + deviation)
}

#' A probabilistic sample
#'
#' A dataset of 500 individuals extracted with simple random sampling from a simulated fictitious population of 50,000 individuals. Further details on the generation of the dataset can be found in Ferri-García and Rueda (2018). The variables present in the dataset are the following:
#' \itemize{
#'   \item education_primaria. A binary variable indicating if the highest academic level achieved by the individual is Primary Education.
#'   \item education_secundaria. A binary variable indicating if the highest academic level achieved by the individual is Secondary Education.
#'   \item education_terciaria. A binary variable indicating if the highest academic level achieved by the individual is Tertiary Education.
#'   \item age. A numeric variable, with values ranging from 18 to 100, indicating the age of the individual.
#'   \item sex. A binary variable indicating if the individual is a man.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name sampleP
#' @usage sampleP
#' @references Ferri-García, R., & Rueda, M. (2018). \emph{Efficiency of propensity score adjustment and calibration on the estimation from non-probabilistic online surveys}. SORT-Statistics and Operations Research Transactions, 1(2), 159-162.
"sampleP"

#' A non-probabilistic sample
#'
#' A dataset of 1000 individuals extracted from the subpopulation of individuals with internet access in a simulated fictitious population of 50,000 individuals. This sample attempts to reproduce a case of nonprobability sampling with selection bias, as there are important differences between the potentially covered population, the covered population and the full target population. Further details on the generation of the dataset can be found in Ferri-García and Rueda (2018). The variables present in the dataset are the following:
#' \itemize{
#'   \item vote_gen. A binary variable indicating if the individual vote preferences are for Party 1. This variable is related to gender.
#'   \item vote_pens. A binary variable indicating if the individual vote preferences are for Party 2. This variable is related to age.
#'   \item vote_pir. A binary variable indicating if the individual vote preferences are for Party 3. This variable is related to age and internet access.
#'   \item education_primaria. A binary variable indicating if the highest academic level achieved by the individual is Primary Education.
#'   \item education_secundaria. A binary variable indicating if the highest academic level achieved by the individual is Secondary Education.
#'   \item education_terciaria. A binary variable indicating if the highest academic level achieved by the individual is Tertiary Education.
#'   \item age. A numeric variable, with values ranging from 18 to 100, indicating the age of the individual.
#'   \item sex. A binary variable indicating if the individual is a man.
#'   \item language. A binary variable indicating if the individual is a native.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name sampleNP
#' @usage sampleNP
#' @references Ferri-García, R., & Rueda, M. (2018). \emph{Efficiency of propensity score adjustment and calibration on the estimation from non-probabilistic online surveys}. SORT-Statistics and Operations Research Transactions, 1(2), 159-162.
"sampleNP"

#' A full population
#'
#' A dataset of a simulated fictitious population of 50,000 individuals. Further details on the generation of the dataset can be found in Ferri-García and Rueda (2018). The variables present in the dataset are the following:
#' \itemize{
#'   \item education_primaria. A binary variable indicating if the highest academic level achieved by the individual is Primary Education.
#'   \item education_secundaria. A binary variable indicating if the highest academic level achieved by the individual is Secondary Education.
#'   \item education_terciaria. A binary variable indicating if the highest academic level achieved by the individual is Tertiary Education.
#'   \item age. A numeric variable, with values ranging from 18 to 100, indicating the age of the individual.
#'   \item sex. A binary variable indicating if the individual is a man.
#'   \item language. A binary variable indicating if the individual is a native.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name population
#' @usage population
#' @references Ferri-García, R., & Rueda, M. (2018). \emph{Efficiency of propensity score adjustment and calibration on the estimation from non-probabilistic online surveys}. SORT-Statistics and Operations Research Transactions, 1(2), 159-162.
"population"
