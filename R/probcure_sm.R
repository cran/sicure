#' Smoothed version of the nonparametric estimator of the conditional
#' probability of cure
#'
#' This function computes a smoothed version of the nonparametric estimator
#' of the probability of cure proposed by Xu and Peng (2014) and deeply studied
#' by López-Cheda et al. (2017). The smoothing is performed using the \code{\link[caTools]{runmean}}
#' function, which computes a moving average of the estimated probabilities in a window determined
#' by a radius \code{r}.
#' The non-smoothed version is implemented in the \code{\link[npcure]{probcure}} function of the \code{npcure} package (López-Cheda et al., 2021).
#'
#' @param x A numeric vector giving the covariate values.
#' @param time A numeric vector giving the observed times.
#' @param delta A numeric vector giving the values of the uncensoring indicator, where 1 indicates that the event of interest has been observed and 0 indicates that the observation is censored.
#' @param logh1 The logarithm of the bandwidth for smoothing the covariate.
#' @param r Radius of moving window.
#'
#'
#' @return A list with two components:
#' \itemize{
#' \item A vector containing the cross-validation estimations of the probability of cure.
#' \item The previous vector smoothed with \code{\link[caTools]{runmean}} with a moving window of \eqn{k = 2r + 1}.
#' }
#'
#' @references
#' López-Cheda, A., Cao, R., Jácome, M. A., Van Keilegom, I. (2017). Nonparametric incidence estimation and bootstrap bandwidth selection in mixture cure models. Computational Statistics & Data Analysis, 105, 144–165. \doi{10.1016/j.csda.2016.08.002}.
#'
#' López-Cheda, A., Jácome, M. A., López-de-Ullibarri, I. (2021). The R Journal, 13(1), 21-41. \doi{10.32614/RJ-2021-027}.
#'
#' Xu, J., Peng, Y. (2014). Nonparametric cure rate estimation with covariates. The Canadian Journal of Statistics, 42, 1-17. \doi{10.1002/cjs.11197}.
#'
#'
#' @export
#'
#' @seealso \code{\link[npcure]{probcure}}
#'
#' @examples
#' # Some artificial data
#' set.seed(123)
#' n <- 50
#' x <- runif(n, -2, 2) # Covariate values
#' y <- rweibull(n, shape = 0.5 * (x + 4)) # True lifetimes
#' c <- rexp(n) # Censoring values
#' p <- exp(2*x)/(1 + exp(2*x)) # Probability of being susceptible
#' u <- runif(n)
#' t  <- ifelse(u < p, pmin(y, c), c) # Observed times
#' d  <- ifelse(u < p, ifelse(y < c, 1, 0), 0) # Uncensoring indicator
#' data <- data.frame(x = x, t = t, d = d)
#'
#' # Smoothed nonparametric estimates of cure probability with bandwidth=2
#' q1 <- probcure.sm(x, t, d, logh1 = log(2), r=2)[[2]]
#' plot(sort(x), q1[order(x)], type = "l", xlab = "Covariate", ylab = "Cure probability",
#'      ylim = c(0, 1))

probcure.sm <- function(x, time, delta, logh1, r){
  n <- length(x)
  probcure_CV <- numeric(n)   # Initialize a vector to store cross-validated cure probability estimates
  for(i in 1:n){ # Estimate the cure probability using leave-one-out cross-validation
    probcure_CV[i] <- npcure::probcure(x, time, delta, data.frame(x,time, delta)[-i,],
                              x0=x[i], local=FALSE, h=exp(logh1))$q[[1]]
  }

  # Smooth, with runmean, the CV cure probability estimates
  probcure_CV_ord <- probcure_CV[order(data.frame(x,time,delta)[,1])]
  probcure_CV_smooth <- caTools::runmean(probcure_CV_ord, k=2*r+1)
  probcure_CV_smooth <- probcure_CV_smooth[order(order(data.frame(x,time,delta)[,1]))]
  return(list(probcure_CV, probcure_CV_smooth))
}
