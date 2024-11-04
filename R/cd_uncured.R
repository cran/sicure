#' Cross-validation conditional density of the susceptible population
#'
#' This function implements a nonparametric cross-validation estimator for the conditional density  of the susceptible population, as proposed by Piñeiro-Lamas (2024) (see Equation (3.5)).
#' A leave-one-out cross-validation approach is considered to ensure that the sample used to construct the estimator and the point at which it is evaluated
#' are independent.
#'
#'
#' @param x A numeric vector giving the covariate values.
#' @param time A numeric vector giving the observed times.
#' @param delta A numeric vector giving the values of the uncensoring indicator, where 1 indicates that the event of interest has been observed and 0 indicates that the observation is censored.
#' @param logh3 The logarithm of the bandwidth for smoothing the time variable.
#' @param logh4 The logarithm of the bandwidth for smoothing the covariate.
#'
#' @return A vector containing the cross-validation conditional density of the
#' susceptible population for each observation \eqn{(X_i, T_i)}.
#'
#' @references
#' Piñeiro-Lamas, B. (2024). High dimensional single-index mixture cure models [PhD thesis]. Universidade da Coruña. Available at \url{https://ruc.udc.es/dspace/handle/2183/37035}
#'
#' @export
#'
#' @seealso \code{\link[sicure]{cd.sm.uncured}}
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
#' suppressWarnings(cd.uncured(x, t, d, log(0.5), log(0.5)))
cd.uncured <- function(x, time, delta, logh3, logh4){
  n = length(x)
  matrix.difS0 = matrix(0, nrow=n-1, ncol=n)
  for(i in 1:n){
    dat.minus.i <- data.frame(x, time, delta)[-i,]
    S0.minus.i <- npcure::latency(x, time, delta, dat.minus.i, x0=x[i],
                          h=exp(logh4), local=FALSE, testimate=time)$S[[1]][[1]]
    S0.minus.i <- c(1, S0.minus.i)
    S0.minus.i.dif <- -diff(S0.minus.i)[-i]
    matrix.difS0[,i] <- S0.minus.i.dif
  }

  matTime <- matrix(rep(time, n), byrow=TRUE, ncol=n)
  matTimedif <- matTime-t(matTime)
  matrix_epa <- k.epa(exp(logh3),matTimedif)
  diag(matrix_epa) <- NA
  matrix_epa <- matrix(matrix_epa[which(!is.na(matrix_epa))], nrow=n-1, ncol=n)
  conditional_density_CV_vector <- colSums(matrix.difS0*matrix_epa)
  return(conditional_density_CV_vector)
}
