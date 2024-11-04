#' K Nearest Neighbors with Mahalanobis Distance
#'
#' This function computes the \code{k} nearest neighbors for a given set of data points,
#' where each observation is a pair of the form \eqn{(X, T)}, with \eqn{X} representing a covariate and \eqn{T} the observed time.
#' The distance between each pair of points is computed using the Mahalanobis distance:
#' \deqn{ d_M((X_i, T_i), (X_j, T_j)) = \sqrt{ \left( \begin{pmatrix} X_i \\ T_i \end{pmatrix} - \begin{pmatrix} X_j \\ T_j \end{pmatrix} \right)^t \Sigma^{-1} \left( \begin{pmatrix} X_i \\ T_i \end{pmatrix} - \begin{pmatrix} X_j \\ T_j \end{pmatrix} \right) }, }
#' where \eqn{\Sigma} is the variance-covariance matrix of the joint distribution of \eqn{(X, T)}.
#'
#' @param x A numeric vector of length \eqn{n} giving the covariate values.
#' @param time A numeric vector giving the observed times.
#' @param k The number of nearest neighbors to search.
#'
#' @return A matrix with \eqn{n} rows and \code{k} columns. Each row represents
#' each pair \eqn{(X_i, T_i)}. The values in each row give the index of the
#' \code{k} nearest neighbors considering Mahalanobis distance.
#'
#' @references
#' Mahalanobis, P. C. (1936). On the generalised distance in statistics. Proceedings of the National Institute of Sciences of India, 2, 49-55.
#'
#' @export
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
#' kNN.Mahalanobis(x=x, time=t, k=5)
kNN.Mahalanobis <- function(x, time, k){
  # Compute the Mahalanobis distance matrix between all pairs of (x, time) points
  Mahalanobis_dist_matrix <- StatMatch::mahalanobis.dist(data.frame(x, time))
  # Find the indices of the k-nearest neighbors for each point based on the Mahalanobis distance
  matrix_kNN_Mahalanobis_index <- t(apply(Mahalanobis_dist_matrix, 1, doBy::which.minn, n=k))
  return(matrix_kNN_Mahalanobis_index) # Return the matrix of k-nearest neighbors indices
}


#' Smoothed cross-validation conditional density estimator of the susceptible population
#'
#' This function implements a smoothed version of the nonparametric cross-validation estimator for the conditional density  of the susceptible population proposed by Piñeiro-Lamas (2024) (see Equation (3.5)).
#' Smoothing is done using the \code{k} nearest neighbors based on Mahalanobis distance.
#' The Mahalanobis distance between two points \eqn{(X_i, T_i)} and \eqn{(X_j, T_j)} is defined as:
#' \deqn{ d_M((X_i, T_i), (X_j, T_j)) = \sqrt{ \left( \begin{pmatrix} X_i \\ T_i \end{pmatrix} - \begin{pmatrix} X_j \\ T_j \end{pmatrix} \right)^t \Sigma^{-1} \left( \begin{pmatrix} X_i \\ T_i \end{pmatrix} - \begin{pmatrix} X_j \\ T_j \end{pmatrix} \right) }, }
#' where \eqn{\Sigma} is the covariance matrix of the joint distribution of \eqn{(X, T)}.
#'
#' @param x A numeric vector giving the covariate values.
#' @param time A numeric vector giving the observed times.
#' @param delta A numeric vector giving the values of the uncensoring indicator, where 1 indicates that the event of interest has been observed and 0 indicates that the observation is censored.
#' @param logh3 The logarithm of the bandwidth for smoothing the time variable.
#' @param logh4 The logarithm of the bandwidth for smoothing the covariate.
#' @param k The number of nearest neighbors used to smooth.
#'
#' @return A vector containing the cross-validation conditional density of the susceptible population for each observation \eqn{(X_i, T_i)}, smoothed by considering the \code{k} nearest neighbors with Mahalanobis distance.
#'
#' @references
#' Mahalanobis, P. C. (1936). On the generalised distance in statistics. Proceedings of the National Institute of Sciences of India, 2, 49-55.
#'
#' Piñeiro-Lamas, B. (2024). High dimensional single-index mixture cure models [PhD thesis]. Universidade da Coruña. Available at \url{https://ruc.udc.es/dspace/handle/2183/37035}
#'
#' @export
#'
#' @seealso \code{\link[sicure]{cd.uncured}}, \code{\link[sicure]{kNN.Mahalanobis}}
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
#' suppressWarnings(cd.sm.uncured(x, t, d, log(0.5), log(0.5), k=10))
cd.sm.uncured <- function(x, time, delta, logh3, logh4, k=10){
  n <- length(x)

  # Compute the initial conditional density estimates
  conditional_density_CV_vector <- cd.uncured(x, time, delta, logh3, logh4)

  # Check if any estimates are zero for uncensored data (delta == 1)
  # (it is a problem for camputing the loglikelihood)
  num0 <- sum(conditional_density_CV_vector[delta==1]==0, na.rm = TRUE)
  if (num0==0){ # If there are no zero estimates, no adjustment is needed
    conditional_density_CV_kNN_Mahalanobis_vector <- conditional_density_CV_vector
  } else {  # If there are zero estimates, use kNN Mahalanobis adjustment
    kNN_Mahalanobis_matrix <- kNN.Mahalanobis(x, time, k) # matrix with the indices of the k nearest neighbors
    conditional_density_kNN <- matrix(conditional_density_CV_vector[kNN_Mahalanobis_matrix], n, k)
    # Calculate the mean of the conditional densities from the k-nearest neighbors
    conditional_density_CV_kNN_Mahalanobis_vector <- apply(conditional_density_kNN, 1, mean)

    # Increase k, if necessary, until there are no zero estimates
    while(sum(conditional_density_CV_kNN_Mahalanobis_vector[delta==1]==0, na.rm = TRUE)!=0){
      k <- k+1
      kNN_Mahalanobis_matrix <- kNN.Mahalanobis(x, time, k)
      conditional_density_kNN <- matrix(conditional_density_CV_vector[kNN_Mahalanobis_matrix], n, k)
      conditional_density_CV_kNN_Mahalanobis_vector <- apply(conditional_density_kNN, 1, mean)
    }

    # Replace the zero estimates for uncensored data with the adjusted estimates
    conditional_density_CV_kNN_Mahalanobis_vector_complete <- conditional_density_CV_kNN_Mahalanobis_vector
    conditional_density_CV_kNN_Mahalanobis_vector <-  conditional_density_CV_vector
    indices_to_replace <- which(delta == 1 & conditional_density_CV_vector == 0)
    conditional_density_CV_kNN_Mahalanobis_vector[indices_to_replace] <- conditional_density_CV_kNN_Mahalanobis_vector_complete[indices_to_replace]
  }
  return(conditional_density_CV_kNN_Mahalanobis_vector) # Return the adjusted (or smoothed) conditional density vector
}


