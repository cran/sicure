#' Logarithm of the likelihood of a single-index mixture cure model
#'
#' This function computes the logarithm of the likelihood of a single-index mixture cure model.
#'
#' @param x A numeric vector giving the covariate values.
#' @param time A numeric vector giving the observed times.
#' @param delta A numeric vector giving the values of the uncensoring indicator, where 1 indicates that the event of interest has been observed and 0 indicates that the observation is censored.
#' @param logh1 The logarithm of the bandwidth used to smooth the covariate in the nonparametric estimation of the probability of cure.
#' @param logh2 The logarithm of the bandwidth used to smooth the covariate in the nonparametric estimation of the latency.
#' @param logh3 The logarithm of the bandwidth used to smooth the time variable in the nonparametric estimation of the conditional density of susceptible individuals.
#' @param logh4 The logarithm of the bandwidth used to smooth the covariate in the nonparametric estimation of the conditional density of susceptible individuals.
#' @param r The radius of the moving window used to smooth the nonparametric estimation of the probability of cure.
#' @param k The number of nearest neighbors used to smooth the nonparametric estimation of the conditional density of susceptible individuals.
#'
#' @return A list with two components:
#' \itemize{
#' \item The value of the negative log-likelihood.
#' \item The \eqn{n} terms that contribute to the negative log-likelihood.
#' }
loglik.simcm <- function(x, time, delta, logh1, logh2, logh3, logh4, r, k=10){
  n <- length(x)  # Number of individuals
  term=numeric(n) # Initialize a vector to store log-likelihood terms

  # Estimate the probability of being susceptible (pi)
  pi_est_vector=1-probcure.sm(x, time, delta, logh1, r)[[2]]

  # Estimate the conditional density of susceptible individuals using k-Nearest Neighbors with Mahalanobis distance
  conditional_density_CV_vector_kNN_Mahalanobis=cd.sm.uncured(x,time, delta, logh3, logh4, k)

  for (i in 1:n){
    dat_i <- data.frame(x, time, delta)[-i,]

    # Get the estimated probability of being susceptible for the i-th individual
    pi_est=pi_est_vector[i]; pi_est <- max(pi_est,0)  # Ensure pi_est is non-negative


    if (data.frame(x, time, delta)[i,3]==0){
      # If the event is censored (delta = 0), estimate S0 for the i-th individual
      S0_est=as.numeric(npcure::latency(x, time, delta, dat_i, x0=data.frame(x, time, delta)[i,1], h=exp(logh2), local=FALSE, testimate=data.frame(x, time, delta)[i,2])$S[[1]])
      term[i]=max(log(1-pi_est + pi_est * S0_est), log(10^(-323))) # Avoid log(0)
    } else {
      # If the event is uncensored (delta = 1), estimate f0 for the i-th individual
      f0_est=conditional_density_CV_vector_kNN_Mahalanobis[i]
      f0_est=max(0, f0_est) # Ensure f0_est is non-negative
      term[i]=max(log(pi_est), log(10^(-323))) + max(log(f0_est), log(10^(-323)))
    }

  }
  # Return the negative sum of loglikelihood terms and the terms themselves
  return(list(-sum(term), term))
}


#' Objective function
#'
#' This function computes the negative log-likelihood for a given set of parameters.
#'
#' @param theta_h A numeric vector containing the initial iterant for the vector of parameters and the initial bandwidths.
#' @param x_cov A matrix or data frame giving the covariate values. Each row represents an individual and each column corresponds to a variable.
#' @param time A numeric vector giving the observed times.
#' @param delta A numeric vector giving the values of the uncensoring indicator, where 1 indicates that the event of interest has been observed and 0 indicates that the observation is censored.
#'
#' @return
#' The value of the negative log-likelihood.
fun.opt <- function(theta_h, x_cov, time, delta){
  # theta_h = (thetaini[-1], logh1, logh2, logh3, logh4)
  # d is the number of covariates minus 1
  d = ncol(x_cov)-1

  x_LC = colSums(c(1,theta_h[1:d])*t(x_cov))  # LC = linear combination ; x_LC has dimension 1

  dat_xtd <- as.data.frame(cbind(x_LC, time, delta))
  ot_dat_xtd <- dat_xtd[order(dat_xtd[, 2]), ]     # order the data by time
  time <- ot_dat_xtd[,2]; delta <- ot_dat_xtd[,3];  x <- ot_dat_xtd[,1]

  # Compute -loglikelihood using loglik.simcm function
  loglik=loglik.simcm(x, time, delta, theta_h[d+1], theta_h[d+2], theta_h[d+3], theta_h[d+4], r=2, k=10)[[1]]

  return(loglik)
}



#' Estimation of the vector of parameters in a single-index mixture cure model
#' with a vector of covariates
#'
#' This function provides the estimation of the vector of parameters (\eqn{\boldsymbol{\theta_0}})  in a single-index mixture cure model with a vector of covariates (see Piñeiro-Lamas, 2024, Section 3.1, pages 37-38).
#'
#' @param x_cov A matrix or data frame giving the covariate values. Each row represents an individual and each column corresponds to a variable.
#' @param time A numeric vector giving the observed times.
#' @param delta A numeric vector giving the values of the uncensoring indicator, where 1 indicates that the event of interest has been observed and 0 indicates that the observation is censored.
#' @param randomsearch A logical value, \code{FALSE} by default, specifying whether a random search of the initial iterant is considered.
#'
#' @details
#' The vector of parameters, \eqn{\boldsymbol{\theta_0}}, is estimated by maximum likelihood, with the likelihood function also depending on four bandwidths.
#' Since the parameters cannot be obtained directly from the likelihood, the estimation is performed using numerical optimization with the Nelder-Mead method.
#' To begin the optimization, an initial iterant is required. To select a starting point for \eqn{\boldsymbol{\theta_0}}, a logistic regression model is fitted using the uncensoring indicator \code{delta} as the response variable and the covariates as predictors.
#' Regarding the initial bandwidths, for \eqn{h_3} the rule of thumb bandwidth selector used to estimate the density of the time variable is considered. For \eqn{h_1}, \eqn{h_2} and \eqn{h_4}, the rule of thumb bandwidth used to estimate the density of \eqn{\boldsymbol{\theta_{0,ini}}' \boldsymbol{X}}, where \eqn{\boldsymbol{\theta_{0,ini}}} is the initial estimate of \eqn{\boldsymbol{\theta_0}} and \eqn{\boldsymbol{X}} is the vector of covariates.
#' If \code{randomsearch = TRUE}, nine additional starting points are generated by adding a Gaussian noise (mean 0, standard deviation 0.25) to the logistic regression-based iterant.
#' Each of these ten starting points is used to estimate the vector of parameters of the single-index mixture cure model, and the one that gives place to the lowest negative log-likelihood value is selected.
#'
#' @return A list with the following components:
#' \item{par}{A numeric vector of the estimated parameters. The last four correspond to the logarithms of the bandwidths.}
#' \item{value}{The value of the objective function (negative log-likelihood) at the estimated parameters.}
#' \item{si}{The estimated single-index variable.}
#'
#' @references
#' Piñeiro-Lamas, B. (2024). High dimensional single-index mixture cure models [PhD thesis]. Universidade da Coruña. Available at \url{https://ruc.udc.es/dspace/handle/2183/37035}
#'
#' @export
#'
#' @seealso \code{\link[sicure]{sicure.f}}, \code{\link[sicure]{sicure.vf}}
#'
#' @examples
#' # Some artificial data
#' set.seed(123)
#' n <- 50
#' mix1a<-rnorm(n,mean=0,sd=1); mix1b<-rnorm(n,mean=0.25,sd=sqrt(2)); alf1<-rbinom(n,1,0.2)
#' mix2a<-rnorm(n,mean=0,sd=1); mix2b<-rnorm(n,mean=0.25,sd=sqrt(2)); alf2<-rbinom(n,1,0.2)
#' mix1<-alf1*mix1a+(1-alf1)*mix1b; mix2<-alf2*mix2a+(1-alf2)*mix2b
#' x_cov<-array(c(mix1,mix2),dim=c(n,2)) # Matrix of covariate values
#' theta<-c(1,1.2)
#' Z<-colSums(theta*t(x_cov))
#' y<-Z+rnorm(n,sd=sqrt(abs(Z))) # True lifetimes
#' # Probability of being susceptible
#' p_fun <- function(x){ 0.55 * exp(1.5*x+1.5)/(1+exp(1.5*x+1.5)) + 0.001 }
#' for (i in 1:n){
#'    w <- runif(1)
#'    if (w > p_fun(Z[i])) y[i] <- Inf
#' }
#' c<-rexp(n,rate=0.98) # Censoring values
#' t<-pmin(y,c) # Observed times
#' d = 1 * (y<=c) # Uncensoring indicator
#' \donttest{
#' suppressWarnings(sicure.v(x_cov, t, d))
#' }
#' @importFrom stats binomial rnorm
sicure.v <- function(x_cov, time, delta, randomsearch=FALSE){

  x_cov <- as.data.frame(scale(x_cov))
  dat <- cbind(x_cov, delta)

  # Selection of the initial iterant for theta:
  model <- stats::glm(delta ~., data = dat, family = binomial)
  # Omit the intercept and force the first component to be one
  initial.it <- summary(model)$coef[-1,1] / summary(model)$coef[2,1]
  theta <- initial.it

  # Omit individuals with NA values
  ind <- stats::complete.cases(x_cov)
  x_cov <- x_cov[ind,]
  Z <- colSums(theta*t(x_cov)) # theta'X_cov (theta = initial value)
  time <- time[ind]
  delta <- delta[ind]
  n <- length(Z)

  if(randomsearch==FALSE){

    # Calculate initial (log)bandwidths
    #print('Initial iterant:')
    density_Z_bw_log <- log(stats::density(Z)$bw)
    density_time_bw_log <- log(stats::density(time)$bw)
    #print(c(initial.it[-1], density_Z_bw_log, density_Z_bw_log, density_time_bw_log, density_Z_bw_log))

    # Optimize the parameters
    res=stats::optim(c(initial.it[-1], density_Z_bw_log, density_Z_bw_log, density_time_bw_log, density_Z_bw_log), fun.opt, x_cov = x_cov, time = time, delta = delta, method='Nelder-Mead', control=list(reltol=1e-3, maxit=100))
    si <- colSums(c(1,as.numeric(res$par[1:(length(res$par) - 4)]))*t(x_cov))
    listt <- list(res$par, res$value, si)
    names(listt) <- c("par", "value", "si")
    return(listt)

  }

  if(randomsearch==TRUE){

    matrix.noise <- matrix(rnorm(length(initial.it[-1])*9,mean=0, sd=0.25), ncol=9)
    matrix.it <- matrix.noise + initial.it[-1]
    matrix.it <- cbind(initial.it[-1], matrix.it)
    matrix.it <- rbind(rep(1,10), matrix.it)

    dx <- ncol(x_cov)

    matrix.res <- matrix(NA, nrow=10, ncol=dx+4)
    for (i in 1:10){
      Z <- colSums(matrix.it[,i]*t(x_cov))
      density_Z_bw_log <- log(stats::density(Z)$bw)
      density_time_bw_log <- log(stats::density(time)$bw)
      res=stats::optim(c(matrix.it[-1,i], density_Z_bw_log, density_Z_bw_log, density_time_bw_log, density_Z_bw_log), fun.opt, x_cov = x_cov, time = time, delta = delta, method='Nelder-Mead', control=list(reltol=1e-3, maxit=100))
      matrix.res[i,] <- c(res$par, res$value)
    }
    theta_est <- c(1,matrix.res[which.min(matrix.res[,dx+4]), 1:(dx-1)])
    si <- colSums(theta_est*t(x_cov))
    listt <- list(matrix.res[which.min(matrix.res[,dx+4]), 1:(dx+3)], matrix.res[which.min(matrix.res[,dx+4]), dx+4], si)
    names(listt) <- c("par", "value", "si")
    return(listt)

  }

}


#' Estimation of the vector of parameters in a single-index mixture cure model
#' with a functional covariate
#'
#' This function provides the estimation of the vector of parameters in a single-index mixture cure model
#' with a functional covariate (see Piñeiro-Lamas, 2024, Section 4.1, pages 83-84).
#' A Functional Principal Components Analysis (FPCA) representation that explains at least the (\code{propvar}\eqn{*100})\%
#' of the variability of the data is considered (for more details, see Ramsay and Silverman, 2005).
#'
#' @param x_cov A matrix or data frame \eqn{n} x \eqn{m} giving the functional covariate values. Each row represents an individual (a curve); \eqn{m} is the number of observed points in each curve.
#' @param time A numeric vector giving the observed times.
#' @param delta A numeric vector giving the values of the uncensoring indicator, where 1 indicates that the event of interest has been observed and 0 indicates that the observation is censored.
#' @param propvar Minimum proportion of explained variability with the FPCA representation.
#' @param randomsearch A logical value, \code{FALSE} by default, specifying whether a random search of the initial iterant is considered.
#'
#' @details
#' The infinite-dimensional nature of the functional data is reduced via FPCA. This basis representation is then truncated, reducing the dimension to \eqn{K}, where each functional observation is summarized into a vector of scores, \eqn{(\xi_1, \xi_2, \dots, \xi_K)}.
#' After this reduction, the model can be treated similarly to the vector covariate case.
#' For more details on the estimation process and the specific arguments, see \code{\link[sicure]{sicure.v}} function, which focuses on single-index mixture cure models with a vector of covariates.
#'
#' @return A list with the following components:
#' \item{par}{A numeric vector of the estimated parameters. The last four correspond to the logarithms of the bandwidths.}
#' \item{value}{The value of the objective function (negative log-likelihood) at the estimated parameters.}
#' \item{si}{The estimated single-index variable.}
#'
#' @references
#' Piñeiro-Lamas, B. (2024). High dimensional single-index mixture cure models [PhD thesis]. Universidade da Coruña. Available at \url{https://ruc.udc.es/dspace/handle/2183/37035}
#'
#' Ramsay, J. O., and Silverman, B. W. (2005). Functional Data Analysis, 2nd ed., Springer, New York.
#'
#' @export
#'
#' @seealso \code{\link[sicure]{sicure.v}}, \code{\link[sicure]{sicure.vf}}
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
#' delta  <- ifelse(u < p, ifelse(y < c, 1, 0), 0) # Uncensoring indicator
#' # Number of individuals (rows)
#' n <- 50
#' # Numbers of observations per individual (columns)
#' m <- 100
#' # Observation times (between 0 and 1)
#' x <- seq(0, 1, length.out = m)
#' # Auxiliar function to simulate the other functions by adding some noise
#' # Shift controls the horizontal displacement of the functions
#' sim_func <- function(x, shift, sd_noise) {
#'   # positive-negative-negative waves
#'   sin(2*pi*(x + shift))+sin(4*pi*(x + shift))-sin(6*pi*(x + shift))+rnorm(m, 0, sd_noise)
#' }
#' # Simulated functions
#' data_matrix <- matrix(NA, nrow=n, ncol=m)
#' for (i in 1:n) {
#'   shift <- runif(1, -0.05, 0.05)
#'   data_matrix[i, ] <- sim_func(x, shift, sd_noise = 0.03)
#' }
#' matplot(x, t(data_matrix), type = "l", lty = 1, ylab='f(x)')
#' \donttest{
#' suppressWarnings(sicure.f(data_matrix, t, delta, 0.9))
#' }
#' @importFrom stats binomial rnorm
sicure.f <- function(x_cov, time, delta, propvar=0.9, randomsearch=FALSE){

  # Perform FPCA with an initial high number of components
  fun_pca <- fda::pca.fd(fda::fd(t(as.matrix(x_cov))), nharm=20)
  var_exp <- fun_pca$varprop ; cum_var_exp <- cumsum(var_exp)

  # Determine the optimal number of components to explain at least (propvar*100)% of variability
  nharm_opt <- max(which(cum_var_exp >= propvar)[1],2)
  message("The optimal number of basis components that explain at least ", propvar * 100, "% of the variability of the data is: ", nharm_opt)

  # Perform PCA again with the optimal number of components
  fun_pca <- fda::pca.fd(fda::fd(t(as.matrix(x_cov))), nharm=nharm_opt)
  x_cov <- fun_pca$scores
  x_cov <- as.data.frame(scale(x_cov))
  dat <- cbind(x_cov, delta)

  # Selection of the initial iterant for theta:
  model <- stats::glm(delta ~., data = dat, family = binomial)
  # Omit the intercept and force the first component to be one
  initial.it <- summary(model)$coef[-1,1] / summary(model)$coef[2,1]
  theta <- initial.it

  # Omit individuals with NA values
  ind <- stats::complete.cases(x_cov)
  x_cov <- x_cov[ind,]
  Z <- colSums(theta*t(x_cov)) # theta'X_cov (theta = initial value)
  time <- time[ind]
  delta <- delta[ind]
  n <- length(Z)

  if(randomsearch==FALSE){

    # Calculate initial (log)bandwidths
    #print('Initial iterant:')
    density_Z_bw_log <- log(stats::density(Z)$bw)
    density_time_bw_log <- log(stats::density(time)$bw)
    #print(c(initial.it[-1], density_Z_bw_log, density_Z_bw_log, density_time_bw_log, density_Z_bw_log))

    # Optimize the parameters
    res=stats::optim(c(initial.it[-1], density_Z_bw_log, density_Z_bw_log, density_time_bw_log, density_Z_bw_log), fun.opt, x_cov = x_cov, time = time, delta = delta, method='Nelder-Mead', control=list(reltol=1e-3, maxit=100))
    si <- colSums(c(1,as.numeric(res$par[1:(length(res$par) - 4)]))*t(x_cov))
    listt <- list(res$par, res$value, si)
    return(listt)

  }

  if(randomsearch==TRUE){

    matrix.noise <- matrix(rnorm(length(initial.it[-1])*9,mean=0, sd=0.25), ncol=9)
    matrix.it <- matrix.noise + initial.it[-1]
    matrix.it <- cbind(initial.it[-1], matrix.it)
    matrix.it <- rbind(rep(1,10), matrix.it)

    dx <- ncol(x_cov)

    matrix.res <- matrix(NA, nrow=10, ncol=dx+4)
    for (i in 1:10){
      Z <- colSums(matrix.it[,i]*t(x_cov))
      density_Z_bw_log <- log(stats::density(Z)$bw)
      density_time_bw_log <- log(stats::density(time)$bw)
      res=stats::optim(c(matrix.it[-1,i], density_Z_bw_log, density_Z_bw_log, density_time_bw_log, density_Z_bw_log), fun.opt, x_cov = x_cov, time = time, delta = delta, method='Nelder-Mead', control=list(reltol=1e-3, maxit=100))
      matrix.res[i,] <- c(res$par, res$value)
    }
    theta_est <- c(1,matrix.res[which.min(matrix.res[,dx+4]), 1:(dx-1)])
    si <- colSums(theta_est*t(x_cov))
    listt <- list(matrix.res[which.min(matrix.res[,dx+4]), 1:(dx+3)], matrix.res[which.min(matrix.res[,dx+4]), dx+4], si)
    names(listt) <- c("par", "value", "si")
    return(listt)


  }

}


#' Estimation of the vector of parameters in a single-index mixture cure model
#' with a vector and a functional covariate
#'
#' This function provides the estimation of the vector of parameters in a single-index mixture cure model
#' with a vector and a functional covariate (see Piñeiro-Lamas, 2024, Section 5.1, page 99).
#' A Functional Principal Components Analysis (FPCA) representation that explains at least the (\code{propvar}\eqn{*100})\%
#' of the variability of the functional data is considered (for more details, see Ramsay and Silverman, 2005).
#'
#' @param x_cov_v A matrix or data frame giving the vector covariate values. Each row represents an individual and each column corresponds to a variable.
#' @param x_cov_f A matrix or data frame \eqn{n} x \eqn{m} giving the functional covariate values. Each row represents an individual (a curve); \eqn{m} is the number of observed points in each curve.
#' @param time A numeric vector giving the observed times.
#' @param delta A numeric vector giving the values of the uncensoring indicator, where 1 indicates that the event of interest has been observed and 0 indicates that the observation is censored.
#' @param propvar Minimum proportion of explained variability with the FPCA representation.
#' @param randomsearch A logical value, \code{FALSE} by default, specifying whether a random search of the initial iterant is considered.
#'
#' @details
#' The infinite-dimensional nature of the functional data is reduced via FPCA. This basis representation is then truncated, reducing the dimension to \eqn{K}, where each functional observation is summarized into a vector of scores, \eqn{(\xi_1, \xi_2, \dots, \xi_K)}.
#' Once this reduction is performed, if the vector covariate has dimension \eqn{d}, a combined joint vector variable can be constructed by considering both the vector covariate and the functional scores, resulting in a total dimension of \eqn{d + K}.
#' This joint variable can then be analyzed within the framework of a single-index mixture cure model with a vector of covariates.
#' For more details on the estimation process and the specific arguments, see \code{\link[sicure]{sicure.v}} function, which focuses on single-index mixture cure models with a vector of covariates.
#'
#' @return A list with the following components:
#' \item{par}{A numeric vector of the estimated parameters. The last four correspond to the logarithms of the bandwidths.}
#' \item{value}{The value of the objective function (negative log-likelihood) at the estimated parameters.}
#' \item{si}{The estimated single-index variable.}
#'
#' @references
#' Piñeiro-Lamas, B. (2024). High dimensional single-index mixture cure models [PhD thesis]. Universidade da Coruña. Available at \url{https://ruc.udc.es/dspace/handle/2183/37035}
#'
#' Ramsay, J. O., and Silverman, B. W. (2005). Functional Data Analysis, 2nd ed., Springer, New York.
#'
#' @export
#'
#' @seealso \code{\link[sicure]{sicure.v}}, \code{\link[sicure]{sicure.f}}
#'
#' @examples
#' # Some artificial data
#' set.seed(123)
#' n <- 50
#' mix1a<-rnorm(n,mean=0,sd=1); mix1b<-rnorm(n,mean=0.25,sd=sqrt(2)); alf1<-rbinom(n,1,0.2)
#' mix2a<-rnorm(n,mean=0,sd=1); mix2b<-rnorm(n,mean=0.25,sd=sqrt(2)); alf2<-rbinom(n,1,0.2)
#' mix1<-alf1*mix1a+(1-alf1)*mix1b; mix2<-alf2*mix2a+(1-alf2)*mix2b
#' x_cov_v<-array(c(mix1,mix2),dim=c(n,2)) # Matrix of covariate values
#' theta<-c(1,1.2)
#' Z<-colSums(theta*t(x_cov_v))
#' y<-Z+rnorm(n,sd=sqrt(abs(Z))) # True lifetimes
#' # Probability of being susceptible
#' p_fun <- function(x){ 0.55 * exp(1.5*x+1.5)/(1+exp(1.5*x+1.5)) + 0.001 }
#' for (i in 1:n){
#'    w <- runif(1)
#'    if (w > p_fun(Z[i])) y[i] <- Inf
#' }
#' c<-rexp(n,rate=0.98) # Censoring values
#' t<-pmin(y,c) # Observed times
#' d = 1 * (y<=c) # Uncensoring indicator
#' # Functional covariate:
#' # Number of individuals (rows)
#' n <- 50
#' # Numbers of observations per individual (columns)
#' m <- 100
#' # Observation times (between 0 and 1)
#' x <- seq(0, 1, length.out = m)
#' # Auxiliar function to simulate the other functions by adding some noise
#' # Shift controls the horizontal displacement of the functions
#' sim_func <- function(x, shift, sd_noise) {
#'   # positive-negative-negative waves
#'   sin(2*pi*(x + shift))+sin(4*pi*(x + shift))-sin(6*pi*(x + shift))+rnorm(m, 0, sd_noise)
#' }
#' # Simulated functions
#' data_matrix <- matrix(NA, nrow=n, ncol=m)
#' for (i in 1:n) {
#'   shift <- runif(1, -0.05, 0.05)
#'   data_matrix[i, ] <- sim_func(x, shift, sd_noise = 0.03)
#' }
#' matplot(x, t(data_matrix), type = "l", lty = 1, ylab='f(x)')
#' \donttest{
#' suppressWarnings(sicure.vf(x_cov_v, data_matrix, t, d, 0.9))
#' }
#' @importFrom stats binomial rnorm
sicure.vf <- function(x_cov_v, x_cov_f, time, delta, propvar=0.9, randomsearch=FALSE){

  # Vector covariate:
  x_cov_v <- as.data.frame(scale(x_cov_v))
  colnames(x_cov_v) <- paste0("V_", seq_len(ncol(x_cov_v)))

  # Functional covariate:
  # Perform FPCA with an initial high number of components
  fun_pca <- fda::pca.fd(fda::fd(t(as.matrix(x_cov_f))), nharm=20)
  var_exp <- fun_pca$varprop ; cum_var_exp <- cumsum(var_exp)

  # Determine the optimal number of components to explain at least (propvar*100)% of variability
  nharm_opt <- max(which(cum_var_exp >= propvar)[1],2)
  message("The optimal number of basis components that explain at least ", propvar * 100, "% of the variability of the data is: ", nharm_opt)

  # Perform PCA again with the optimal number of components
  fun_pca <- fda::pca.fd(fda::fd(t(as.matrix(x_cov_f))), nharm=nharm_opt)
  x_cov_f <- fun_pca$scores
  x_cov_f <- as.data.frame(scale(x_cov_f))
  colnames(x_cov_f) <- paste0("F_", seq_len(ncol(x_cov_f)))

  x_cov <- cbind(x_cov_v, x_cov_f)
  dat <- cbind(x_cov, delta)

  # Selection of the initial iterant for theta:
  model <- stats::glm(delta ~., data = dat, family = binomial)
  # Omit the intercept and force the first component to be one
  initial.it <- summary(model)$coef[-1,1] / summary(model)$coef[2,1]
  theta <- initial.it

  # Omit individuals with NA values
  ind <- stats::complete.cases(x_cov)
  x_cov <- x_cov[ind,]
  Z <- colSums(theta*t(x_cov)) # theta'X_cov (theta = initial value)
  time <- time[ind]
  delta <- delta[ind]
  n <- length(Z)

  if(randomsearch==FALSE){

    # Calculate initial (log)bandwidths
    #print('Initial iterant:')
    density_Z_bw_log <- log(stats::density(Z)$bw)
    density_time_bw_log <- log(stats::density(time)$bw)
    #print(c(initial.it[-1], density_Z_bw_log, density_Z_bw_log, density_time_bw_log, density_Z_bw_log))

    # Optimize the parameters
    res=stats::optim(c(initial.it[-1], density_Z_bw_log, density_Z_bw_log, density_time_bw_log, density_Z_bw_log), fun.opt, x_cov = x_cov, time = time, delta = delta, method='Nelder-Mead', control=list(reltol=1e-3, maxit=100))
    si <- colSums(c(1,as.numeric(res$par[1:(length(res$par) - 4)]))*t(x_cov))
    listt <- list(res$par, res$value, si)
    return(listt)

  }

  if(randomsearch==TRUE){

    matrix.noise <- matrix(rnorm(length(initial.it[-1])*9,mean=0, sd=0.25), ncol=9)
    matrix.it <- matrix.noise + initial.it[-1]
    matrix.it <- cbind(initial.it[-1], matrix.it)
    matrix.it <- rbind(rep(1,10), matrix.it)

    dx <- ncol(x_cov)

    matrix.res <- matrix(NA, nrow=10, ncol=dx+4)
    for (i in 1:10){
      Z <- colSums(matrix.it[,i]*t(x_cov))
      density_Z_bw_log <- log(stats::density(Z)$bw)
      density_time_bw_log <- log(stats::density(time)$bw)
      res=stats::optim(c(matrix.it[-1,i], density_Z_bw_log, density_Z_bw_log, density_time_bw_log, density_Z_bw_log), fun.opt, x_cov = x_cov, time = time, delta = delta, method='Nelder-Mead', control=list(reltol=1e-3, maxit=100))
      matrix.res[i,] <- c(res$par, res$value)
    }
    theta_est <- c(1,matrix.res[which.min(matrix.res[,dx+4]), 1:(dx-1)])
    si <- colSums(theta_est*t(x_cov))
    listt <- list(matrix.res[which.min(matrix.res[,dx+4]), 1:(dx+3)], matrix.res[which.min(matrix.res[,dx+4]), dx+4], si)
    names(listt) <- c("par", "value", "si")
    return(listt)


  }

}
