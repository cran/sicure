% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sicure.R
\name{sicure.f}
\alias{sicure.f}
\title{Estimation of the vector of parameters in a single-index mixture cure model
with a functional covariate}
\usage{
sicure.f(x_cov, time, delta, propvar = 0.9, randomsearch = FALSE)
}
\arguments{
\item{x_cov}{A matrix or data frame \eqn{n} x \eqn{m} giving the functional covariate values. Each row represents an individual (a curve); \eqn{m} is the number of observed points in each curve.}

\item{time}{A numeric vector giving the observed times.}

\item{delta}{A numeric vector giving the values of the uncensoring indicator, where 1 indicates that the event of interest has been observed and 0 indicates that the observation is censored.}

\item{propvar}{Minimum proportion of explained variability with the FPCA representation.}

\item{randomsearch}{A logical value, \code{FALSE} by default, specifying whether a random search of the initial iterant is considered.}
}
\value{
A list with the following components:
\item{par}{A numeric vector of the estimated parameters. The last four correspond to the logarithms of the bandwidths.}
\item{value}{The value of the objective function (negative log-likelihood) at the estimated parameters.}
\item{si}{The estimated single-index variable.}
}
\description{
This function provides the estimation of the vector of parameters in a single-index mixture cure model
with a functional covariate (see Piñeiro-Lamas, 2024, Section 4.1, pages 83-84).
A Functional Principal Components Analysis (FPCA) representation that explains at least the (\code{propvar}\eqn{*100})\%
of the variability of the data is considered (for more details, see Ramsay and Silverman, 2005).
}
\details{
The infinite-dimensional nature of the functional data is reduced via FPCA. This basis representation is then truncated, reducing the dimension to \eqn{K}, where each functional observation is summarized into a vector of scores, \eqn{(\xi_1, \xi_2, \dots, \xi_K)}.
After this reduction, the model can be treated similarly to the vector covariate case.
For more details on the estimation process and the specific arguments, see \code{\link[sicure]{sicure.v}} function, which focuses on single-index mixture cure models with a vector of covariates.
}
\examples{
# Some artificial data
set.seed(123)
n <- 50
x <- runif(n, -2, 2) # Covariate values
y <- rweibull(n, shape = 0.5 * (x + 4)) # True lifetimes
c <- rexp(n) # Censoring values
p <- exp(2*x)/(1 + exp(2*x)) # Probability of being susceptible
u <- runif(n)
t  <- ifelse(u < p, pmin(y, c), c) # Observed times
delta  <- ifelse(u < p, ifelse(y < c, 1, 0), 0) # Uncensoring indicator
# Number of individuals (rows)
n <- 50
# Numbers of observations per individual (columns)
m <- 100
# Observation times (between 0 and 1)
x <- seq(0, 1, length.out = m)
# Auxiliar function to simulate the other functions by adding some noise
# Shift controls the horizontal displacement of the functions
sim_func <- function(x, shift, sd_noise) {
  # positive-negative-negative waves
  sin(2*pi*(x + shift))+sin(4*pi*(x + shift))-sin(6*pi*(x + shift))+rnorm(m, 0, sd_noise)
}
# Simulated functions
data_matrix <- matrix(NA, nrow=n, ncol=m)
for (i in 1:n) {
  shift <- runif(1, -0.05, 0.05)
  data_matrix[i, ] <- sim_func(x, shift, sd_noise = 0.03)
}
matplot(x, t(data_matrix), type = "l", lty = 1, ylab='f(x)')
\donttest{
suppressWarnings(sicure.f(data_matrix, t, delta, 0.9))
}
}
\references{
Piñeiro-Lamas, B. (2024). High dimensional single-index mixture cure models [PhD thesis]. Universidade da Coruña. Available at \url{https://ruc.udc.es/dspace/handle/2183/37035}

Ramsay, J. O., and Silverman, B. W. (2005). Functional Data Analysis, 2nd ed., Springer, New York.
}
\seealso{
\code{\link[sicure]{sicure.v}}, \code{\link[sicure]{sicure.vf}}
}
