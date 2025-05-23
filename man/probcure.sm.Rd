% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probcure_sm.R
\name{probcure.sm}
\alias{probcure.sm}
\title{Smoothed version of the nonparametric estimator of the conditional
probability of cure}
\usage{
probcure.sm(x, time, delta, logh1, r)
}
\arguments{
\item{x}{A numeric vector giving the covariate values.}

\item{time}{A numeric vector giving the observed times.}

\item{delta}{A numeric vector giving the values of the uncensoring indicator, where 1 indicates that the event of interest has been observed and 0 indicates that the observation is censored.}

\item{logh1}{The logarithm of the bandwidth for smoothing the covariate.}

\item{r}{Radius of moving window.}
}
\value{
A list with two components:
\itemize{
\item A vector containing the cross-validation estimations of the probability of cure.
\item The previous vector smoothed with \code{\link[caTools]{runmean}} with a moving window of \eqn{k = 2r + 1}.
}
}
\description{
This function computes a smoothed version of the nonparametric estimator
of the probability of cure proposed by Xu and Peng (2014) and deeply studied
by López-Cheda et al. (2017). The smoothing is performed using the \code{\link[caTools]{runmean}}
function, which computes a moving average of the estimated probabilities in a window determined
by a radius \code{r}.
The non-smoothed version is implemented in the \code{\link[npcure]{probcure}} function of the \code{npcure} package (López-Cheda et al., 2021).
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
d  <- ifelse(u < p, ifelse(y < c, 1, 0), 0) # Uncensoring indicator
data <- data.frame(x = x, t = t, d = d)

# Smoothed nonparametric estimates of cure probability with bandwidth=2
q1 <- probcure.sm(x, t, d, logh1 = log(2), r=2)[[2]]
plot(sort(x), q1[order(x)], type = "l", xlab = "Covariate", ylab = "Cure probability",
     ylim = c(0, 1))
}
\references{
López-Cheda, A., Cao, R., Jácome, M. A., Van Keilegom, I. (2017). Nonparametric incidence estimation and bootstrap bandwidth selection in mixture cure models. Computational Statistics & Data Analysis, 105, 144–165. \doi{10.1016/j.csda.2016.08.002}.

López-Cheda, A., Jácome, M. A., López-de-Ullibarri, I. (2021). The R Journal, 13(1), 21-41. \doi{10.32614/RJ-2021-027}.

Xu, J., Peng, Y. (2014). Nonparametric cure rate estimation with covariates. The Canadian Journal of Statistics, 42, 1-17. \doi{10.1002/cjs.11197}.
}
\seealso{
\code{\link[npcure]{probcure}}
}
