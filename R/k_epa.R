#' Rescaled Epanechnikov kernel
#'
#' This function computes the rescaled Epanechnikov kernel for a given value \code{s} and bandwidth \code{g}:
#' \deqn{K\left(\frac{s}{g}\right) = \frac{3}{4g} \left( 1 - \left(\frac{s}{g}\right)^2 \right) I\left(\left|\frac{s}{g}\right| \leq 1\right),}
#' where \eqn{I} is the indicator function (it takes the value 1 if the condition is true and 0 otherwise).
#'
#' @param g A numeric value or vector which contains the bandwidth(s).
#' @param s A numeric value or vector with the values of the variable in which the kernel will be evaluated.
#'
#' @return Rescaled Epanechnikov kernel for the given value \code{s} and bandwidth \code{g}.
#' @export
#'
#' @examples
#' k.epa(g=5,s=2)
#' k.epa(g=5,s=c(2,1.5))
#' k.epa(g=c(5,6),s=c(2,1.5))
k.epa <- function(g,s){
  return(3/(4*g) *(1-s^2/g^2)*c((s/g)<=1)*c((s/g)>=(-1))) # Epanechnikov
}

