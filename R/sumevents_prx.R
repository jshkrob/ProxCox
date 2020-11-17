
#' helper function for Partial Cox Likelihood calculation for proximal method
#' @param i current index
#' @param X feature matrix
#' @param Ftime matrix of values of the basis for certain time points
#' @param theta current parameters of partial Cox likelihood
#' @param n dimension
#' @param p number of variables
#' @param eventtimes indexes of individuals that had event of interest
#'
sumevents_prx <- function (i, X, Ftime, theta, n, p, eventtimes) {
    xj <- matrix(X[i:n, ], nrow = (n - i + 1))
    Fi <- Ftime[eventtimes == i, ] # subset basis for those with the time (t_i)
    rr <- exp(xj %*% theta %*% Fi) # denominator of loss
    hi <- 1/sum(rr)
    xmean <- hi * t(xj) %*% rr # calculate expectation for X_i!
    list(hi, -xmean %*% t(Fi)) # hessian !
}
