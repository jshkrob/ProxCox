###########################################################################################
# computes proximal for sparse-smooth penalty given coefficients uu, penalty matrix, time-step
# and sparse penalty lambda_1
###########################################################################################
prox <- function(uu, PENALTY, t, lambda1) {
    prxfunc <- function(u) {
        #print(1 - ((t*lambda1) / sqrt( t(u) %*% PENALTY %*% u  )))
        p <- max(0, 1 - ((t*lambda1) / sqrt( t(u) %*% PENALTY %*% u  ))) * u
        return(p)
    }
    prx <- apply(uu, 1, FUN = prxfunc)
    prxx <- as.vector(prx)
    return(prxx)
}
