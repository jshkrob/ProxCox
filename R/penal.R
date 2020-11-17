
#' Computes matrices used for weighted group penalty
#'
#' This function takes the number of examples, a sequence of observed failure times t and penalization parameter
#' lambda2. Note that lambda2 controls the smoothness of the splines since we are regularizing
#' the second derivatives of B-spline basis functions. The output will be a t x t matrix of "weights"
#'
#' @param n size of sample data set
#' @param t survival times
#' @param lambda2 smoothness parameter
#' @return penalization matrices
#' @export
penal <- function(n, t, lambda2) {
    K <- length(t)
    or <- t[order(t)]
    h <- apply(cbind(or[2:K], or[1:(K-1)]), MARGIN = 1, function(x) x[1] - x[2])
    A_11 <- matrix(data = rep(0, K^2), ncol = K, nrow = K)
    A_12 <- matrix(data = rep(0, K^2), ncol = K, nrow = K)
    A_22 <- matrix(data = rep(0, K^2), ncol = K, nrow = K)
    A_11[1,1] <- h[1]/3; A_11[K,K] <- h[K-1]/3;
    A_12[1,1] <- (h[1]^3)/(-45); A_12[K,K] <- -( (h[K-1]^3) / 45)
    A_22[1,1] <- (4/315)*(h[1]^5); A_22[K,K] <- -(h[K-1]^5) * (4/315)
    for (i in 2:(K-1)) {
        A_11[i,i] <- (h[i-1]/3) + (h[i]/3)
        A_12[i,i] <- -(h[i-1]^3/45) - (h[i]^3/45)
        A_22[i,i] <- (4/315) * (h[i-1]^5 + h[i]^5)
    }
    for (i in 1:K) {
        A_11[i,i-1] <- (h[i-1]/6); A_11[i-1,i] <- (h[i-1]/6)
        A_12[i, i-1] <- (-7/360) * (h[i-1]^3); A_12[i-1, i] <- (-7/360) * (h[i-1]^3);
        A_22[i, i-1] <- (31/15120) * (h[i-1]^5); A_22[i-1, i] <- (31/15120) * (h[i-1]^5);
    }
    # Mapping matrix F
    FF <- matrix(data = rep(0, K^2), ncol = K, nrow = K)
    B <- matrix(data = rep(0, (K-2)^2), ncol = K-2, nrow = K-2)
    D <- matrix(data = rep(0, (K-2)*K), ncol = K, nrow = K-2)
    for(i in 1:(K-2)) {
        D[i,i] <- (1/h[i])
        D[i, i+1] <- -(1/h[i]) - (1/h[i+1])
        D[i, i+2] <- (1/h[i+1])
        B[i,i] <- (h[i] + h[i+1]) / 3
        if (i < (K-2)){ B[i, i+1] <- h[i+1]/6; B[i+1, i] <- h[i+1]/6 }
    }

    BB <- solve(B)
    FF = rbind(rep(0, K), BB %*% D, rep(0, K))
    PEN1 <- A_11 + (t(FF) %*% t(A_12)) + (A_12 %*% FF) + (t(FF) %*% A_22 %*% FF)
    PEN2 <- t(D) %*% BB %*% D
    return((PEN1 + lambda2*PEN2))
}
