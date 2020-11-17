###########################################################################################
# Based off of the CoxRidge package, computes non-regularized minimized for B-spline paramatrized
# by formulation in Yi Yang, 2020 using a Newton-Raphson descent method
###########################################################################################
#' @export
reg.fit <- function (time, death, X, Ftime, theta, t_s, iter.max) {
    epsilon <- 5; iter <- 0; at.iter <- 0;
    likelihood_iter <- rep(0, iter.max); epsilon_iter <- rep(0, iter.max);
    ord <- order(time)
    Ftime <- Ftime[ord, ]
    n <- nrow(X); d <- death; p <- ncol(X); q <- ncol(Ftime)
    if(missing(theta)){
        theta <- matrix(rep(0, p * q), ncol = q) # creates p X q matrix of 0s
        BTHETA <- theta
    }
    X <- apply(X, 2, function(x) { x - mean(x)}) # normalize
    eventtimes <- (1:n)[d == 1]
    k <- length(eventtimes)
    JJ <- Ftime
    Ftime <- Ftime[eventtimes, ]
    while (epsilon > 1e-06 && iter < iter.max) {
        th <- theta
        iter <- iter + 1
        # current theta(s-1): gradient calculations on theta(s-1)
        tmp <- sapply(eventtimes, sumevents_reg, X, Ftime, theta, n, p, eventtimes)
        hi <- unlist(tmp[1, ])
        lik <- -(sum(log(hi)) + sum(X[eventtimes, ] * t(theta %*% t(Ftime))))
        scores <- matrix(unlist(tmp[2, ]), ncol = k) %*% rep.int(1, k) + matrix(t(X[eventtimes, ]) %*% Ftime, ncol = 1)
        hess <- round(-(matrix(matrix(unlist(tmp[3, ]), ncol = k) %*% rep(1, k), ncol = p * q, byrow = FALSE)), digits = 6 )
        theta <- matrix(matrix(theta, ncol = 1) - t_s * (ginv(hess) %*% matrix(t(scores), ncol = 1)), ncol = q)
        epsilon <- max(abs(theta - th))
        epsilon_iter[iter] <- epsilon
        likelihood_iter[iter] <- lik
        cat('Likelihood:', lik, "\n")
        cat("Epsilon:", epsilon, "\nIteration:", iter, "\n")
        if (iter > iter.max) {
            stop("Likelihood did not converge after ", iter, " iterations")
        }

    }
    # sterr <- round(matrix(sqrt(abs(diag(solve(hess)))), nrow = p), 5)
    sterr <- 0
    outlist <- list(time = time, X = X, Ft = JJ, Ftime = Ftime, death = death, theta = theta,
                    se = sterr, lik = lik, iter = iter, hi = hi,
                    likelihood_iter = likelihood_iter, epsilon_iter = epsilon_iter)
    outlist
}
