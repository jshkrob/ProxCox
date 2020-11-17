###########################################################################################
# Proximal gradient descent algorithm for Cox regression with time-varying effect
# time - failure times
# death - censoring indicator (1 = noncensored,  0 = censored)
# X - covariate matrix
# Ftime - basis of functions evaluated at each time point of `time` vector
# knots - knots of spline `Ftime`
# theta - optional initial value for theta parameters
# iter.max - number of iterations of algorithm before breaking
# lambda1, lambda2 - penalization parameters
# ttime - initial time-step
# acceleration - controls for acceleration (FISTA)
# delta - backtracking stepping size for proximal gd step
###########################################################################################

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param time PARAM_DESCRIPTION
#' @param death PARAM_DESCRIPTION
#' @param X PARAM_DESCRIPTION
#' @param names PARAM_DESCRIPTION, Default: NULL
#' @param Ftime PARAM_DESCRIPTION
#' @param knots PARAM_DESCRIPTION
#' @param theta PARAM_DESCRIPTION
#' @param iter.max PARAM_DESCRIPTION
#' @param lambda1 PARAM_DESCRIPTION
#' @param lambda2 PARAM_DESCRIPTION
#' @param ttime PARAM_DESCRIPTION
#' @param acceleration PARAM_DESCRIPTION, Default: FALSE
#' @param delta PARAM_DESCRIPTION, Default: 0.5
#' @param verbal PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname prox.fit
#' @export
prox.fit <- function (time, death, X, names=NULL, Ftime, knots, theta, iter.max, lambda1, lambda2, ttime, acceleration = FALSE, delta = 0.5, verbal = FALSE) {
    if(!is.null(names)) colnames(X) <- names
    survt <- "Surv(time, death)"
    ff <- as.formula(paste(survt, paste(names, collapse=" + "), sep=" ~ "))

    epsilon <- 5; iter <- 0; t <- ttime; BTHETA <- NULL; BLIKL <- Inf; BPEN <- NULL; at.iter <- 0;
    likelihood_iter <- rep(0, iter.max);
    unpenlikelihood_iter <- rep(0, iter.max)
    epsilon_iter <- rep(0, iter.max);
    ord <- order(time)
    Ftime <- Ftime[ord, ]
    n <- nrow(X); d <- death; p <- ncol(X); q <- ncol(Ftime)
    if(missing(theta)){
        theta <- matrix(rep(0, p * q), ncol = q) # creates p X q matrix of 0s
        # vv <- survival::coxph(ff, data = as.data.frame(cbind(X, time, death)))
        # coef <- vv$coefficients
        # theta[,1]<-coef
        BTHETA <- theta
        # print(theta)
    } else {
        theta <- matrix(rep(0, p * q), ncol = q) # creates p X q matrix of 0s
        BTHETA <- theta + matrix(rnorm(p*q, mean = 0, sd = 0.01), ncol = q)
    }

    if (acceleration == TRUE) { theta_1 <- matrix(rep(0, p*q), ncol = q) }
    X <- apply(X, 2, function(x) { x - mean(x)}) # normalize
    eventtimes <- (1:n)[d == 1]
    k <- length(eventtimes)
    JJ <- Ftime
    Ftime <- Ftime[eventtimes, ]

    # Fts should have knots
    penalty <- penal(n = n, knots, lambda2 = lambda2) # OMEGA_1, OMEGA_2 for each cov (list of matrices)
    if (verbal) { cat("PENALTY: \n"); print(penalty) }

    tt <- numeric(); tt[1] <- 1; for(i in 2:(iter.max+1)) {tt[i] <- (1+sqrt(1 + (4 * tt[i-1]^2)))/2}


    while (epsilon > 1e-06 && iter < iter.max) {
        th <- theta
        LHS <- 1; RHS <- 0;
        iter <- iter + 1
        t_s <- t


        if (iter == 1) {
            # current theta(s-1): gradient calculations on theta(s-1)
            tmp <- sapply(eventtimes, sumevents_c, X, Ftime, theta, n, p, eventtimes)
            hi <- unlist(tmp[1, ])
            PEN <- lambda1 * sum(apply( theta, 1, FUN = function(th_p) sqrt(t(th_p) %*% penalty %*% th_p )))
            lik <- -(sum(log(hi)) + sum(X[eventtimes, ] * t(theta %*% t(Ftime)))) + PEN # previous likelihood

            # minimizing negative log lik
            scores <- matrix(unlist(tmp[2, ]), ncol = k) %*% rep.int(1, k) + matrix(t(X[eventtimes, ]) %*% Ftime, ncol = 1)
            thetanew <- matrix(matrix(theta, ncol = 1) + (t_s * matrix(t(scores), ncol = 1)), ncol = q)
            grad <- matrix(t(scores), ncol = 1)
        } else {
            # minimizing negative log lik
            lik <- -lik_n + PEN
            scores <- matrix(unlist(tmp[2, ]), ncol = k) %*% rep.int(1, k) + matrix(t(X[eventtimes, ]) %*% Ftime, ncol = 1)
            thetanew <- matrix(matrix(theta, ncol = 1) + (t_s * matrix(t(scores), ncol = 1)), ncol = q)
            grad <- matrix(t(scores), ncol = 1)
        }


        #grad <- matrix(t(scores), ncol = 1) # gradient w.r.t 11, 12, ..., 1q, 21, 22, ..., 2q, ...., p*q
        start.time <- Sys.time(); iter.back <- 0;
        while (round(LHS, 6) > round(RHS, 6)) {
            iter.back <- iter.back + 1;
            #thetanew <- matrix(matrix(th, ncol = 1) + (t_s * grad), ncol = q)
            prx <- prox((matrix(matrix(th, ncol = 1) + (t_s * grad), ncol = q)), penalty, t_s, lambda1);
            G <- (1/t_s) * (matrix(t(th), ncol = 1) -  prx);
            theta <- matrix(prx, byrow = TRUE, ncol=q)
            tmp <- sapply(eventtimes, sumevents_c, X, Ftime, theta, n, p, eventtimes)

            # backtracking line search
            hi <- unlist(tmp[1, ])
            lik_n <- sum(log(hi)) + sum(X[eventtimes, ] * t(theta %*% t(Ftime))) # current likelihood)
            LHS <- - lik_n
            RHS <- (lik - PEN) + (t_s  * (t(grad) %*% G)) + ( (t_s / 2)*( norm(G, "2")^2) )
            cat(LHS, ">", RHS, "\n")
            t_s <- delta * t_s
        }
        end.time <- Sys.time()
        if (verbal) cat("Time taken in backtracking:", (end.time - start.time), "\n")

        if (acceleration == TRUE) {
            theta_acc <- theta + ((tt[iter]-1)/tt[iter+1]) * (theta - theta_1) # (x_k - x_k-1)
            theta_1 <- theta # x_k-1
            tmp <- sapply(eventtimes, sumevents_c, X, Ftime, theta_acc, n, p, eventtimes)
            hi <- unlist(tmp[1, ])
            lik_n <- sum(log(hi)) + sum(X[eventtimes, ] * t(theta_acc %*% t(Ftime))) # current likelihood)
            PEN <- lambda1 * sum(apply( theta_acc, 1, FUN = function(th_p) sqrt(t(th_p) %*% penalty %*% th_p )))
        } else {
            PEN <- lambda1 * sum(apply( theta, 1, FUN = function(th_p) sqrt(t(th_p) %*% penalty %*% th_p )))
        }
        if (verbal == TRUE) cat("Penalty Term:", PEN,  "\n")

        if (-lik_n + PEN < BLIKL) {
            if (verbal) print("new best likelihood")
            BTHETA <- theta; BLIKL <- (-lik_n + PEN); BPEN <- PEN; at.iter <- iter
        }

        if (!acceleration) {
            epsilon <- max(abs(theta - th));

        } else {
            epsilon <- max(abs(theta_acc - th));

        }
        unpenlikelihood_iter[iter] <- -lik_n
        likelihood_iter[iter] <- -lik_n + PEN
        epsilon_iter[iter] <- epsilon
        if (verbal == TRUE) {
            cat("Epsilon:", epsilon, "\nIteration:", iter, "\n")
            cat("Backtracking Step Size:", t_s, "\n")
            cat("# of Backtracking Steps:", iter.back, "\n")
            cat('Neg Log Likelihood:', -lik_n, "\n")
            cat('Neg Log Likelihood w/ Penal:', -lik_n + PEN, "\n")
            cat("------------------------------ \n")
        }

        if (iter > iter.max) { stop("Likelihood did not converge after ", iter, " iterations") }
    }
    sterr <- 0
    cat("Best likelihood obs: ", BLIKL, "at iteration ", at.iter, "\n")
    if (!acceleration) {
        # Hessian calculation at final estimate:
        #tmp <- sapply(eventtimes, sumevents_reg, X, Ftime, BTHETA, n, p, eventtimes)
        #hi <- unlist(tmp[1, ])
        #lik <- -(sum(log(hi)) + sum(X[eventtimes, ] * t(BTHETA %*% t(Ftime))))
        #scores <- matrix(unlist(tmp[2, ]), ncol = k) %*% rep.int(1, k) + matrix(t(X[eventtimes, ]) %*% Ftime, ncol = 1)
        # hess <- (matrix(matrix(unlist(tmp[3, ]), ncol = k) %*% rep(1, k), ncol = p * q, byrow = FALSE))
        hess <- NULL
        outlist <- list(time = time, X = X, Ft = JJ, Ftime = Ftime, death = death, theta = theta, hess = hess,
                        btheta = BTHETA, se = sterr, lik = lik + PEN, blik = BLIKL, unpen_lik = lik_n,
                        unpen_blik = BLIKL - BPEN, iter = iter, hi = hi, PEN = penalty,
                        likelihood_iter = likelihood_iter, unpenlikelihood_iter =  unpenlikelihood_iter, epsilon_iter = epsilon_iter)
    }
    else {
        outlist <- list(time = time, X = X, Ft = JJ, Ftime = Ftime, death = death, theta = theta_acc,
                        btheta = BTHETA, se = sterr, lik = lik + PEN, blik = BLIKL, unpen_lik = lik_n,
                        unpen_blik = BLIKL - BPEN, iter = iter, hi = hi, PEN = penalty,
                        likelihood_iter = likelihood_iter, epsilon_iter = epsilon_iter)
    }

    outlist
}
