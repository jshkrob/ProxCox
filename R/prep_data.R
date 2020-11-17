
#' @title prep_data
#' @description outputs list of processes survival data necissary for prox.fit and reg.fit
#' @param time survival time (either inorder or unordered)
#' @param cens censored indicator (binary vector)
#' @param X covariates of interest
#' @param num_knots number of knots to use in spline
#' @param center centering, Default: FALSE
#' @return list containg processed data (time, cens, X) and basis functions evaluated at each time point, knots used
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[mgcv]{smooth.construct.cr.smooth.spec}}
#' @rdname prep_data
#' @export
#' @importFrom mgcv smooth.construct.cr.smooth.spec
#' @importFrom mgcv s
prep_data <- function(time, cens, X, num_knots, center = FALSE) {
    cens = cens[order(time)]
    time <- time[order(time)]
    names.X <- colnames(X)
    if (center){
        X <- X[order(time), ]
        X <- apply(X, 2, function(x) {(x - mean(x)) / sqrt(var(x))})
        #D <- cbind(time = time, cens = cens, apply(X, 2, function(x) {(x - mean(x)) / sqrt(var(x))}))
    } else {
        X <- X[order(time), ]
        #D <- cbind(time = time, cens = cens, X = X)
    }
    colnames(X) <- c(names.X)
    knots <- quantile(time, seq(0.0, 1.0, 1/(num_knots-1)))
    objc <- mgcv::s(time, bs = "cr", k = num_knots)
    obj <- mgcv::smooth.construct.cr.smooth.spec(objc, list(time = time), knots = list( time = knots))
    li <- list(time = time, cens = cens, data = X, Fts = obj$X, knots = knots)
    return(li)
}
