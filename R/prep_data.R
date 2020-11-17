
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param time PARAM_DESCRIPTION
#' @param cens PARAM_DESCRIPTION
#' @param X PARAM_DESCRIPTION
#' @param num_knots PARAM_DESCRIPTION
#' @param center PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
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
    objc <- s(time, bs = "cr", k = num_knots)
    obj <- mgcv::smooth.construct.cr.smooth.spec(objc, list(time = time), knots = list( time = knots))
    li <- list(time = time, cens = cens, data = X, Fts = obj$X, knots = knots)
    return(li)
}
