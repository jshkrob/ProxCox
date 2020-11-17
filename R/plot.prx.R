#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param obj PARAM_DESCRIPTION
#' @param alpha PARAM_DESCRIPTION, Default: 0.05
#' @param xlab PARAM_DESCRIPTION, Default: 'time'
#' @param ylab PARAM_DESCRIPTION, Default: 'X effect'
#' @param all.terms PARAM_DESCRIPTION, Default: TRUE
#' @param variable PARAM_DESCRIPTION
#' @param continue PARAM_DESCRIPTION, Default: FALSE
#' @param col PARAM_DESCRIPTION
#' @param b.theta PARAM_DESCRIPTION, Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname plot.prx
#' @export plot.prx
plot.prx <- function(obj,alpha=0.05,xlab="time",ylab="X effect",all.terms=TRUE,variable, continue = FALSE, col, b.theta=TRUE) {
    dm <- ncol(obj$X) * nrow(obj$X) # p*q
    sq <- seq(1,dm,by=ncol(obj$Ft)) # q = nubmer of basis functions (1, q+1, 2*q + 1, )
    dim <- variable
    if (b.theta == TRUE) {f2 <- obj$btheta[dim,]%*%t(obj$Ft)}
    else{         f2 <- obj$theta[dim,]%*%t(obj$Ft)}
    if (continue == TRUE) {
        lines(obj$time,f2,lwd=3, col = col)
    } else {
        plot(obj$time,f2,'l',xlab=xlab,ylab=ylab,lwd=3, ylim = c(-3, 3), xlim = c(0, max(obj$time)))
        #plot(obj$time,f2,'l',xlab=xlab,ylab=ylab,lwd=3, col=col)
        #lines(obj$time,f2,lwd=3,lty=2)
    }
}
