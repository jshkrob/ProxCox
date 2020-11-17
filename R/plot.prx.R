#' @export
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
