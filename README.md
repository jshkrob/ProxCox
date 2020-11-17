# Prox Cox: R package

This is an R package for finding time-varying effects of covariates in high-dimensional Cox models.
Here is an example of use of the R package:

```R
library(ProxCox)
D <- cbind(time = rfst, cens = GBSG$cens, apply(GBSG[, c("age", "grade", "prm", "posnodal", "tumsize")], 2, function(x){(x-mean(x))/sqrt(var(x))}))
knots_GBSG <- quantile(rfst, seq(0.0, 1.0, 0.2))
obj_GBSG <- s(rfst, bs = "cr", k = 6, m = 2)
objc_GBSG <- smooth.construct.cr.smooth.spec(obj_GBSG, list(rfst = rfst), knots = list( rfst = knots_GBSG))
Fts.GBSG <- cbind(rep(1, length(rfst)), objc_GBSG$X) # attatched const term
prox.fit.GBSG <- ProxCox::prox.fit(rfst, cens, D[, c("age", "grade", "prm", "posnodal", "tumsize")], names<-c("age", "grade", "prm", "posnodal", "tumsize"),
                          objc_GBSG$X, knots_GBSG, lambda1 = 10, lambda2 = 5, ttime = 0.5, iter.max = 20, delta = 0.05, acceleration = FALSE, verbal = TRUE)
prox.fit.GBSG2 <- prox.fit(rfst, cens, D[, c("grade", "age", "prm", "posnodal", "tumsize")], names<-c("grade", "age", "prm", "posnodal", "tumsize"),
                          objc_GBSG$X, knots_GBSG, lambda1 = 1, lambda2 = 5, ttime = 1, iter.max = 1000, delta = 0.05, acceleration = FALSE, verbal = TRUE)

reg.fit.GBSG <- reg.fit(time = rfst, death = cens, X = GBSG[, c("age", "grade", "prm", "posnodal", "tumsize")], Ftime = Fts.GBSG, t_s = 1, iter.max = 1000)

ProxCox::plot.prx(prox.fit.GBSG2, variable = 5, b.theta = FALSE)
```
