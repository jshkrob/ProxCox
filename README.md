# Prox Cox: R package

This is an R package for finding time-varying effects of covariates in high-dimensional Cox models.
Here is an example of use of the R package:

```R
library(CoxRidge)
data("GBSG")
attach(GBSG)
listg <- ProxCox::prep_data(rfst, GBSG$cens, GBSG[, c("grade", "age", "prm", "posnodal", "tumsize")], num_knots = 6, center = TRUE)
prox.fit.GBSG2 <- ProxCox::prox.fit(listg$time, listg$cens, listg$data, names = colnames(listg$data), listg$Fts, listg$knots, lambda1 = 3, lambda2 = 5, ttime = 1, iter.max = 1000, delta = 0.05, acceleration = FALSE, verbal = TRUE)
ProxCox::plot.prx(prox.fit.GBSG2, variable = 1, b.theta = FALSE)
```
