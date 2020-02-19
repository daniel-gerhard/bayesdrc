print.bdrm <- function(x, ...){
  est <- apply(x$psamples, 2, mean)
  sds <- apply(x$psamples, 2, sd)
  quants <- t(apply(x$psamples, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))
  res <- data.frame(est, sds, quants)
  names(res) <- c("Mean", "SD", "2.5%", "25%", "50%", "75%", "97.5%")
  print(res, digits=3)
}