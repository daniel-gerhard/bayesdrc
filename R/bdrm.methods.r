print.bdrm <- function(x, ...){
  # rhat
  m <- dim(x$psamples)[3]
  n <- dim(x$psamples)[1]
  B <- n/(m-1) * apply((apply(x$psamples, c(2,3), mean) - apply(x$psamples, 2, mean))^2, 1, sum)
  W <- apply(apply(x$psamples, c(2,3), var), 1, mean)
  vpy <- (n-1)/n * W + B/n
  rhat <- sqrt(vpy/W)
  
  # summary stats
  est <- apply(x$psamples, 2, mean)
  sds <- apply(x$psamples, 2, sd)
  quants <- t(apply(x$psamples, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))
  
  # print
  res <- data.frame(est, sds, quants, rhat)
  names(res) <- c("Mean", "SD", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat")
  print(res, digits=3)
}