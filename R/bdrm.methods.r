print.bdrm <- function(x, ...){
  # remove burnin
  pm <- x$psamples[-(1:x$burnin),,,drop=FALSE]
  # rhat
  m <- dim(pm)[3]
  n <- dim(pm)[1]
  B <- n/(m-1) * apply((apply(pm, c(2,3), mean) - apply(pm, 2, mean))^2, 1, sum)
  W <- apply(apply(pm, c(2,3), var), 1, mean)
  vpy <- (n-1)/n * W + B/n
  rhat <- sqrt(vpy/W)
  
  # summary stats
  est <- apply(pm, 2, mean)
  sds <- apply(pm, 2, sd)
  quants <- t(apply(pm, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))
  
  # print
  res <- data.frame(est, sds, quants, rhat)
  names(res) <- c("Mean", "SD", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat")
  print(res, digits=3)
}



