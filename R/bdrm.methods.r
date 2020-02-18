print.bdrm <- function(x, ...){
  est <- apply(mod$psamples, 2, mean)
  sds <- apply(mod$psamples, 2, sd)
  quants <- t(apply(mod$psamples, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))
  rdat <- data.frame(est, sds, quants)
  names(rdat) <- c("Mean", "SD", "2.5%", "25%", "50%", "75%", "97.5%")
  print(rdat, digits=3)
}