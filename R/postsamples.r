#' Extract posterior samples from a bdrm object
#' 
#' @param x an object of class bdrm.
#' @param which index of the parameter of interest; default are all parameters.
#' @param include_burnin should samples from the warmup phase be included.
#' 
#' @return a data.frame containing samples from the posterior distribution with columns: chain, iteration, parameter and sample.


postsamples <- function(x, which=1:length(x$names), include_burnin=FALSE){
  if (include_burnin == TRUE){
    psa <- x$psamples[,which,,drop=FALSE]
    it <- x$iter
  } else {
    psa <- x$psamples[-c(1:x$burnin),which,,drop=FALSE]
    it <- x$iter-x$burnin
  }
  pmd <- as.data.frame(apply(psa, 2, function(x) x))
  pmd$iteration <- rep(1:it, times=x$chains)
  pmd$chain <- factor(rep(1:x$chains, each=it))
  pmdl <- pmd %>% gather(key="parameter", value="sample", -chain, -iteration)
  return(pmdl)
}