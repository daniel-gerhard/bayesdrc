#' traceplot for bdrm objects
#' 
#' @param x an object of class bdrm.
#' @param which index of the parameter of interest; default are all parameters.
#' @param include_burnin should samples from the warmup phase be included.
#' 
#' @returns a ggplot2 plot.

traceplot <- function(x, which=1:length(x$names), include_burnin=FALSE){
  pmdl <- postsamples(x, which=which, include_burnin=include_burnin)
  if (include_burnin == TRUE){
    ggplot(pmdl, aes(x=iteration, y=sample, colour=chain)) +
      geom_line(alpha=0.75) +
      facet_wrap(~ parameter, scales="free_y") +
      geom_vline(xintercept=x$burnin, linetype=2)
  } else {
    ggplot(pmdl, aes(x=iteration, y=sample, colour=chain)) +
      geom_line(alpha=0.75) +
      facet_wrap(~ parameter, scales="free_y")
  }
}