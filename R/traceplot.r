#' traceplot for bdrm objects
#' 
#' @param x an object of class bdrm.
#' @param which index of the parameter of interest; default are all parameters.
#' 
#' @returns a ggplot2 plot.

traceplot <- function(x, which=1:length(x$names)){
  pmdl <- postsamples(x, which=which)
  ggplot(pmdl, aes(x=iteration, y=sample, colour=chain)) +
    geom_line(alpha=0.75) +
    facet_wrap(~ parameter, scales="free_y")
}