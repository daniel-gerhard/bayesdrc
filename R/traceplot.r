traceplot <- function(x, which=1:length(x$names)){
  pmdl <- postsamples(x, which=which)
  ggplot(pmdl, aes(x=iteration, y=sample, colour=chain)) +
    geom_line() +
    facet_wrap(~ parameter, scales="free_y")
}