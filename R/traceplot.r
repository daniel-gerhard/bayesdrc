traceplot <- function(x, which=1:length(x$names)){
  pmdl <- postsamples(x, which=which)
  ggplot(pmdl, aes(x=iteration, y=sample, colour=chain)) +
    geom_line(alpha=0.75) +
    facet_wrap(~ parameter, scales="free_y")
}