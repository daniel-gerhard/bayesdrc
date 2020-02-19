traceplot <- function(x){
  pm <- x$psamples     
  pmd <- as.data.frame(apply(pm, 2, function(x) x))
  pmd$iteration <- rep(1:(x$iter-x$burnin), times=x$chains)
  pmd$chain <- factor(rep(1:x$chains, each=x$iter-x$burnin))
  pmdl <- pmd %>% gather(key="parameter", value="sample", -chain, -iteration)
  ggplot(pmdl, aes(x=iteration, y=sample, colour=chain)) +
    geom_line() +
    facet_wrap(~ parameter, scales="free_y")
}