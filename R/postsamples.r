postsamples <- function(x, which=1:length(x$names)){
  pmd <- as.data.frame(apply(x$psamples[,which,,drop=FALSE], 2, function(x) x))
  pmd$iteration <- rep(1:(x$iter-x$burnin), times=x$chains)
  pmd$chain <- factor(rep(1:x$chains, each=x$iter-x$burnin))
  pmdl <- pmd %>% gather(key="parameter", value="sample", -chain, -iteration)
  return(pmdl)
}