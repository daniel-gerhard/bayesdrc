bdrm <- function(formula, data, model, linfct=NULL, fixed=NULL, lwr=NULL, upr=NULL, prior.mu, prior.sd, atau=0.001, btau=0.001, chains=4, iter=10000, burnin=9000, adapt=10000, startval=NULL){
  UseMethod("bdrm")
}

bdrm.formula <- function(formula, data, model, linfct=NULL, fixed=NULL, lwr=NULL, upr=NULL, prior.mu, prior.sd, atau=0.001, btau=0.001, chains=4, iter=10000, burnin=9000, adapt=10000, startval=NULL){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)[,-1]
  
  if (is.null(fixed)) fixed <- rep(NA, model$p)
  if (is.null(lwr)) lwr <- rep(-Inf, model$p)
  if (is.null(upr)) upr <- rep(Inf, model$p)
  
  # linfct
  if (is.null(linfct)) linfct <- lapply(1:model$p, function(i) matrix(1))
  attr(linfct, which="lfid") <- rep(1:model$p, lapply(linfct, ncol))
  
  pm <- bdrm.fit(x, y, model=model, linfct=linfct, fixed=fixed, lwr=lwr, upr=upr, prior.mu=prior.mu, prior.sd=prior.sd, atau=atau, btau=btau, chains=chains, iter=iter, burnin=burnin, adapt=adapt, startval=startval)
  return(pm)
}





