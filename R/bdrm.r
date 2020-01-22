bdrm <- function(formula, data, family){
  UseMethod("bdrm")
}

bdrm.formula <- function(formula, data, family, model){
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
  return(list(y, x))
}


test <- data.frame(response = c(1, 2, 3, 4), dose=c(4, 3, 2, 1))
bdrm(response ~ dose, data=test, family="gaussian")


logistic <- function(fixed=c(NA, NA, NA, NA, NA), 
                     prior.mu=c(),
                     prior.sd=c(),
                     lwr=c(-Inf, -Inf, -Inf, -Inf, -Inf),
                     upr=c(Inf, Inf, Inf, Inf, Inf)){
  
}
