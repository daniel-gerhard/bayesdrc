bdrm <- function(formula, data, model, iter){
  UseMethod("bdrm")
}

bdrm.formula <- function(formula, data, model, iter){
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


test <- data.frame(response = c(1, 2, 3, 4), 
                   dose = c(4, 3, 2, 1))
bdrm(response ~ dose, data=test)


logistic <- function(fixed=c(NA, NA, NA, NA, NA), 
                     prior.mu=c(),
                     prior.sd=c(),
                     lwr=c(-Inf, -Inf, -Inf, -Inf, -Inf),
                     upr=c(Inf, Inf, Inf, Inf, Inf)){

  fct <- function(x, beta) beta[2] + (beta[3]-beta[2])/(1 + exp(-beta[1]*(x-beta[4]))^beta[5])  
  loglik <- function(x, y, beta, variance){
    mu <- fct(x, beta)
    LL <- sum(dnorm(y, mu, 1/sqrt(variance), log=TRUE))
    return(LL)
  }
  ss <- function(x, y, beta){
    mu <- richards(x, beta)
    SS <- sum((y-mu)^2)
    return(SS)
  }
  mod <- list(fct, prior.mu, prior.sd, lwr, upr, loglik, ss)
  return(mod)
}


bdrm.fit <- function(x, y, model, iter){
  p <- length(model$fixed)
  bmu <- model$prior.mu
  bsd <- model$prior.sd
  clwr <- model$lwr
  cupr <- model$upr
  loglik <- model$loglik
  for (i in 1:iter){
    for (j in 1:p){
      bn[j] <- rtruncnorm(1, a=clwr[j], b=cupr[j], mean=bo[j], sd=0.5)
      logR <- loglik(x, y, bn) - 
        loglik(x, y, bo) +
        log(dtruncnorm(bn[j], a=clwr[j], b=cupr[j], mean=bmu[j], sd=bsd[j])) - 
        log(dtruncnorm(bo[j], a=clwr[j], b=cupr[j], mean=bmu[j], sd=bsd[j])) -
        log(dtruncnorm(bn[j], a=clwr[j], b=cupr[j], mean=bo[j], sd=0.5)) +
        log(dtruncnorm(bo[j], a=clwr[j], b=cupr[j], mean=bn[j], sd=0.5))
      logU <- log(runif(1, 0, 1))
      if (logR > logU){
        bo <- bn
        #accept[1] <- accept[1] + 1
      } else {
        bn <- bo
      }
    }
  }
}
