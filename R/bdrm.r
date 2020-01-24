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


logistic <- function(fixed=c(NA, NA, NA, NA, NA), 
                     prior.mu=c(10, 0, 1, 0.5, 1),
                     prior.sd=c(1, 1, 1, 1, 1),
                     lwr=c(-Inf, -Inf, -Inf, -Inf, -Inf),
                     upr=c(Inf, Inf, Inf, Inf, Inf)){

  fct <- function(x, beta) beta[2] + (beta[3]-beta[2])/(1 + exp(-beta[1]*(x-beta[4]))^beta[5])  
  loglik <- function(x, y, beta, variance){
    mu <- fct(x, beta)
    LL <- sum(dnorm(y, mu, 1/sqrt(variance), log=TRUE))
    return(LL)
  }
  ss <- function(x, y, beta){
    mu <- fct(x, beta)
    SS <- sum((y-mu)^2)
    return(SS)
  }
  mod <- list(fct=fct, fixed=fixed, prior.mu=prior.mu, prior.sd=prior.sd, lwr=lwr, upr=upr, loglik=loglik, ss=ss)
  return(mod)
}

library(truncnorm)
dose = seq(0, 1, length=25)
model <- logistic()
test <- data.frame(response = rnorm(length(dose), model$fct(dose, model$prior.mu), 0.1), 
                   dose=dose)
plot(response ~ dose, data=test)
bdrm(response ~ dose, data=test)
iter <- 100
startval <- c(10, 0, 1, 0.5, 1)
x <- test$dose
y <- test$response

chains <- 4

bdrm.fit <- function(x, y, model, chains, iter, startval){
  p <- length(model$fixed)
  bmu <- model$prior.mu
  bsd <- model$prior.sd
  clwr <- model$lwr
  cupr <- model$upr
  loglik <- model$loglik
  ss <- model$ss
  jsd <- sqrt(bsd^2 * (2.4/sqrt(p))^2)
  # adaptive phase
  apm <- array(NA, dim=c(iter, length(startval), chains))
  avar <- matrix(NA, nrow=iter, ncol=chains)
  apm[1,,] <- startval
  avar[1,] <- 1
  atau <- 0.001
  btau <- 0.001
  
  accept <- matrix(0, nrow=chains, ncol=p)
  
  for (k in 1:chains){
    for (i in 2:iter){
      bn <- bo <- apm[i-1,,k]
      av <- avar[i-1,k]
      for (j in 1:p){
        bn[j] <- rtruncnorm(1, a=clwr[j], b=cupr[j], mean=bo[j], sd=jsd[j])
        logR <- loglik(x, y, bn, av) - 
          loglik(x, y, bo, av) +
          log(dtruncnorm(bn[j], a=clwr[j], b=cupr[j], mean=bmu[j], sd=bsd[j])) - 
          log(dtruncnorm(bo[j], a=clwr[j], b=cupr[j], mean=bmu[j], sd=bsd[j])) -
          log(dtruncnorm(bn[j], a=clwr[j], b=cupr[j], mean=bo[j], sd=jsd[j])) +
          log(dtruncnorm(bo[j], a=clwr[j], b=cupr[j], mean=bn[j], sd=jsd[j]))
        logU <- log(runif(1, 0, 1))
        if (logR > logU){
          bo <- bn
          accept[k,j] <- accept[k,j] + 1
        } else {
          bn <- bo
        }
      }
      avar[i,k] <- rgamma(1, atau + length(y)/2, btau + 0.5 * ss(x, y, bo))
      apm[i,,k] <- bn
    }
  }
  #### burnin
  #### sampling
}
