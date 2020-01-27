bdrm <- function(formula, data, model, chains, iter, startval=NULL){
  UseMethod("bdrm")
}

bdrm.formula <- function(formula, data, model, chains, iter, startval=NULL){
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
  pm <- bdrm.fit(x, y, model=model, chains=chains, iter=iter, startval=startval)
  return(pm)
}


logistic <- function(fixed=c(NA, NA, NA, NA, NA), 
                     prior.mu=c(NA, NA, NA, NA, NA),
                     prior.sd=c(NA, NA, NA, NA, NA),
                     lwr=c(-Inf, -Inf, -Inf, -Inf, -Inf),
                     upr=c(Inf, Inf, Inf, Inf, Inf)){

  fct <- function(x, beta) beta[2] + (beta[3])/(1 + exp(-beta[1]*(x-beta[4]))^beta[5])  
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


bdrm.fit <- function(x, y, model, chains, iter, startval){
  fixed <- model$fixed
  p <- length(fixed)
  bmu <- model$prior.mu
  bsd <- model$prior.sd
  clwr <- model$lwr
  cupr <- model$upr
  loglik <- model$loglik
  ss <- model$ss
  jsd <- sqrt(bsd^2 * (2.4/sqrt(p))^2)
  # adaptive phase
  apm <- array(NA, dim=c(aiter, length(startval), chains))
  avar <- matrix(NA, nrow=aiter, ncol=chains)
  if (is.null(startval)) startval <- matrix(rtruncnorm(p*chains, a=clwr, b=cupr, mean=bmu, sd=bsd), nrow=p, ncol=chains)
  startval <- matrix(startval, nrow=p, ncol=chains)
  if (any(!is.na(fixed))){
    wst <- which(!is.na(fixed))
    for (w in wst){
      startval[w,] <- fixed[w] 
    }
  }
  apm[1,,] <- startval
  avar[1,] <- 1
  atau <- 0.001
  btau <- 0.001
  
  accept <- matrix(0, nrow=chains, ncol=p)
  
  for (k in 1:chains){
    for (i in 2:aiter){
      bn <- bo <- apm[i-1,,k]
      av <- avar[i-1,k]
      for (j in 1:p){
        if (is.na(fixed[j])){
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
          if (i == 200){
            vest <- sqrt(var(apm[50:199,j,k]) * (2.4/sqrt(p))^2)
            if (vest == 0) vest <- jsd[j] * 0.5
            jsd[j] <- vest
          }
          if (i %in% seq(300, aiter, by=100)){ 
            if (accept[k,j]/100 < 0.4) jsd[j] <- jsd[j]*0.95
            if (accept[k,j]/100 > 0.55) jsd[j] <- jsd[j]*1.05
            accept[k,j] <- 0
          }
        } 
      }
      avar[i,k] <- rgamma(1, atau + length(y)/2, btau + 0.5 * ss(x, y, bo))
      apm[i,,k] <- bn
    }
  }

  #### burnin
  pm <- array(NA, dim=c(iter, length(startval), chains))
  var <- matrix(NA, nrow=iter, ncol=chains)
  pm[1,,] <- apm[aiter,,]
  var[1,] <- avar[aiter,]
  accept <- matrix(0, nrow=chains, ncol=p)
  
  for (k in 1:chains){
    for (i in 2:iter){
      bn <- bo <- pm[i-1,,k]
      av <- var[i-1,k]
      for (j in 1:p){
        if (is.na(fixed[j])){
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
      }
      var[i,k] <- rgamma(1, atau + length(y)/2, btau + 0.5 * ss(x, y, bo))
      pm[i,,k] <- bn
    }
  }

  #### sampling
  return(pm)
}


######################################################################
#######################################################################
library(truncnorm)
dose = seq(0, 1, length=25)
model <- logistic()
mus <- c(10, 15, 2, 0.5, 1)
test <- data.frame(response = rnorm(length(dose), model$fct(dose, mus), 0.1), 
                   dose=dose)
plot(response ~ dose, data=test)
bdrm(response ~ dose, data=test)
aiter <- 10000
iter <- 10000
startval <- c(10, 0, 2, 0.5, 1)
x <- test$dose
y <- test$response

chains <- 4

pm <- bdrm(response ~ dose, data=test, 
           model=logistic(prior.mu=c(10, 15, 2, 0.5, 1), 
                          prior.sd=c(10, 10, 10, 10, 10), 
                          lwr=c(0, -Inf, 0, 0, 0), 
                          fixed=c(NA, NA, NA, NA, 1)), 
           chains=4, 
           iter=10000, 
           startval=c(1, 0, 2, 0.5, 1))

par(mfrow=c(2,2))
plot(pm[,1,1], type="l", ylim=c(min(pm[,1,]), max(pm[,1,])))
lines(pm[,1,2], col="red")
lines(pm[,1,3], col="blue")
lines(pm[,1,4], col="green")

plot(pm[,2,1], type="l", ylim=c(min(pm[,2,]), max(pm[,2,])))
lines(pm[,2,2], col="red")
lines(pm[,2,3], col="blue")
lines(pm[,2,4], col="green")

plot(pm[,3,1], type="l", ylim=c(min(pm[,3,]), max(pm[,3,])))
lines(pm[,3,2], col="red")
lines(pm[,3,3], col="blue")
lines(pm[,3,4], col="green")

plot(pm[,4,1], type="l", ylim=c(min(pm[,4,]), max(pm[,4,])))
lines(pm[,4,2], col="red")
lines(pm[,4,3], col="blue")
lines(pm[,4,4], col="green")


pairs(pm[,1:4,1], pch=15)
