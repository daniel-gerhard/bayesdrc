bdrm.fit <- function(x, y, model, linfct, fixed, lwr, upr, prior.mu, prior.sd, atau, btau, chains, iter, burnin, adapt, startval){
  # linfct mu fct
  mufct <- function(x, beta, model, linfct){
    b <- sapply(1:length(fixed), function(i){
      linfct[[i]] %*% beta[attr(linfct, which="lfid") == i]
    })
    return(model$fct(x, b))
  }
  
  aiter <- adapt
  p <- length(attr(linfct, which="lfid"))
  bmu <- prior.mu
  bsd <- prior.sd
  clwr <- lwr
  cupr <- upr
  loglik <- model$loglik
  ss <- model$ss
  jsd <- sqrt(bsd^2 * (2.4/sqrt(p))^2)
  # adaptive phase
  apm <- array(NA, dim=c(aiter, p, chains))
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
  
  accept <- matrix(0, nrow=chains, ncol=p)
  
  for (k in 1:chains){
    for (i in 2:aiter){
      bn <- bo <- apm[i-1,,k]
      av <- avar[i-1,k]
      for (j in 1:p){
        if (is.na(fixed[j])){
          bn[j] <- rtruncnorm(1, a=clwr[j], b=cupr[j], mean=bo[j], sd=jsd[j])
          mun <- mufct(x, bn, model, linfct)
          muo <- mufct(x, bo, model, linfct)          
          logR <- loglik(y, mun, av) - 
            loglik(y, muo, av) +
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
      muo <- mufct(x, bo, model, linfct) 
      avar[i,k] <- rgamma(1, atau + length(y)/2, btau + 0.5 * ss(y, muo))
      apm[i,,k] <- bn
    }
  }
  
  #### burnin
  pm <- array(NA, dim=c(iter, p, chains))
  var <- matrix(NA, nrow=iter, ncol=chains)
  pm[1,,] <- apm[aiter,,]
  var[1,] <- avar[aiter,]
  accept <- matrix(0, nrow=chains, ncol=p)
  
  for (k in 1:chains){
    bn <- bo <- pm[1,,k]
    av <- var[1,k]   
    for (i in 2:iter){
      for (j in 1:p){
        if (is.na(fixed[j])){
          bn[j] <- rtruncnorm(1, a=clwr[j], b=cupr[j], mean=bo[j], sd=jsd[j])
          mun <- mufct(x, bn, model, linfct)
          muo <- mufct(x, bo, model, linfct)              
          logR <- loglik(y, mun, av) - 
            loglik(y, muo, av) +
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
      muo <- mufct(x, bo, model, linfct) 
      av <- var[i,k] <- rgamma(1, atau + length(y)/2, btau + 0.5 * ss(y, muo))
      pm[i,,k] <- bn
    }
  }
  
  #### sampling
  return(pm[-(1:burnin),,])
}
