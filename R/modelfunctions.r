logistic <- function(names=c("b1", "b2", "b3", "b4", "b5")){
  fct <- function(x, beta){
    beta[2] + beta[3]/(1 + exp(-beta[1]*(x-beta[4]))^beta[5])  
  }
  loglik <- function(y, mu, variance){
    LL <- sum(dnorm(y, mu, 1/sqrt(variance), log=TRUE))
    return(LL)
  }
  ss <- function(y, mu){
    SS <- sum((y-mu)^2)
    return(SS)
  }
  mod <- list(fct=fct, 
              loglik=loglik, 
              ss=ss, 
              p=5,
              names=names)
  return(mod)
}


