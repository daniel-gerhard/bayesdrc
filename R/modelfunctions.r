#' Logistic Dose-Response Model
#' 
#' Assuming a logistic dose-response relationship in a bdrm model.
#' 
#' @param names a vector with names for each model parameter
#' 
#' @return a list containing at least the following components:
#' \describe{
#'   \item{fct}{the logistic function with a dose x and a vector of parameters beta}
#'   \item{loglik}{the log-likelihood function}
#'   \item{ss}{sum-of-squares for observed and fitted values}
#'   \item{p}{number of model parameters}
#'   \item{names}{a vector with parameter names}
#' }
#' 
#' @details 
#' The logistic function is defined as 
#' \deqn{f(x, \beta) = \beta_2 + \frac{\beta_3}{(1 + \exp(-\beta_1 (x-\beta_4))^{\beta_5})}}
#' 
#' @keywords models

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


