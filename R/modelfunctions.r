#' Logistic Dose-Response Model
#' 
#' Assuming a logistic dose-response relationship in a bdrm model.
#' 
#' @param names a vector with names for each model parameter
#' 
#' @return a list containing at least the following components:
#' \describe{
#'   \item{fct}{the logistic function with a dose x and a vector of parameters beta}
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
  mod <- list(fct=fct, 
              p=5,
              names=names)
  return(mod)
}


#' Weibull Dose-Response Model
#' 
#' Assuming a Weibull dose-response relationship in a bdrm model.
#' 
#' @param names a vector with names for each model parameter
#' 
#' @return a list containing at least the following components:
#' \describe{
#'   \item{fct}{the logistic function with a dose x and a vector of parameters beta}
#'   \item{p}{number of model parameters}
#'   \item{names}{a vector with parameter names}
#' }
#' 
#' @details 
#' The logistic function is defined as 
#' \deqn{f(x, \beta) = \beta_2 + \beta_3 \exp(- \exp( -\beta_1 (\log(x)-\log(\beta_4)))) }
#' there exists a second parameterisation that is implemented in function \code{weibull2}
#' \deqn{f(x, \beta) = \beta_2 + \beta_3 (1 - \exp(- \exp( \beta_1 (\log(x)-\log(\beta_4))))) }
#' 
#' @keywords models

weibull1 <- function(names=c("b1", "b2", "b3", "b4")){
  fct <- function(x, beta){
    beta[2] + beta[3]*exp(-exp(-beta[1]*(log(x)-log(beta[4]))))  
  }
  mod <- list(fct=fct, 
              p=4,
              names=names)
  return(mod)
}

#' @rdname weibull1
weibull2 <- function(names=c("b1", "b2", "b3", "b4")){
  fct <- function(x, beta){
    beta[2] + beta[3]*(1 - exp(-exp(beta[1]*(log(x)-log(beta[4])))))  
  }
  mod <- list(fct=fct, 
              p=4,
              names=names)
  return(mod)
}

#' Asymptotic Regression Model
#' 
#' Assuming an asymptotic regression model in bdrm
#' 
#' @param names a vector with names for each model parameter
#' 
#' @return a list containing at least the following components:
#' \describe{
#'   \item{fct}{the logistic function with a dose x and a vector of parameters beta}
#'   \item{p}{number of model parameters}
#'   \item{names}{a vector with parameter names}
#' }
#' 
#' @details 
#' The asymptotic regression model is defined as 
#' \deqn{f(x, \beta) = \beta_1 + \beta_2 (1 - \exp(\frac{-x}{\beta_3})) }
#' 
#' @keywords models

asyreg <- function(names=c("b1", "b2", "b3")){
  fct <- function(x, beta){
    beta[1] + beta[2]*(1 - exp(-x/beta[3]))  
  }
  mod <- list(fct=fct, 
              p=3,
              names=names)
  return(mod) 
}


#' Log-Normal Dose-Response Model
#' 
#' Assuming a log-normal dose-response relationship in a bdrm model.
#' 
#' @param names a vector with names for each model parameter
#' 
#' @return a list containing at least the following components:
#' \describe{
#'   \item{fct}{the logistic function with a dose x and a vector of parameters beta}
#'   \item{p}{number of model parameters}
#'   \item{names}{a vector with parameter names}
#' }
#' 
#' @details 
#' The log-normal model is defined as 
#' \deqn{f(x, \beta) = \beta_2 + \beta_3 \Phi( \beta_1 (\log(x)-\beta_4)) }
#' 
#' @keywords models

lognormal <- function(names=c("b1", "b2", "b3", "b4")){
  fct <- function(x, beta){
    beta[2] + beta[3]*pnorm(beta[1]*(log(x)-beta[4]))  
  }
  mod <- list(fct=fct, 
              p=4,
              names=names)
  return(mod)
}
