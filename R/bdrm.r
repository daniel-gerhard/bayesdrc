#' Fitting Bayesian Dose-Response Models
#' 
#' Sampling from the posterior distribution of parameters of a sigmoidal dose-response curve. \code{bdrm} is comparable to function \code{drm} in the drc package, but there are slight differences between the user-interface and of course prior distributions need to be provided.
#' 
#' @param formula an object of class \code{"\link{formula}"} specifying the name of dose and response.
#' @param data an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{bdrm} is called.
#' @param model a bdrm model function, specifying the functional dose-response relationship.
#' @param linfct either a list with a rhs \code{"\link{formula}"} for each of the parameters of the nonlinear dose-response function, conditioning on additional predictor variables, or a list of model matrices, see \code{"\link{model.matrix}"}. The default NULL assumes only an intercept for each of the parameters.
#' @param fixed a vector with fixed values for each of the model parameters. When NA a prior distribution is assumed. The default NULL specifies NA for each of the parameters.
#' @param lwr a vector with lower bounds for each parameter; NULL assumes -Inf for all parameters.
#' @param upr a vector with upper bounds for each parameter; NULL assumes Inf for all parameters. 
#' @param prior.mu a vector with prior means for each parameter. Also needs to include a value for each fixed parameter as a placeholder.
#' @param prior.sd a vector with prior standard deviations for each parameter. Also needs to include a value for each fixed parameter as a placeholder.
#' @param atau first parameter of the inverse gamma prior for the residual variance.
#' @param btau second parameter of the inverse gamma prior for the residual variance.
#' @param chains number of chains for the MCMC algorithm.
#' @param iter number of iterations for each chain.
#' @param burnin number of warm-up samples.
#' @param adapt number of samples in an adaptation phase.
#' @param startval a vector with starting values for each parameter. If NULL, the starting values are sampled from the prior distributions.  
#' 
#' @return An object of class \code{"bdrm"} is a list containing at least the following components:
#' \describe{
#'  \item{model}{a bdrm model object.}
#'  \item{linfct}{a list with model matrices.}
#'  \item{fixed}{a vector with fixed values for each of the model parameters.}
#'  \item{lwr}{a vector with lower bounds for each parameter.}
#'  \item{upr}{a vector with upper bounds for each parameter.}
#'  \item{prior.mu}{a vector with prior means for each parameter.}
#'  \item{prior.sd}{a vector with prior standard deviations for each parameter.}
#'  \item{atau}{first parameter of the inverse gamma prior for the residual variance.}
#'  \item{btau}{second parameter of the inverse gamma prior for the residual variance.}
#'  \item{chains}{number of chains.}
#'  \item{iter}{number of samples from the posterior distribution.}
#'  \item{burnin}{number of burnin samples.}
#'  \item{adapt}{number of iterations in an adaptation phase.}
#'  \item{startval}{a vector with starting values}
#'  \item{psamples}{an array with samples from the posterior distribution of the model parameters.}
#'  \item{vsamples}{a matrix with posterior samples of the residual variance.}
#'  \item{names}{a vector with names for each parameter.}
#' }
#' 
#' @seealso \code{\link[drc]{drm}}
#' 
#' @keywords models

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
  if (is.null(linfct)) linfct <- lapply(1:model$p, function(i) model.matrix(~ 1, data=data))
  if (all(unlist(lapply(linfct, function(i) inherits(i, what="formula"))))) linfct <- lapply(linfct, function(i) model.matrix(i, data=data))
  attr(linfct, which="lfid") <- rep(1:model$p, lapply(linfct, ncol))
  attr(linfct, which="name") <- paste(model$names[attr(linfct, which="lfid")], unlist(sapply(linfct, colnames)), sep=":")
  
  postsamp <- bdrmfit(x, y, model=model, linfct=linfct, fixed=fixed, lwr=lwr, upr=upr, prior.mu=prior.mu, prior.sd=prior.sd, atau=atau, btau=btau, chains=chains, iter=iter, burnin=burnin, adapt=adapt, startval=startval)
  
  dimnames(postsamp$pm)[[2]] <- attr(linfct, which="name")
  dimnames(postsamp$pm)[[3]] <- paste("chain", 1:chains, sep="")
  
  result <- list(model=model,
                 linfct=linfct,
                 fixed=fixed,
                 lwr=lwr, upr=upr,
                 prior.mu=prior.mu, prior.sd=prior.sd,
                 atau=atau, btau=btau,
                 chains=chains, iter=iter, burnin=burnin, adapt=adapt,
                 startval=startval,
                 psamples=postsamp$pm,
                 vsamples=postsamp$vm,
                 names=attr(linfct, which="name"))
  class(result) <- "bdrm"
  return(result)
}





