######################################################################
#######################################################################
library(truncnorm)
library(ggplot2)
library(tidyr)
#
library(bayesdrc)

##################################

dose = seq(0, 1, length=25)
model <- logistic()
mus <- c(10, 15, 2, 0.5, 1)
test <- data.frame(response = rnorm(length(dose), model$fct(dose, mus), 0.1), 
                   dose=dose)
plot(response ~ dose, data=test)
startval <- c(10, 0, 2, 0.5, 1)

###
dose = seq(0, 1, length=25)
model <- weibull1()
mus <- c(10, 15, 2, 0.5)
w1 <- function(x, beta) beta[2] + beta[3]*exp(-exp(-beta[1]*(log(x)-log(beta[4])))) 
test <- data.frame(response = rnorm(length(dose), w1(dose, mus), 0.1), 
                   dose=dose)
plot(response ~ dose, data=test)
startval <- c(10, 0, 2, 0.5)


###############3
mod1 <- bdrm(response ~ dose, data=test, 
             model=logistic(), 
             fixed=c(NA, NA, NA, NA, 1),
             lwr=c(0, -Inf, 0, 0, 0),
             prior.mu=c(10, 15, 2, 0.5, 1), 
             prior.sd=c(10, 10, 10, 1, 0.5), 
             atau=0.001,
             btau=0.001,
             iter=10000, burnin=8000, adapt=20000)



mod2 <- bdrm(response ~ dose, data=test, 
             model=weibull1(), 
             fixed=c(NA, NA, NA, NA),
             lwr=c(0, -Inf, 0, 0),
             prior.mu=c(10, 15, 2, 0.5), 
             prior.sd=c(10, 10, 10, 1), 
             atau=0.001,
             btau=0.001,
             iter=10000, burnin=8000, adapt=20000)



################################################


library(drcData)
data(ryegrass)
y <- ryegrass$rootl
x <- log(ryegrass$conc + 0.25)
plot(x, y)

rt <- data.frame(x, y)

mod <- bdrm(y ~ x, data=rt, 
            model=logistic(), 
            fixed=c(NA, NA, NA, NA, NA),
            prior.mu=c(-10, 8, 0, 1, 1), 
            prior.sd=c(10, 10, 10, 1, 0.5), 
            upr=c(0, Inf, Inf, Inf, Inf),
            lwr=c(-Inf, 0, 0, -Inf, 0),
            atau=0.001,
            btau=0.001,
            iter=15000, burnin=10000, adapt=20000)

pm <- mod$psamples

library(coda)
dim(pm)
pml <- lapply(1:4, function(i) as.mcmc(pm[,,i]))
mcl <- as.mcmc.list(pml)
summary(mcl)
#plot(mcl)
gelman.diag(mcl, autoburnin=FALSE, multivariate=FALSE)


xc <- seq(-2, 4, length=100)
yc <- logistic()$fct(xc, summary(mcl)$statistics[,1])
par(mfrow=c(1,1))
plot(x, y)

pr <- apply(pm, c(1, 3), function(b) logistic()$fct(xc, b))
qprm <- apply(pr, 1, quantile, probs=0.5)
qpr1 <- apply(pr, 1, quantile, probs=0.025)
qpr2 <- apply(pr, 1, quantile, probs=0.975)

lines(xc, qprm)
lines(xc, qpr1, lty=2)
lines(xc, qpr2, lty=2)


####################################

library(drcData)
data(spinach)
spinach$ldose <- log(spinach$DOSE + 0.005)

mod <- bdrm(SLOPE ~ ldose, data=spinach,
            model=logistic(), 
            linfct=c(~ 0 + HERBICIDE, ~1, ~ 0 + HERBICIDE, ~ 0 + HERBICIDE, ~ 0 + HERBICIDE),
            fixed=c(NA, NA, 0, NA, NA, NA, NA, NA, NA),
            prior.mu=c(-10, -10, 0, 2, 2, 0, 0, 1, 1), 
            prior.sd=c(10, 10, 10, 10, 10, 1, 1, 0.5, 0.5), 
            upr=c(0, 0, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
            lwr=c(-Inf, -Inf, 0, 0, 0, -Inf, -Inf, 0, 0),
            atau=0.001,
            btau=0.001,
            iter=15000, burnin=10000, adapt=20000)

traceplot(mod, which=1:9)


xc <- seq(-5, 5, length=100)

pra <- apply(mod$psamples, c(1, 3), function(b) logistic()$fct(xc, b[c(1, 3, 4, 6, 8)]))
qprma <- apply(pra, 1, quantile, probs=0.5)
qpr1a <- apply(pra, 1, quantile, probs=0.025)
qpr2a <- apply(pra, 1, quantile, probs=0.975)

prb <- apply(mod$psamples, c(1, 3), function(b) logistic()$fct(xc, b[c(2, 3, 5, 7, 9)]))
qprmb <- apply(prb, 1, quantile, probs=0.5)
qpr1b <- apply(prb, 1, quantile, probs=0.025)
qpr2b <- apply(prb, 1, quantile, probs=0.975)

plot(SLOPE ~ ldose, data=spinach, col=c("red2", "blue2")[spinach$HERBICIDE])
lines(xc, qprma, col="red2")
lines(xc, qpr1a, lty=2, col="red2")
lines(xc, qpr2a, lty=2, col="red2")
lines(xc, qprmb, col="blue2")
lines(xc, qpr1b, lty=2, col="blue2")
lines(xc, qpr2b, lty=2, col="blue2")
