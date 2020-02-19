######################################################################
#######################################################################
library(truncnorm)
dose = seq(0, 1, length=25)
model <- logistic()
mus <- c(10, 15, 2, 0.5, 1)
test <- data.frame(response = rnorm(length(dose), model$fct(dose, mus), 0.1), 
                   dose=dose)
plot(response ~ dose, data=test)
startval <- c(10, 0, 2, 0.5, 1)
x <- test$dose
y <- test$response


mod <- bdrm(response ~ dose, data=test, 
            model=logistic(), 
            fixed=c(NA, NA, NA, NA, 1),
            lwr=c(0, -Inf, 0, 0, 0),
            prior.mu=c(10, 15, 2, 0.5, 1), 
            prior.sd=c(10, 10, 10, 1, 0.5), 
            atau=0.001,
            btau=0.001,
            iter=10000, burnin=8000, adapt=20000)

pm <- mod$psamples

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


library(coda)
dim(pm)
pml <- lapply(1:4, function(i) as.mcmc(pm[,,i]))
mcl <- as.mcmc.list(pml)
summary(mcl)
#plot(mcl)
gelman.diag(mcl, autoburnin=FALSE, multivariate=FALSE)




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

plot(SLOPE ~ ldose, data=spinach)

xm <- model.matrix(~ 0 + HERBICIDE, data=spinach)
mm <- model.matrix(~ 1, data=spinach)
lfct <- list(mm, mm, xm, xm, mm)


mod <- bdrm(SLOPE ~ ldose, data=spinach,
            model=logistic(), 
            linfct=lfct,
            fixed=c(NA, NA, NA, NA, NA, NA, 1),
            prior.mu=c(-10, 0, 2, 2, 0, 0, 1), 
            prior.sd=c(10, 10, 10, 10, 1, 1, 0.5), 
            upr=c(0, Inf, Inf, Inf, Inf, Inf, Inf),
            lwr=c(-Inf, 0, 0, 0, -Inf, -Inf, 0),
            atau=0.001,
            btau=0.001,
            iter=15000, burnin=10000, adapt=20000)

pm <- mod$psamples

par(mfrow=c(3,2))
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

plot(pm[,5,1], type="l", ylim=c(min(pm[,5,]), max(pm[,5,])))
lines(pm[,5,2], col="red")
lines(pm[,5,3], col="blue")
lines(pm[,5,4], col="green")

plot(pm[,6,1], type="l", ylim=c(min(pm[,6,]), max(pm[,6,])))
lines(pm[,6,2], col="red")
lines(pm[,6,3], col="blue")
lines(pm[,6,4], col="green")



library(coda)
dim(pm)
pml <- lapply(1:4, function(i) as.mcmc(pm[,,i]))
mcl <- as.mcmc.list(pml)
summary(mcl)
#plot(mcl)
gelman.diag(mcl, autoburnin=FALSE, multivariate=FALSE)

