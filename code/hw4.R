rm(list=ls())
library(classInt)
library(fields)
library(maps)
library(sp)
library(gstat)
library(geoR)
library(mvtnorm)
library(MCMCpack)
library(coda)

#############################
## California temperatures ##
#############################
setwd("C:\\Users\\SAMSUNG\\Desktop\\3학년\\2학기\\베이즈통계")
load("CAtemps.RData")

ploteqc <- function(spobj, z, breaks, ...){
  pal <- tim.colors(length(breaks)-1)
  fb <- classIntervals(z, n = length(pal), 
                       style = "fixed", fixedBreaks = breaks)
  col <- findColours(fb, pal)
  plot(spobj, col = col, ...)
  image.plot(legend.only = TRUE, zlim = range(breaks), col = pal)
}

## Plotting

range(CAtemp$avgtemp)
breaks <- seq(40, 75, by = 5)
ploteqc(CAtemp, CAtemp$avgtemp, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Average Annual Temperatures, 1961-1990, Degrees F")

range(CAgrid$elevation)
breaks <- seq(-100, 3600, by = 100)
ploteqc(CAgrid, CAgrid$elevation, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Elevations at prediction locations, m")

## Prior parameters
linmod <- lm(avgtemp~lon+lat+elevation, data = CAtemp)
summary(linmod)

m.beta <- rep(0, 4); V.beta <- 1000 * diag(4)
a.s2 <- 0.001; b.s2 <- 0.001
a.t2 <- 0.001; b.t2 <- 0.001

rhoseq <- seq(0.01, 300, length = 100)
a.rho <- 2; b.rho <- 50
plot(rhoseq, dgamma(rhoseq, shape = a.rho, scale = b.rho), type = "l") # prior for rho

## Setup, storage, and starting values

y <- CAtemp$avgtemp
n <- nrow(CAtemp); m <- nrow(CAgrid)
d <- rdist.earth(coordinates(CAtemp))
X <- cbind(rep(1, n), CAtemp$lon, CAtemp$lat, CAtemp$elevation)
Xpred <- cbind(rep(1, m), CAgrid$lon, CAgrid$lat, CAgrid$elevation)

B <- 1000

beta.samps <- matrix(NA, nrow = 4, ncol = B)
beta.samps[,1] <- coef(linmod)

s2.samps <- t2.samps <- rho.samps <- rep(NA, B)
s2.samps[1] <- 5
rho.samps[1] <- 100
t2.samps[1] <- 2

eta.obs.samps <- matrix(NA, nrow = n, ncol = B)

v.prop <- 100^2

## MCMC sampler

Gamma <- exp(-d/rho.samps[1]) # initalize Gamma matrix
Ginv <- solve(Gamma)
count=0
for(i in 2:B){
  
  if(i%%100==0) print(i)
  
  ## eta_obs | Rest
  firstmu=(X%*%beta.samps)[,i-1]
  firstsigsq=s2.samps[i-1]*exp(-d/rho.samps[i-1])
  if (i==2){
    eta.obs.samps[,1]<-firstmu
  }
  secondmu=eta.obs.samps[,i-1]
  secondsigsq=t2.samps[i-1]*diag(200)
  V <- solve(solve(firstsigsq)+solve(secondsigsq))
    m <- V%*%solve(firstsigsq)%*%firstmu+V%*%solve(secondsigsq)%*%secondmu
    eta.obs.samps[,i] <- rmvnorm(1, mean = m, sigma = V, method = "svd")
  
  ## beta | Rest
  
  V <- solve(solve(V.beta)+t(X)%*%solve(exp(-d/rho.samps[i-1]))%*%X/s2.samps[i-1])
    m <- V%*%(t(X)%*%solve(exp(-d/rho.samps[i-1]))%*%eta.obs.samps[,i]/s2.samps[i-1])
    beta.samps[,i] <- rmvnorm(1, mean = m, sigma = V, method = "svd")
  
  ## s2 | Rest
  a <- a.s2+(200/2)
    b <- (b.s2+t(eta.obs.samps[,i]-(X%*%beta.samps)[,i])%*%solve(exp(-d/rho.samps[i-1]))%*%(eta.obs.samps[,i]-(X%*%beta.samps[,i])))/2
    s2.samps[i] <- rinvgamma(1,a,scale = 1/b)
  
  ## t2 | Rest
  a <- a.t2+(200/2)
    b <-  b.t2+(t(y-eta.obs.samps[,i])%*%(y-eta.obs.samps[,i]))/2
    t2.samps[i] <- rinvgamma(1, a, scale=1/b)
  
  ## rho | Rest   
  # MH update
  f<-function(rho){
    ans = dmvnorm(eta.obs.samps[,i],X%*%beta.samps[,i],s2.samps[i]*exp(-d/rho),log=TRUE)+dgamma(rho,shape=a.rho,scale=b.rho,log=TRUE)
    return (ans)
  }

  samp_rho=rnorm(1,rho.samps[i-1],20)
  if (log(runif(1))<f(samp_rho)-f(rho.samps[i-1])){
    rho.samps[i]<-samp_rho
    count=count+1
  }
  else{
    rho.samps[i]<-rho.samps[i-1]
  }

}
#acceptance rate
print(count/1000)
par(mfrow=c(2,2))

#traceplot
plot(beta.samps[1,],type='l')
plot(beta.samps[2,],type='l')
plot(beta.samps[3,],type='l')
plot(beta.samps[4,],type='l')

plot(s2.samps,type='l')
plot(t2.samps, type= 'l')
plot(rho.samps, type= 'l')
plot(eta.obs.samps[1,],type='l')

#effective sample size
library(coda)
effectiveSize(beta.samps[1,])
effectiveSize(beta.samps[2,])
effectiveSize(beta.samps[3,])
effectiveSize(beta.samps[4,])

effectiveSize(eta.obs.samps[1,])
effectiveSize(s2.samps)
effectiveSize(t2.samps)
effectiveSize(rho.samps)

#acf
acf(beta.samps[1,], lag.max=1000)
acf(beta.samps[2,], lag.max=1000)
acf(beta.samps[3,], lag.max=1000)
acf(beta.samps[4,], lag.max=1000)

acf(eta.obs.samps[1,],lag.max=1000)
acf(s2.samps,lag.max=1000)
acf(t2.samps,lag.max=1000)
acf(rho.samps,lag.max=1000)

#Hpd interval
library(HDInterval)
hdi(beta.samps[1,],0.95)
hdi(beta.samps[2,],0.95)
hdi(beta.samps[3,],0.95)
hdi(beta.samps[4,],0.95)

hdi(eta.obs.samps[1,],0.95)
hdi(s2.samps,0.95)
hdi(t2.samps,0.95)
hdi(rho.samps,0.95)

#mean
mean(beta.samps[1,])
mean(beta.samps[2,])
mean(beta.samps[3,])
mean(beta.samps[4,])

mean(eta.obs.samps[1,])
mean(s2.samps)
mean(t2.samps)
mean(rho.samps)

## Prediction

dcross <- rdist.earth(coordinates(CAtemp), coordinates(CAgrid))
dpred <- rdist.earth(coordinates(CAgrid))


eta.pred <- matrix(NA, nrow = nrow(CAgrid), ncol = B)

for(j in 1:1000){
  print(j)
  
  m <- Xpred%*%beta.samps[,j]+t(exp(-dcross/rho.samps[j]))%*%solve(exp(-d/rho.samps[j]))%*%(eta.obs.samps[,j]-X%*%beta.samps[,j])
    V <- s2.samps[j]*(exp(-dpred/rho.samps[j])-t(exp(-dcross/rho.samps[j]))%*%solve(exp(-d/rho.samps[j]))%*%exp(-dcross/rho.samps[j]))
    eta.pred[,j] <- rmvnorm(1, m, V, method = "svd")
}


saveavgtemp<-c()

par(mfrow=c(1,1))
## Plotting
for (idx in 1:664){
  saveavgtemp<-c(saveavgtemp,mean(eta.pred[idx,]))
}

range(saveavgtemp)
breaks <- seq(30, 75, by = 5)
ploteqc(CAgrid, saveavgtemp, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Average Annual Temperatures, 1961-1990, Degrees F")

savestdtemp<-c()
## Plotting
for (idx in 1:664){
  savestdtemp<-c(savestdtemp,sqrt(var(eta.pred[idx,])))
}
range(savestdtemp)
breaks <- seq(0.02, 0.2, by = 0.04)
ploteqc(CAgrid, savestdtemp, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Average Annual Temperatures Std, 1961-1990, Degrees F")
