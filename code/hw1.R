a=5
b=6
x=seq(0,100,1)
#prior
plot(dgamma(x,shape=a,scale=b),type='l')

#likelihood
xx=seq(0,100,1)
data=52
plot(dpois(data,lambda=xx),type="l")

#posterior
plot(dgamma(x,shape=data+a,scale=b/(b+1)),type='l')

library(optim.functions)
fn<-function(theta){
  ans=theta^56*exp((-7/6)*theta)
  return(-ans)
}
optim(0,fn,method="Brent",lower=0,upper=100)

#posterior mean
(data+a)*b/(b+1)

#95% symmetrical density
qgamma(c(0.025,0.975),57, scale=6/7)

#95% hpd interval
library(HDInterval)
tst <- rgamma(1e6, 57, scale=6/7)
hdi(tst, credMass=0.95)

#change the prior
a=8
b=5

#prior
plot(dgamma(x,shape=a,scale=b),type='l')

#likelihood
xx=seq(0,100,1)
data=52
plot(dpois(data,lambda=xx),type="l")

#posterior
plot(dgamma(x,shape=data+a,scale=b/(b+1)),type='l')

#library(optim.functions)
fn<-function(theta){
  ans=theta^59*exp((-6/5)*theta)
  return(-ans)
}
optim(0,fn,method="Brent",lower=0,upper=100)

#posterior mean
(data+a)*b/(b+1)

#95% symmetrical density
qgamma(c(0.025,0.975),60,scale=5/6)

#95% hpd interval
tst <- rgamma(1e7, 60, scale=5/6)
hdi(tst, credMass=0.95)

