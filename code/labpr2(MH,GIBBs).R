#construct the gibbs sampler.(1000 numbers of X, Y  are given)
library(invgamma)

#initial
b_0init = 0.1
b_1init = 0.2
sigsqinit = 0.001

#prior params (normal, normal, invgamma)
mu=0
n=1000
tau=100
rate=0.01
shape=0.01

b0=c()
b1=c()
sigsq=c()
for (i in 1:1000){
  #target b0 given others
  M=sum(Y-b_1init*X)/sigsqinit+mu/tau
  V=n/sigsqinit+1/tau
  nextb0=rnorm(1,M/V,sqrt(1/V))
  b_0init=nextb0
  b0<-c(b0,nextb0)
  
  #target b1 given others
  M=sum((Y-b_0init)*X)/sigsqinit+mu/tau
  V=sum(X^2)/sigsqinit+1/tau
  nextb1=rnorm(1,M/V,sqrt(1/V))
  b_1init=nextb1
  b1<-c(b1,nextb1)
  
  #target sigsq given others
  newrate=sum((Y-b_0init-b_1init*X)**2)/2+rate
  nextsigsq=rinvgamma(1,n/2+shape,newrate)
  sigsqinit=nextsigsq
  sigsq<-c(sigsq,nextsigsq)
}
  
#traceplot
plot(b0,type='l')
plot(b1,type='l')
plot(sigsq,type='l')

#compare it with lm function
lm(Y~X)
plot(X,Y)

#estimated formula in frequentist is Y = -0.9936X+0.9838
mean(b0)
mean(b1)
#we can see that the result is almost same!(freq vs bayes)

### now we have to deal with it by MH algorithm.###

prior = function(param){
  a = param[1]
  b = param[2]
  c = param[3]
  aprior = dnorm(a, sd = sqrt(tau), log = T)
  bprior = dnorm(b, sd = sqrt(tau), log = T)
  cprior = dinvgamma(c, 0.01, 0.01, log = T)
  return(aprior+bprior+cprior)
}

likelihood <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  
  pred = a*X + b
  singlelikelihoods = dnorm(Y, mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)   
}

proposalfunction <- function(param){
  return(rnorm(3,mean = param, sd= c(0.05,0.1,0.005)))
}

posterior <- function(param){
  return (likelihood(param) + prior(param))
}

startvalue = c(0,0.5,0.8)
acceptcount=0
saveb0=c()
saveb1=c()
savesigsq=c()
for (i in 1:10000){
  #generate parameter from proposal
  propparam=proposalfunction(startvalue)
  
  #Accept rate via MH algorithm
  if (log(runif(1,0,1))<posterior(propparam)-posterior(startvalue)){
    startvalue<-propparam
    saveb0<-c(saveb0,startvalue[2])
    saveb1<-c(saveb1,startvalue[1])
    savesigsq<-c(savesigsq,startvalue[3])
    acceptcount=acceptcount+1
  }
}
mean(saveb0)
mean(saveb1)
mean(savesigsq)
acceptcount
plot(saveb0,type='l')
plot(saveb1,type='l')
plot(savesigsq,type='l')
