library(mvtnorm)
library(invgamma)
#generate random number X and Y
identmat=diag(10)
mean=rep(0,10)
n<-10000

X=rmvnorm(n,mean,identmat)
Truebeta = matrix(c(0.5,-0.5,1,-1,0.7,0,0,0,0,0),10,1)
Y=rbinom(n,1,exp(X%*%Truebeta)/(1+exp(X%*%Truebeta)))


#sum of log-likelihood function
likelihood = function(B){
  XB<-X%*%B
  lambda<-exp(XB)/(exp(XB)+1)
  sumlikelihood = sum(dbinom(Y,1,lambda,log=TRUE))
  return(sumlikelihood)
}

# log beta prior 
priorb = function(B,indicator,sigmasq1,sigmasq2){
  res=sum(indicator*dnorm(B,mean=0, sd=sqrt(sigmasq1), log=TRUE)+(1-indicator)*dnorm(B,mean=0,sd=sqrt(sigmasq2), log=TRUE))
  return (res)
}

#logsum of indicator prior
priorind = function(indicator){
  res=sum(dbinom(indicator,1,0.5,log=TRUE))
  return(res)
}

#log sigmasq1 prior
priorsigsq1 = function(sigmasq1){
  res=dinvgamma(sigmasq1,1,20,log=TRUE)
  return(res)
}

#log simgasq2 prior
priorsigsq2 =function(sigmasq2){
  res=dgamma(sigmasq2,1,20, log=TRUE)
  return(res)
}

posterior = function(B,indicator,sigmasq1,sigmasq2){
  res = likelihood(B)+priorb(B,indicator,sigmasq1,sigmasq2)+priorind(indicator)+priorsigsq1(sigmasq1)+priorsigsq2(sigmasq2)
  return(res)
}

#initial values
initsigmasq1<-rinvgamma(1,1,20)
initsigmasq2<-rgamma(1,1,20)
initindicator<-c(1,0,1,0,1,0,1,0,1,0)
initbetas<-rep(0,10)
iter<-10000
p<-length(initbetas)

beta.samps <- ind.samps <-matrix(NA, nrow = p, ncol = iter)
beta.samps[,1] <- initbetas
ind.samps[,1] <- initindicator

s1.samps <- s2.samps <- rep(NA, iter)
s1.samps[1] <- initsigmasq1
s2.samps[1] <- initsigmasq2

#proposal sd
propsdbeta<-0.015
propsdind<-0.5
propsdsigsq1<-20
propsdsigsq2<-0.02

#checking acceptance rate
count1<-count2<-count3<-0


for(i in 2:iter){
  #propose betas
  newbeta<-c()
  for(j in 1:10){
    newbeta=c(newbeta,rnorm(1, mean=beta.samps[j,i-1], sd=propsdbeta))
  }
  if(posterior(newbeta,ind.samps[,i-1],s1.samps[i-1],s2.samps[i-1])-posterior(beta.samps[,i-1],ind.samps[,i-1],s1.samps[i-1],s2.samps[i-1])>log(runif(1))){
       beta.samps[,i]<-newbeta
       count1=count1+1
  }else{beta.samps[,i]<-beta.samps[,i-1]}
  
  #propose indicator
  newindicator<-c()
  for(k in 1:10){
    a<-rnorm(1,mean=ind.samps[k,i-1], sd = propsdind)
    if(a>=0.5){
      newindicator=c(newindicator,1)
    }else{newindicator=c(newindicator,0)}
  }
  if(posterior((beta.samps[,i]),newindicator,s1.samps[i-1],s2.samps[i-1])-posterior((beta.samps[,i]),ind.samps[,i-1],s1.samps[i-1],s2.samps[i-1])>log(runif(1))){
    ind.samps[,i]<-newindicator
  }else{ind.samps[,i]<-ind.samps[,i-1]}
  
  #propose sigmasq1
  newsigmasq1<-rnorm(1,mean=s1.samps[i-1],sd=propsdsigsq1)
  if(newsigmasq1<0){
    s1.samps[i]<-s1.samps[i-1]
  }else{
  if(posterior(beta.samps[,i],ind.samps[,i],newsigmasq1,s2.samps[i-1])-posterior(beta.samps[,i],ind.samps[,i],s1.samps[i-1],s2.samps[i-1])>log(runif(1))){
    s1.samps[i]<-newsigmasq1
    count2=count2+1
  }else{s1.samps[i]<-s1.samps[i-1]}}
  
  #propose sigmasq2
  newsigmasq2<-rnorm(1,mean=s2.samps[i-1],sd=propsdsigsq2)
  if(newsigmasq2<0){
    s2.samps[i]<-s2.samps[i-1]
  }else{
  if(posterior(beta.samps[,i],ind.samps[,i],s1.samps[i],newsigmasq2)-posterior(beta.samps[,i],ind.samps[,i],s1.samps[i],s2.samps[i-1])>log(runif(1))){
    s2.samps[i]<-newsigmasq2
    count3=count3+1
  }else{s2.samps[i]<-s2.samps[i-1]}}
}

#traceplot
par(mfrow=c(3,4))
plot(beta.samps[1,],type='l')
plot(beta.samps[2,],type='l')
plot(beta.samps[3,],type='l')
plot(beta.samps[4,],type='l')
plot(beta.samps[5,],type='l')
plot(beta.samps[6,],type='l')
plot(beta.samps[7,],type='l')
plot(beta.samps[8,],type='l')
plot(beta.samps[9,],type='l')
plot(beta.samps[10,],type='l')
plot(s1.samps,type='l')
plot(s2.samps,type='l')

#density plot
plot(density(beta.samps[1,]))
plot(density(beta.samps[2,]))
plot(density(beta.samps[3,]))
plot(density(beta.samps[4,]))
plot(density(beta.samps[5,]))
plot(density(beta.samps[6,]))
plot(density(beta.samps[7,]))
plot(density(beta.samps[8,]))
plot(density(beta.samps[9,]))
plot(density(beta.samps[10,]))
plot(density(s1.samps))
plot(density(s2.samps))

#burnin
burninbeta=beta.samps[,1000:iter]
burninsigsq1=s1.samps[1000:iter]
burninsigsq2=s2.samps[1000:iter]

#now I calculate hpd interval and posterior mean based on burnin data
#hpd interval
library(HDInterval)
hdi(burninbeta[1,],0.95)
hdi(burninbeta[2,],0.95)
hdi(burninbeta[3,],0.95)
hdi(burninbeta[4,],0.95)
hdi(burninbeta[5,],0.95)
hdi(burninbeta[6,],0.95)
hdi(burninbeta[7,],0.95)
hdi(burninbeta[8,],0.95)
hdi(burninbeta[9,],0.95)
hdi(burninbeta[10,],0.95)
hdi(burninsigsq1,0.95)
hdi(burninsigsq2,0.95)

#posterior mean
mean(burninbeta[1,])
mean(burninbeta[2,])
mean(burninbeta[3,])
mean(burninbeta[4,])
mean(burninbeta[5,])
mean(burninbeta[6,])
mean(burninbeta[7,])
mean(burninbeta[8,])
mean(burninbeta[9,])
mean(burninbeta[10,])
mean(burninsigsq1)
mean(burninsigsq2)

#acceptance rate of each parameter
count1/iter         #beta
count2/iter         #sigmasq1
count3/iter         #sigmasq2


#effective sample size
library(coda)
effectiveSize(burninbeta[1,])
effectiveSize(burninbeta[2,])
effectiveSize(burninbeta[3,])
effectiveSize(burninbeta[4,])
effectiveSize(burninbeta[5,])
effectiveSize(burninbeta[6,])
effectiveSize(burninbeta[7,])
effectiveSize(burninbeta[8,])
effectiveSize(burninbeta[9,])
effectiveSize(burninbeta[10,])
effectiveSize(burninsigsq1)
effectiveSize(burninsigsq2)

#posterior inclusion probability
for(i in 1:p){
  print(mean(ind.samps[i,]))
}




###use nimble package to solve it
library(nimble)
Xmat=X


model_glm<-nimbleCode({
  
  #Data Model
  for (i in 1:n){
    XB[i]<-inprod(X[i,1:10],beta[1:10])
    lambda[i] <-exp(XB[i])/(exp(XB[i])+1)
    Y[i] ~ dbinom(prob=lambda[i], size=1)
  }
  #parameter Model
  slabprior ~ dinvgamma(1,20)
  spikeprior ~ dgamma(1,20)
  
  for (i in 1:p){
    indprior[i]~ dbern(0.5)
    beta[i]~dnorm(0,sd=sqrt(slabprior)*indprior[i]+sqrt(spikeprior)*(1-indprior[i]))
    }
  
})

niter=10000
consts<-list(n=dim(Xmat)[1], p=dim(Xmat)[2],X=Xmat)
dat<-list(Y=Y)
inits<-list(slabprior=100,spikeprior=0.4, beta=rep(0,dim(Xmat)[2]),indprior=c(1,0,1,0,1,0,1,0,1,0))

#Run MCMC
pt<-proc.time()
samples_glm<-nimbleMCMC(model_glm,data=dat,inits=inits,
                        constants = consts,
                        monitors=c("beta","indprior",'slabprior','spikeprior'),
                        samplesAsCodaMCMC=TRUE,WAIC=FALSE,summary=FALSE,
                        niter=niter,nburnin=0,thin=1,nchains=1)
ptFinal_glm<-proc.time()-pt
ptFinal_glm

#head(samples_glm)

#traceplot of betas, slabprior(sigmasq1), spikeprior(sigmasq2)
par(mfrow=c(3,4))
plot(samples_glm[1:10000,1],type='l')
plot(samples_glm[1:10000,2],type='l')
plot(samples_glm[1:10000,3],type='l')
plot(samples_glm[1:10000,4],type='l')
plot(samples_glm[1:10000,5],type='l')
plot(samples_glm[1:10000,6],type='l')
plot(samples_glm[1:10000,7],type='l')
plot(samples_glm[1:10000,8],type='l')
plot(samples_glm[1:10000,9],type='l')
plot(samples_glm[1:10000,10],type='l')
plot(samples_glm[1:10000,21],type='l')
plot(samples_glm[1:10000,22],type='l')

#density plot of betas, slabprior, spikeprior
plot(density(samples_glm[1:10000,1]))
plot(density(samples_glm[1:10000,2]))
plot(density(samples_glm[1:10000,3]))
plot(density(samples_glm[1:10000,4]))
plot(density(samples_glm[1:10000,5]))
plot(density(samples_glm[1:10000,6]))
plot(density(samples_glm[1:10000,7]))
plot(density(samples_glm[1:10000,8]))
plot(density(samples_glm[1:10000,9]))
plot(density(samples_glm[1:10000,10]))
plot(density(samples_glm[1:10000,21]))
plot(density(samples_glm[1:10000,22]))

#in nimble case, betas are converged very well so I do not proceed burnin
#hpd interval
library(HDInterval)
hdi(samples_glm[1:10000,1],0.95)
hdi(samples_glm[1:10000,2],0.95)
hdi(samples_glm[1:10000,3],0.95)
hdi(samples_glm[1:10000,4],0.95)
hdi(samples_glm[1:10000,5],0.95)
hdi(samples_glm[1:10000,6],0.95)
hdi(samples_glm[1:10000,7],0.95)
hdi(samples_glm[1:10000,8],0.95)
hdi(samples_glm[1:10000,9],0.95)
hdi(samples_glm[1:10000,10],0.95)
hdi(samples_glm[1:10000,21],0.95)
hdi(samples_glm[1:10000,22],0.95)

#posterior mean
mean(samples_glm[1:10000,1])
mean(samples_glm[1:10000,2])
mean(samples_glm[1:10000,3])
mean(samples_glm[1:10000,4])
mean(samples_glm[1:10000,5])
mean(samples_glm[1:10000,6])
mean(samples_glm[1:10000,7])
mean(samples_glm[1:10000,8])
mean(samples_glm[1:10000,9])
mean(samples_glm[1:10000,10])
mean(samples_glm[1:10000,21])
mean(samples_glm[1:10000,22])

#acceptance prpobability
#acceptance probability of betas
pp<-0
for (i in 1:9999){
  if (samples_glm[,1][i]==samples_glm[,1][i+1]){
    pp=pp+1
  }
}
acceptratebeta<-(niter-pp)/niter

#acceptance probability of sigmasq1
ppp<-0
for (i in 1:9999){
  if (samples_glm[,21][i]==samples_glm[,21][i+1]){
    ppp=ppp+1
  }
}
acceptratesigmasq1<-(niter-ppp)/niter

#acceptance probability of sigmasq2
pppp<-0
for (i in 1:9999){
  if (samples_glm[,22][i]==samples_glm[,22][i+1]){
    pppp=pppp+1
  }
}
acceptratesigmasq2<-(niter-pppp)/niter
print(c(acceptratebeta,acceptratesigmasq1,acceptratesigmasq2))

#effective sample size
library(coda)
effectiveSize(samples_glm[,1])
effectiveSize(samples_glm[,2])
effectiveSize(samples_glm[,3])
effectiveSize(samples_glm[,4])
effectiveSize(samples_glm[,5])
effectiveSize(samples_glm[,6])
effectiveSize(samples_glm[,7])
effectiveSize(samples_glm[,8])
effectiveSize(samples_glm[,9])
effectiveSize(samples_glm[,10])
effectiveSize(samples_glm[,21])
effectiveSize(samples_glm[,22])

#posterior inclusion probability
for(i in 11:20){
  print(mean(samples_glm[1:10000,i]))
}

