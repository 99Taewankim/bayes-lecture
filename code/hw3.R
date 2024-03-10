#1. Do mcmc by own code
library(mvtnorm)
identmat=diag(4)
mean=c(0,0,0,0)
X=rmvnorm(1000,mean,identmat)
Truebeta = matrix(c(0.5,-0.5,0,1),4,1)
Y=rpois(1000,exp(X%*%Truebeta))



runMCMC <- function(startvalue, iterations,p){
  acceptcount=0
  chain = array(dim = c(iterations+1,p))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,],p)
    
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
      acceptcount=acceptcount+1
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  print(acceptcount)
  return(list(ch=chain,cnt=acceptcount))
}


#sum of log-likelihood function
likelihood = function(B){
  XB<-X%*%B
  lambda<-exp(XB)
  singlelikelihoods = sum(log(exp(-lambda)*lambda**Y/factorial(Y)))
  sumlikelihood = sum(singlelikelihoods)
  return(sumlikelihood)
}

# Prior distribution(log)
prior = function(B){
  res=sum(dnorm(B,mean=0, sd=sqrt(10), log=TRUE))
  return (res)
}
  

#posterior distribution(log target)
posterior = function(B){
  return (likelihood(B) + prior(B))
}

#propose beta from normal with arbitrary sd
proposalsd=c(0.01,0.03,0.02,0.05)

proposalfunction <- function(B,p){
  return(rnorm(p,mean = B, sd= proposalsd))
}

#initial Beta
startB=c(0,0,0,0)

#run MCMC 
iteration=10000
chain= runMCMC(startB,iteration,4)

acceptcount=chain$cnt
chain=chain$ch

#acceptance probability
acceptcount/iteration

#effective sample size
library(coda)
effectiveSize(chain)

#####overcome autocorr problem!
#First way is thinning
iterations <- 20 * 10000 
chain <- runMCMC(startB, iterations,4)
totcount<-chain$cnt
chain=chain$ch
# Select values every 4 iterations for the final chain
cfinal <- matrix(NA, ncol = 4, nrow = 50000)
for (i in 1:50000){
  cfinal[i, ] <- chain[i*4,]
}

#effective sample size
effectiveSize(cfinal)

#acceptance probability
totcount/iterations

#second way is deleting front items(I don't use this method in report)
# Remove the first 5000 values of the chain
#burnIn <- 5000
#effectiveSize(chain[-(1:burnIn),1])
#effectiveSize(chain[-(1:burnIn),2])
#effectiveSize(chain[-(1:burnIn),3])
#effectiveSize(chain[-(1:burnIn),4])

#traceplot
#par(mfrow=c(2,2))
plot(cfinal[,1],type='l')
plot(cfinal[,2],type='l')
plot(cfinal[,3],type='l')
plot(cfinal[,4],type='l')

#histogram
hist(cfinal[,1])
hist(cfinal[,2])
hist(cfinal[,3])
hist(cfinal[,4])

#HPD Interval
library(HDInterval)
hdi(cfinal[,1],0.95)
hdi(cfinal[,2],0.95)
hdi(cfinal[,3],0.95)
hdi(cfinal[,4],0.95)

#posterior mean(BjÀÇ »çÈÄ Æò±Õ)
mean(cfinal[,1])
mean(cfinal[,2])
mean(cfinal[,3])
mean(cfinal[,4])

###2. using nimble
library(nimble)

identmat=diag(4)
mean=c(0,0,0,0)
Xmat=rmvnorm(1000,mean,identmat)
Truebeta = matrix(c(0.5,-0.5,0,1),4,1)
Y=rpois(1000,exp(Xmat%*%Truebeta))

model_glm<-nimbleCode({
  
  #Data Model
  for (i in 1:n){
    lambda[i] <-exp(XB[i])
    Y[i] ~ dpois(lambda=lambda[i])
  }
  XB[1:n]<-X[1:n,1:p]%*%beta[1:p]
  
  #parameter Model
  beta[1:p]~dmnorm(mean=M[1:p],cov=Cov[1:p,1:p])
})

niter=50000
consts<-list(n=dim(Xmat)[1], p=dim(Xmat)[2],X=Xmat , M=rep(0,dim(Xmat)[2]),Cov=(diag(dim(Xmat)[2]))*sqrt(10))
dat<-list(Y=Y)
inits<-list(beta=rep(0,dim(Xmat)[2]))

#Run MCMC
pt<-proc.time()
samples_glm<-nimbleMCMC(model_glm,data=dat,inits=inits,
                        constants = consts,
                        monitors=c("beta"),
                        samplesAsCodaMCMC=TRUE,WAIC=FALSE,summary=FALSE,
                        niter=niter,nburnin=0,thin=1,nchains=1)
ptFinal_glm<-proc.time()-pt
ptFinal_glm

#traceplot

plot(samples_glm[,1],type='l')
plot(samples_glm[,2],type='l')
plot(samples_glm[,3],type='l')
plot(samples_glm[,4],type='l')

#histogram
hist(samples_glm[,1])
hist(samples_glm[,2])
hist(samples_glm[,3])
hist(samples_glm[,4])

#HPD Interval
library(HDInterval)
hdi(samples_glm[,1],0.95)
hdi(samples_glm[,2],0.95)
hdi(samples_glm[,3],0.95)
hdi(samples_glm[,4],0.95)

#posterior mean(BjÀÇ »çÈÄ Æò±Õ)
mean(samples_glm[,1])
mean(samples_glm[,2])
mean(samples_glm[,3])
mean(samples_glm[,4])

#effective sample size
library(coda)
effectiveSize(samples_glm)

#accept probability
pp<-0
for (i in 1:49999){
  if (samples_glm[,1][i]==samples_glm[,1][i+1]){
      pp=pp+1
  }
}
acceptratenimble<-(niter-pp)/niter
acceptratenimble

#3. adaptMCMC
library(adaptMCMC)

identmat=diag(4)
mean=c(0,0,0,0)
Xmat=rmvnorm(1000,mean,identmat)
Truebeta = matrix(c(0.5,-0.5,0,1),4,1)
Y=rpois(1000,exp(Xmat%*%Truebeta))
proposalsd=c(0.01,0.03,0.02,0.05)
#initial Beta
init.pars=c(beta1=0,beta2=0,beta3=0,beta4=0)

#log posterior
logPost = function(pars){
  with(as.list(pars), {
    beta=pars
    lambda<-exp(Xmat%*%beta)
    #likelihood+prior
    sum(dpois(Y,lambda,log=T))+sum(dnorm(beta,mean=0, sd=sqrt(10),log=T))
  })
}

out.mcmc<-MCMC(p=logPost, n=50000, init=init.pars, scale=proposalsd, adapt=TRUE, acc.rate=.3)

#acceptance rate
out.mcmc$acceptance.rate

plot(out.mcmc$samples[,1],type='l')
plot(out.mcmc$samples[,2],type='l')
plot(out.mcmc$samples[,3],type='l')
plot(out.mcmc$samples[,4],type='l')

#histogram
hist(out.mcmc$samples[,1])
hist(out.mcmc$samples[,2])
hist(out.mcmc$samples[,3])
hist(out.mcmc$samples[,4])

#HPD Interval
library(HDInterval)
hdi(out.mcmc$samples[,1],0.95)
hdi(out.mcmc$samples[,2],0.95)
hdi(out.mcmc$samples[,3],0.95)
hdi(out.mcmc$samples[,4],0.95)

#posterior mean(BjÀÇ »çÈÄ Æò±Õ)
mean(out.mcmc$samples[,1])
mean(out.mcmc$samples[,2])
mean(out.mcmc$samples[,3])
mean(out.mcmc$samples[,4])

#effective sample size
library(coda)
effectiveSize(out.mcmc$samples)

#overplot all samples
#beta1
#par(mfrow=c(1,1))
plot(hist(cfinal[,1]),col='red',density=60,xlim=c(0,1),main="MCMC_samples_beta1",xlab='beta1')
hist(samples_glm[,1],col = 'green', density =40, add=T)
hist(out.mcmc$samples[,1],col='blue',add=T, density= 40)
#true beta
abline(v=0.5, col='black',lwd=2)

#beta2
plot(hist(cfinal[,2]),col='red',density=60,xlim=c(-1,0),main="MCMC_samples_beta2",xlab='beta2')
hist(samples_glm[,2],col = 'green', density =40, add=T)
hist(out.mcmc$samples[,2],col='blue',add=T, density= 40)
#true beta
abline(v=-0.5, col='black',lwd=2)

#beta3
plot(hist(cfinal[,3]),col='red',density=60,xlim=c(-1,1),main="MCMC_samples_beta3",xlab='beta3')
hist(samples_glm[,3],col = 'green', density =40, add=T)
hist(out.mcmc$samples[,3],col='blue',add=T, density= 40)
#true beta
abline(v=0, col='black',lwd=2)

#beta4
plot(hist(cfinal[,4]),col='red',density=60,xlim=c(0,2),main="MCMC_samples_beta4",xlab='beta4')
hist(samples_glm[,4],col = 'green', density =40, add=T)
hist(out.mcmc$samples[,4],col='blue',add=T, density= 40)
#true beta
abline(v=1, col='black',lwd=2)
