# below codes are based on r. please go to line 80.
library(mvtnorm)

#generate random number X and Y
identmat=diag(4)
mean=c(0,0,0,0)
n<-300000
#n<-1000
X=rmvnorm(n,mean,identmat)
Truebeta = matrix(c(0.5,-0.5,0,1),4,1)
Y=rbinom(n,1,exp(X%*%Truebeta)/(1+exp(X%*%Truebeta)))

#initial Beta
startB=matrix(c(0,0,0,0),4,1)


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
  lambda<-exp(XB)/(exp(XB)+1)
  sumlikelihood = sum(Y%*%log(lambda)+(1-Y)%*%log(1-lambda))
# sumlikelihood = sum(dbinom(Y,1,lambda,log=TRUE))
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
proposalsd=c(0.05,0.07,0.09,0.05)

proposalfunction <- function(B,p){
  return(rnorm(p,mean = B, sd= proposalsd))
}


#run MCMC 
iteration=10000
chain= runMCMC(startB,iteration,4)

acceptcount=chain$cnt
chain=chain$ch
plot(chain[,1],type='l')
plot(chain[,2],type='l')
plot(chain[,3],type='l')
plot(chain[,4],type='l')



####based on R code, I wrote rcpp function and use it
library(Rcpp)
setwd("C:\\Users\\SAMSUNG\\Desktop\\3학년\\2학기\\베이즈통계")
sourceCpp("hw5rcpp.cpp")

library(mvtnorm)
#generate random number X and Y
identmat=diag(4)
mean=c(0,0,0,0)
n<-300000

X=rmvnorm(n,mean,identmat)
Truebeta = matrix(c(0.5,-0.5,0,1),4,1)
Y=rbinom(n,1,exp(X%*%Truebeta)/(1+exp(X%*%Truebeta)))
Y=matrix(Y)

#log-likelihood
Rcpplikelihood(X,matrix(Y),Truebeta)
#compare with dbinom
lambda<-exp(X%*%Truebeta)/(exp(X%*%Truebeta)+1)
sumlikelihood = sum(dbinom(Y,1,lambda,log=TRUE));sumlikelihood

##log prior
priormean=c(0,0,0,0)
priorsd=c(sqrt(10),sqrt(10),sqrt(10),sqrt(10))         
logprior(Truebeta,priormean,priorsd)
#compare with dnorm
sum(dnorm(Truebeta,priormean,sd=priorsd,log=TRUE))

#proposal sd
proposalsd<-c(0.01,0.01,0.01,0.01)    #arbitrary proposal sd
iter=2000
#initial Beta
startB=matrix(c(0,0,0,0),4,1)
p=length(startB)
#sourceCpp("hw5rcpp.cpp")
data<-runMCMC(X,Y,startB,iter,priorsd,proposalsd,p)

count=0
for(i in 1:iter){
  if (data[i,1]!=data[i+1,1]){
    count=count+1
  }
}

#acceptancerate
print(count/iter)   


#traceplot
par(mfrow=c(2,2))
plot(data[,1],type='l')
plot(data[,2],type='l')
plot(data[,3],type='l')
plot(data[,4],type='l')

#burnin
newdata=data[400:2000,1:4]
#density plot
plot(density(newdata[,1]))
plot(density(newdata[,2]))
plot(density(newdata[,3]))
plot(density(newdata[,4]))

#hpd interval
library(HDInterval)
hdi(newdata[,1],0.95)
hdi(newdata[,2],0.95)
hdi(newdata[,3],0.95)
hdi(newdata[,4],0.95)

#posterior mean
mean(newdata[,1])
mean(newdata[,2])
mean(newdata[,3])
mean(newdata[,4])

#effective sample size
library(coda)
effectiveSize(newdata[,1])
effectiveSize(newdata[,2])
effectiveSize(newdata[,3])
effectiveSize(newdata[,4])


#########pararrel computing
library(parallel)
detectCores()
