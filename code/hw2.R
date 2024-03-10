#HW2 Question1
#Cauchy dist as a target distribution
count=0
#proposal dist : Normal Distribution
init=0
std=15
fn<-function(x,theta=2,eta=0){
  ans=(1+((x-eta)/theta)^2)*theta*pi
  ans=1/ans
  return(ans)
}
save_x<-c()

#MH algorithm
for(i in 1:10000){
  sampx=rnorm(1,mean=init,std)
  sampunif<-runif(1,0,1)
  if (sampunif<fn(sampx)/fn(init)) {
    init<-sampx
    count<-count+1}
  
  else{
    init=init}
  save_x<-c(save_x,init)
  
}
save_x

acceptrate<-count/10000;acceptrate            #acceptance probability
hist(save_x)               #sampled target distribution pdf(histogram)
plot(save_x,type='l')      #traceplot

plot(density(save_x),col='red')
par(new=TRUE)    #or, using points(~~~)
plot(dcauchy(seq(min(save_x),max(save_x),1),0,2),type='l')

#25th percentile
a<-qcauchy(0.25,location=0, scale=2);a
#50th percentile
b<-qcauchy(0.5,location=0, scale=2);b
#75th percentile
c<-qcauchy(0.75,location=0, scale=2);c

#percentile of my sample
sort(save_x)[2500]
sort(save_x)[5000]
sort(save_x)[7500]

x  #I load 1000 samples generated from cauchy distribution stored as x
library(HDInterval)


#HW2-Question2
acceptcounteta=0
acceptcounttheta=0
savetheta<-c()
saveeta<-c()

#arbitrary initial values
initeta=2.5
inittheta=2

#proposed normal dist std
initstdeta=0.1   
initstdtheta=0.2

#sum of log-likelihood function
likelihood = function(param){
  eta = param[1]
  theta = param[2]
  singlelikelihoods = dcauchy(x, location=eta, scale=theta, log = T)
  sumlikelihood = sum(singlelikelihoods)
  return(sumlikelihood)
}

# Prior distribution(log)
prior = function(param){
  a = param[1]
  b = param[2]
  aprior = dnorm(a, sd = 10, log = T)
  bprior = dunif(b, min=0, max=10, log = T)
  
  return(aprior+bprior)
  
}
#posterior distribution(log transformed)
posterior = function(param){
  return (likelihood(param) + prior(param))
}

#my proposal is normal distribution
proposalfunction = function(param){
  return(rnorm(2,mean = param, sd= c(initstdeta,initstdtheta)))
}

for (i in 1:10000){
  #generate parameter from proposal
  propeta=proposalfunction(c(initeta,inittheta))[1]
  
  #Accept rate via MH algorithm
  if (log(runif(1,0,1))<posterior(c(propeta,inittheta))-posterior(c(initeta,inittheta))){
    initeta<-propeta
    saveeta<-c(saveeta,initeta)
    acceptcounteta=acceptcounteta+1
  }
  proptheta=proposalfunction(c(initeta,inittheta))[2]
  #as theta>0, we have to deal with it
  if (proptheta<0){next}
  
  # Accept rate
  if (log(runif(1,0,1))<posterior(c(initeta,proptheta))-posterior(c(initeta,inittheta))){
    inittheta<-proptheta
    savetheta<-c(savetheta,inittheta)
    acceptcounttheta=acceptcounttheta+1
  }
}
acceptcounteta
acceptcounttheta
print(paste("accept rate of eta and theta :", acceptcounteta/10000, acceptcounttheta/10000))
hist(saveeta)
hist(savetheta)
plot(saveeta,type='l')
plot(savetheta,type='l')
mean(saveeta)
mean(savetheta)
hdi(saveeta,0.95)
hdi(savetheta,0.95)
