library(Rcpp)
setwd("C:\\Users\\SAMSUNG\\Desktop\\3학년\\2학기\\베이즈통계")
sourceCpp("Rcpppr.cpp")
Niter<-10000
initmu<-0
initsigma<-5

MCMC_Rcpp<-RcppMCMC(Niter,initmu, initsigma)

hist(MCMC_Rcpp)
