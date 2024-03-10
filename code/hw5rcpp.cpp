#include <RcppArmadillo.h>
#include <omp.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]


NumericMatrix mmult(NumericMatrix m , NumericMatrix v , bool byrow=true ){
	NumericMatrix res(m.nrow(),v.ncol());
  	if( m.ncol() != v.nrow() ) stop("Non-conformable arrays") ;

  	for (int i = 0; i < m.nrow(); i++){
  		for (int j = 0; j < v.ncol(); j++){
  			for (int k = 0; k<m.ncol(); k++){
  				res(i,j)+=m(i,k)*v(k,j);
			  }
		  }
	  }

  return res ;
  }

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

NumericMatrix fff(vec y, int n){
	NumericMatrix m(n,1);
	for (int i = 0; i < n; i++){
		m(i,0)=y(i);
	}
	return m;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double Rcpplikelihood(NumericMatrix X, NumericMatrix Y, NumericMatrix B){
	vec lambda;
	lambda=exp(mmult(X,B))/(exp(mmult(X,B))+ 1);
	NumericMatrix m=fff(log(lambda),300000);
	NumericMatrix mm=fff(log(1-lambda),300000);
	vec likelihood = mmult(transpose(Y),m)+mmult((1-transpose(Y)),mm);
return(sum(likelihood));
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double logprior(vec x, vec means, vec priorsd) {
  int n = x.size();
  vec res(n);
  for(int i = 0; i < n; i++) {
    res(i) =log(1/(priorsd(i)*pow(2*M_PI,0.5))*exp(-pow(x(i)-means(i),2)/(2*priorsd(i)*priorsd(i))));
  }
  
  return sum(res);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

	
vec proposalf(vec B, vec sd, int p){
	vec res(p);
	for(int j = 0; j< p; j++){
		res(j)=sd(j)*randn()+B(j);
	}
	return res;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
mat runMCMC(NumericMatrix X, NumericMatrix Y, NumericMatrix initB, int iterations, vec priorsd, vec proposalsd, int p){
	vec priormean = zeros(p);
	mat chain = zeros(iterations+1,p);
	//putting inital Beta 
	for(int k=0; k<p; k++){
		chain(0,k)=initB(k,0);
	}
	for(int i =0; i<iterations;i++){
		mat B =chain.row(i);
		vec b = zeros(p);
		for(int kk=0; kk<p; kk++){
			b(kk)=B(0,kk);
		}
		vec prop = proposalf(b,proposalsd,p);
		NumericMatrix matB=fff(b,p);
		NumericMatrix matprop=fff(prop,p);
		double prob = logprior(prop,priormean,priorsd)+Rcpplikelihood(X,Y,matprop)-logprior(b,priormean,priorsd)-Rcpplikelihood(X,Y,matB);
		if(  log(randu() ) < prob ){
			for (int idx=0; idx<p;idx++){
				chain(i+1,idx)=prop(idx);
			}
	    }else{
	    	for (int idx=0; idx<p;idx++){
	    		chain(i+1,idx)=chain(i,idx);
			}
        }
	}
	return chain;
	
}
