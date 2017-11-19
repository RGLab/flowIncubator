//
// Copyright (c) 2017 by Greg Finak. All Rights Reserved.
//
#include <Rcpp.h>
#include <stdlib.h>
using namespace Rcpp;

// Rcpp function to compute the ks statistic
//
// [[Rcpp:::plugins(cpp14)]]
// [[Rcpp::export]]
double kstest_c(NumericVector x, NumericVector y) {
	x = na_omit(x);
	y = na_omit(y);
	int nx = x.size();
	int ny = y.size();
	Environment base("package:base");
	Function order = base["order"];
	NumericVector w = NumericVector(x.size()+y.size());
	std::copy(x.begin(),x.end(),w.begin());
	std::copy(y.begin(),y.end(),w.begin()+x.size());
	NumericVector srt = clone(w);
	std::sort(srt.begin(),srt.end(),[](const double& lhs, const double& rhs){
		return std::tie(lhs,lhs) < std::tie(rhs,rhs);
	});
	w = order(w,Named("method","radix"));
	NumericVector z = NumericVector(w.size());
	std::transform(w.begin(),w.end(),z.begin(),[&](double d) {
		if(d<=nx){
			return(1.0/(double)nx);
		}else{
			return(-1.0/(double)ny);
		}
	});
	NumericVector cz = NumericVector(z.size());
	std::partial_sum(z.begin(),z.end(),cz.begin(),[](double a, double b){
		return(a+b);
	});
	NumericVector df = diff(srt);
	std::vector<int> ind(df.size());
	int j=0;
	for(int i=0;i<df.size();i++){
		if(df[i]!=0.0) ind[j++] = i+1;
	}
	ind.resize(j);
	IntegerVector iind = wrap(ind);
	NumericVector I = cz[iind];
	I.push_back(cz[nx+ny]);
	double D = max(abs(I));
	return(D);
}



/*** R
library(microbenchmark)
	set.seed(100)
	x <- c(rnorm(5000),rep(NA_integer_,10),rep(5,1000))
	y <- c(rep(NA_integer_,10),runif(3000),rep(5,1000))
	nx = length(na.omit(x))
	ny = length(na.omit(y))
# ks.test(x,y)$statistic
	
	
	all.equal(kstest_c(x,y),ks.test(x,y)$statistic[[1]])
	microbenchmark(kstest_c(x,y),ks.test(x,y)$statistic[[1]])
*/
