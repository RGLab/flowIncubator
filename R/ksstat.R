#Copyright (c) 2017, Greg Finak
#' kstest
#' An Rcpp implementation returning the ks statistic between two numeric vectors.
#' @useDynLib flowIncubator
#' @param x vector a
#'
#' @param y vector b
#'
#' @importFrom Rcpp sourceCpp
#' @export
ksstat = function(x,y){
	D = kstest_c(x,y)
	while(D>1){
		D = kstest_c(x,y)
	}
	return(D)
}