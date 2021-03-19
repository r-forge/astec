#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector make_average_logbodymass(NumericVector a) {

 // Declare loop counters, and vector sizes
 int i, na = a.size();

 // Create vector filled with 0
 NumericVector av((na-1));

 // Crux of the algorithm
 for(i =0; i < (na-1); i++) {
   av[i] = exp(a[i]) * (exp(a[(na-1)])-1)/a[(na-1)];
 }

 // Return result
 return av;
}

