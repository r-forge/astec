#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector make_average_logbodymass_within(NumericVector a) {

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

// [[Rcpp::export]]
NumericMatrix make_average_logbodymass_matrix(NumericMatrix a) {

 // Declare loop counters, and vector sizes
 int i,j, na = a.nrow(), n=a.ncol();

 // Create vector filled with 0
 NumericMatrix av((na-1),n);

 // Crux of the algorithm
 for (j=0;j<n;j++){
  for(i =0; i < (na-1); i++) {
   av(i,j) = exp(a(i,j)) * (exp(a((na-1),j))-1)/a((na-1),j);
  }
 }

 // Return result
 return av;
}

// [[Rcpp::export]]
NumericMatrix make_average_logbodymass_matrix3(NumericMatrix a) {

 // Declare loop counters, and vector sizes
 int i,j, na = a.nrow(), n=a.ncol();

 // Create vector filled with 0
 NumericMatrix av((na-1),n);

 // Crux of the algorithm
 for (j=0;j<n;j++){
   av(_,j)=make_average_logbodymass_within(a(_,j));
 }

 // Return result
 return av;
}


// [[Rcpp::export]]
void make_average_logbodymass_matrix_void(NumericMatrix a ,NumericMatrix& av) {

 // Declare loop counters, and vector sizes
 int j, n=a.ncol();

 // Crux of the algorithm
 for (j=0;j<n;j++){
   av(_,j)=make_average_logbodymass_within(a(_,j));
 }
}

// [[Rcpp::export]]
void make_average_logbodymass_matrix_void2(NumericMatrix a ,NumericMatrix& av) {

 // Declare loop counters, and vector sizes
 int i,j, na = a.nrow(), n=a.ncol();

 // Crux of the algorithm
 for (j=0;j<n;j++){
  for(i =0; i < (na-1); i++) {
   av(i,j) = exp(a(i,j)) * (exp(a((na-1),j))-1)/a((na-1),j);
  }
 }
}



// [[Rcpp::export]]
void multiply_void(NumericMatrix& a ,double k) {

 // Declare loop counters, and vector sizes
 int i,j, n=a.ncol(), m=a.nrow();

 // Crux of the algorithm
 for (i=0;i<m;i++){
  for (j=0;j<n;j++){
   a(i,j)=k*a(i,j);
  }
 }
}


