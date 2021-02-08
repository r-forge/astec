// Libraries
#include <RcppGSL.h>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <limits>
#include <string>
#include <math.h>
#include <malloc.h>
#include <list>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <R.h>
#include <Rmath.h>

//g++ -O3 -o test_link test_link.cpp -lgsl -lgslcblas -lm 

using namespace std; 
extern "C" {
void metacommunity_model(double *parameters,unsigned int *stats_output){

  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();

  r = gsl_rng_alloc (gsl_rng_default);
  unsigned long int seed = 2018;
  gsl_rng_set (r, seed);
  
  unsigned int K=4;
  unsigned int n=10;
  unsigned int *x;
  x = new unsigned int[4];
  for (int i=0;i<4;i++){
	x[i]=0;
  }
  double *p;
  p = new double[4];
  for (int i=0;i<4;i++){
	p[i]=parameters[i];
  }
  gsl_ran_multinomial(r,K,n,p,x);
  for (int i=0;i<4;i++){
	  stats_output[i]=x[i];
  }
  
  gsl_rng_free (r);
  
}
}

