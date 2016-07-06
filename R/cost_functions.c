#include <R.h>
#include <Rmath.h>
#include <Rinternals.h> // RK addition
#include <R_ext/RS.h>  // RK addition
#include <R_ext/Lapack.h> // RK addition
#include <R_ext/BLAS.h> // RK addition
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
//#include "cost_general_functions.c"
  #define SWAP(a,b)   { int t; t=a; a=b; b=t; }  // Macro for swapping

double mll_nonparametric_ed(double x, double x2, double x3, int n, double shape){
  return(x2-(x*x)/n);
}

double mbic_mean(double x, double x2, double x3, int n, double shape){
  return(x2-(x*x)/n+log(n));
}
