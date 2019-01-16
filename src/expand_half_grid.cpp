#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List expand_half_grid(CharacterVector x){
  int n = x.size();
  int m = n * (n-1) / 2;
  CharacterVector l(m), r(m);
  int k=0;
  for(int i=0; i<n-1; i++) {
    for(int j=i+1; j<n; j++) {
      l[k]=x[i];
      r[k]=x[j];
      k++;
    }
  }
  return List::create(l, r);
}

// [[Rcpp::export]]
List expand_r(NumericMatrix x) {
  int n=x.nrow();
  int m = n * (n-1) / 2;
  NumericVector I(m), J(m), R(m);
  int k=0; // index vector
  for(int i=0; i<n-1; i++) {
    for(int j=i+1; j<n; j++) {
      I(k) = i+1;
      J(k) = j+1;
      R(k) = x(i,j);
      k++;
    }
  }
  return List::create(I, J, R);
}
  
