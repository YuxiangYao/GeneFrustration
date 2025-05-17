// C++ for local frustration

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::IntegerVector LocalInitialFrustration(
    IntegerVector vvv, int Max, int Min){
  int ii,jj;
  int n_sam=vvv.size();
  IntegerVector res(Max-Min+1,0);
  for(ii=0; ii<n_sam; ++ii){
    jj=vvv[ii]-Min;
    res[jj]=res[jj]+1;
  }
  return res;
}


