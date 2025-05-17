// Smooth the PCA-plots scores (landscape).
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector SmoothKNN(IntegerMatrix knn_info,IntegerVector Socres,int K){
  int ii,jj,nsam=knn_info.nrow();
  NumericVector Res(nsam);
  double sum=0;
  for(ii=0;ii<nsam;ii++){
    sum=0;
    for(jj=0;jj<K;jj++){
      sum+=Socres[knn_info(ii,jj)];}
    Res[ii]=sum/K;}
  return Res;
}
