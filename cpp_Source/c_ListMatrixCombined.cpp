// C++ for combining matrices within a list:

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix aftercombined_transisMat(List sublist_matrix){
  int ii,jj;
  IntegerMatrix tmp=sublist_matrix[1];
  int n_matrix=sublist_matrix.size(),
    n_row=tmp.nrow(), n_col=tmp.ncol();
  NumericMatrix resss(n_matrix*n_row, n_col);
  
  int index=0;
  for(ii=0; ii<n_matrix; ++ii){
    IntegerMatrix tmp2=sublist_matrix[ii];
    for(jj=0; jj<n_row; ++jj){
      resss.row(index)=tmp2.row(jj);
      index++;}}
  return resss;
}