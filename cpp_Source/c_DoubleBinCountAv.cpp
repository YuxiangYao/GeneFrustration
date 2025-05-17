// C++ for double-bin for plotting energy geom_tile:

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix DoubleBins(IntegerVector x,IntegerVector y, NumericVector Eng, int MaxX, int MaxY){
  int ii,jj;
  NumericMatrix resss(MaxX,MaxY);
  NumericMatrix count(MaxX,MaxY);
  for(ii=0; ii<Eng.size(); ++ii){
    resss(x[ii],y[ii])=resss(x[ii],y[ii])+Eng[ii];
    count(x[ii],y[ii])=count(x[ii],y[ii])+1;}
  for(ii=0; ii<MaxX; ++ii){
    for(jj=0; jj<MaxY; ++jj){
      if(count(ii,jj)>0){
        resss(ii,jj)=resss(ii,jj)/count(ii,jj);}
      else {
        resss(ii,jj)=0;}}}
  return resss;
}

// [[Rcpp::export]]
NumericMatrix DoubleCount(IntegerVector x,IntegerVector y, NumericVector Eng, int MaxX, int MaxY){
  int ii,jj;
  NumericMatrix count(MaxX,MaxY);
  for(ii=0; ii<Eng.size(); ++ii){
    count(x[ii],y[ii])=count(x[ii],y[ii])+1;}
  return count;
}