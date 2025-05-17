// Provided code adheres to C++ specifications rather than R!
#include <stdio.h>
#include <stdlib.h>// for itoa
#include <iostream>
#include <cstdio>
#include <vector>
#include <set>
#include <iterator>
#include <algorithm>
#include <random>
#include <cmath>
#include <Rcpp.h>// The header offers a seamless integration of R and C++.
using namespace Rcpp;

std::mt19937 mt(12345);
std::uniform_int_distribution<int> cpp_rbinary(0,1);
std::uniform_real_distribution<double> cpp_runif(0,1);
// Sign function.
int Sign(int x){
    if(x>0)return 1;
    else if(0==x)return 0;
    else return -1;
}
// Update judgements (Eq.1 in main text).
int WeightSumLogi(int *topo,int *wigt,int *vecs,int degs,int code){
    int ii,sum=0;
    for(ii=0;ii<degs;ii++){
        sum+=vecs[topo[ii]]*wigt[ii];}
    if(0==sum){// Keep old value.
        return vecs[code];}
    else {
        return (sum>0);}
}
// Different bits between two binary vectors.
int DifferSets(std::vector<int> &a,std::vector<int> &b,int &num){
    int sum=0;
    for(int ii=0;ii<num;ii++){
        sum+=(a[ii]!=b[ii]);}
    return sum;
}
// Deleta one element from Set-class.
void DeleteElementSet(std::set<int> *A,int element){
    auto it=A->find(element);
    if(it != A->end()){// Avoid element not existed in Set a. 
        A->erase(element);}
}

// Class for Boolean network system
class boolnetsys {
    public:
        int size;
        std::vector<int> TF_ind;
        //std::vector<int> Ep_ind;
        std::vector<int> sss;
        std::vector<std::vector<int>> TF_topo;// TF linking
        std::vector<std::vector<int>> TF_wigt;// TF's wight 
        //std::vector<std::vector<int>> Ep_topo;// Epi linking
        //std::vector<std::vector<int>> Ep_wigt;// Epi's wight
        std::vector<std::vector<int>> TF_out;// Pointed node's TF children.
        //std::vector<std::vector<int>> Ep_out;// Pointed node's EPi children.
        std::set<int> candidate;
    // Functions:
        void Initialization(IntegerMatrix AdjMat);
    // Control:
        int OEKD_control(int IDs,int OE_KD);
};
void boolnetsys::Initialization(IntegerMatrix AdjMat){
    int ii,jj,kk,inds,n_gene=AdjMat.nrow();
    // Set numbers.
    size=n_gene;
    sss.resize(n_gene);
    TF_ind.resize(n_gene);
    TF_topo.resize(n_gene);
    TF_wigt.resize(n_gene);
    TF_out.resize(n_gene);
    std::vector<int> tmpt(n_gene);// Which gene.
    std::vector<int> tmpw(n_gene);// What's weight.
    // Set node_ii all input links.
    for(ii=0;ii<n_gene;ii++){
        inds=0;tmpt.clear();tmpw.clear();
        for(jj=0;jj<n_gene;jj++){
            if(0!=AdjMat(ii,jj)){// Non-zero, have regulated links.
                inds++;
                tmpt.push_back(jj);
                tmpw.push_back(AdjMat(ii,jj));}}
        TF_topo[ii]=tmpt;
        TF_wigt[ii]=tmpw;
        TF_ind[ii]=inds;}
    // Set node_ii all output links (for fast algorithm).
    for(ii=0;ii<n_gene;ii++){
        tmpt.clear();
        for(jj=0;jj<n_gene;jj++){
            if(0!=AdjMat(jj,ii)){// Non-zero, have child nodes.
                tmpt.push_back(jj);}}
        TF_out[ii]=tmpt;}
}

int boolnetsys::OEKD_control(int IDs,int OE_KD){// Default: state vector has reached attractors.
    int step=0,tmp_id,index,MAXSTEP=100000;
    int *topo=nullptr,*wigt=nullptr,*sss_vec=sss.data();
    auto flag=candidate.begin();
    candidate.clear();// Clear set of candidate.
    sss[IDs]=OE_KD;// Set fixed points.
    for(const auto& xx : TF_out[IDs]){
        topo=TF_topo[xx].data();
        wigt=TF_wigt[xx].data();
        if(sss[xx]!=WeightSumLogi(topo,wigt,sss_vec,TF_ind[xx],xx)){
            candidate.insert(xx);}
        else {
            if(candidate.find(xx)!=candidate.end()){// Is in the set. [candidate.count(xx) also okay]
                candidate.erase(xx);}}}
    candidate.erase(IDs);// Must remove the pointed gene (no need to change)
    while(candidate.size()>0&&step<MAXSTEP){
        flag=candidate.begin();
        tmp_id=(int)(cpp_runif(mt)*candidate.size());// Random choose a candidate.
        advance(flag,tmp_id);
        index=(*flag);
        sss[index]=(0==sss[index]);// Update: 1->0 or 0->1
        candidate.erase(index);// Remove this unstable element.
        for(const auto& xx : TF_out[index]){
            topo=TF_topo[xx].data();
            wigt=TF_wigt[xx].data();
            if(sss[xx]!=WeightSumLogi(topo,wigt,sss_vec,TF_ind[xx],xx)){
                candidate.insert(xx);}
            else {
                if(candidate.find(xx)!=candidate.end()){// Is in the set. [candidate.count(xx) also okay]
                    candidate.erase(xx);}}}
        candidate.erase(IDs);// Always remove the pointed gene.
        step++;}
    return step;
}

// [[Rcpp::export]]
List OEKD_Testing(IntegerMatrix AdjMat,
    IntegerMatrix SSmat,IntegerVector Markerlist,int RandomSeed){
  mt.seed(RandomSeed);
  int ii,jj,kk,tmp,ncell=SSmat.nrow(),ngene=SSmat.ncol(),nmarker=Markerlist.size();
  Rcpp::IntegerMatrix oe_diff(ncell,ngene),oe_time(ncell,ngene);
  Rcpp::IntegerMatrix kd_diff(ncell,ngene),kd_time(ncell,ngene);
  boolnetsys celloekd;
  celloekd.Initialization(AdjMat);
  std::vector<int> old_sss(ngene);
  Rcpp::IntegerMatrix xoe(ncell,nmarker*ngene),xkd(ncell,nmarker*ngene);
  std::vector<int> marker_id;
  // Obtain the marker's ID.
  for(ii=0;ii<nmarker;ii++){
    marker_id.push_back(Markerlist[ii]-1);}
  // Check each 
  for(ii=0;ii<ncell;ii++){
    for(kk=0;kk<ngene;kk++){// Copy stable state.
      old_sss[kk]=SSmat(ii,kk);}
    for(jj=0;jj<ngene;jj++){
      // for OE
      celloekd.sss=old_sss;
      oe_time(ii,jj)=celloekd.OEKD_control(jj,1);
      oe_diff(ii,jj)=DifferSets(celloekd.sss,old_sss,ngene);
      tmp=jj*nmarker;
      for(kk=0;kk<nmarker;kk++){
        xoe(ii,tmp+kk)=celloekd.sss[marker_id[kk]];}
      // for KD
      celloekd.sss=old_sss;
      kd_time(ii,jj)=celloekd.OEKD_control(jj,0);
      kd_diff(ii,jj)=DifferSets(celloekd.sss,old_sss,ngene);
      tmp=jj*nmarker;
      for(kk=0;kk<nmarker;kk++){
        xkd(ii,tmp+kk)=celloekd.sss[marker_id[kk]];}}}
  Rcpp::List ResList;
  ResList.push_back(oe_diff,"OE_Diff");
  ResList.push_back(oe_time,"OE_Time");
  ResList.push_back(kd_diff,"KD_Diff");
  ResList.push_back(kd_diff,"KD_Time");
  ResList.push_back(xoe,"OE_marker");
  ResList.push_back(xkd,"KD_marker");
  return ResList;
}

// Code is over.