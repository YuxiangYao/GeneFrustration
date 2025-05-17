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
        std::vector<int> sss;
        std::vector<std::vector<int>> TF_topo;// TF linking
        std::vector<std::vector<int>> TF_wigt;// TF's wight 
        std::vector<std::vector<int>> TF_out;// Pointed node's TF children.
        std::set<int> candidate;
    // Functions:
        void Initialization(IntegerMatrix AdjMat);
        void RandomVec();
        void FirstCandidate();
        int PseudoEnergy();
        void OEKD_RandomVec_PointSS(IntegerVector oriSS,IntegerVector gene,IntegerVector vals);
        int OEKD_GotoStable(IntegerVector gene,int n_control,std::set<int> &gene_code);
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
void boolnetsys::RandomVec(){
    for(int ii=0;ii<size;ii++){
        sss[ii]=cpp_rbinary(mt);}
}
int boolnetsys::PseudoEnergy(){
    int ii,jj,sum,tmpr;
    std::vector<int> tmpw,tmpt,sss_spin=sss;
    for(ii=0;ii<sss_spin.size();ii++){
        if(0==sss_spin[ii]){sss_spin[ii]=-1;}}
    sum=0;
    for(ii=0;ii<size;ii++){
        tmpw=TF_wigt[ii];
        tmpt=TF_topo[ii];
        tmpr=sss_spin[ii];
        for(jj=0;jj<TF_ind[ii];jj++){// Get in-degree.
            sum+=tmpw[jj]*tmpr*sss_spin[tmpt[jj]];}}
    return -sum;
}
void boolnetsys::FirstCandidate(){
    int ii,*topo=nullptr,*wigt=nullptr,*sss_vec=sss.data();
    for(ii=0;ii<size;ii++){
        if(TF_ind[ii]>0){// Only need to check non-zero in-degree.
            topo=TF_topo[ii].data();
            wigt=TF_wigt[ii].data();
            // Whether ii-th node update regulated by its upstrem?
            if(sss[ii]!=WeightSumLogi(topo,wigt,sss_vec,TF_ind[ii],ii)){
                candidate.insert(ii);}}}
}

void boolnetsys::OEKD_RandomVec_PointSS(IntegerVector oriSS,IntegerVector gene,IntegerVector vals){
  int ii,g_len=gene.size();
  for(ii=0;ii<size;ii++){
    sss[ii]=oriSS[ii];}
  for(ii=0;ii<g_len;ii++){// Set controller factors with configurated values.
    sss[gene[ii]-1]=vals[ii];}// Watch out the index between R and C++.
} 

int boolnetsys::OEKD_GotoStable(IntegerVector gene,int n_control,std::set<int> &gene_code){
  int step=0,tmp_id,MAXSTEP=100000,index;
  int *topo=nullptr,*wigt=nullptr,*sss_vec=sss.data();
  // std::vector<int> tmp_slot(size);
  // Mostly same as normal stable state finder, but remove check the pointed gene.
  auto flag=candidate.begin();
  for(const auto& xx : gene){// Remove all controllers.
    flag=candidate.find(xx-1);// !!Note: index C++ vs R
    if(flag!=candidate.end()){candidate.erase(xx-1);}}
  // Note: may reach a stable state.
  while(candidate.size()>0&&step<MAXSTEP){
    flag=candidate.begin();
    tmp_id=(int)(cpp_runif(mt)*candidate.size());// Random choose a candidate.
    advance(flag,tmp_id);
    index=(*flag);
    sss[index]=(0==sss[index]);// Update: 1->0 or 0->1
    candidate.erase(index);// Remove this unstable element.
    for(const auto& xx : TF_out[index]){
      flag=gene_code.find(xx);
      if(flag==gene_code.end()){// Check if it is contained in controller set.
        topo=TF_topo[xx].data();
        wigt=TF_wigt[xx].data();
        if(sss[xx]!=WeightSumLogi(topo,wigt,sss_vec,TF_ind[xx],xx)){
          candidate.insert(xx);}
        else {
          if(candidate.find(xx)!=candidate.end()){// Is in the set. [candidate.count(xx) also okay]
            candidate.erase(xx);}}
      }
    }
  step++;}
  return step;
}

// [[Rcpp::export]]
List B_ScanningAttrctors_oekd_pointed(IntegerMatrix AdjMat,
    IntegerMatrix OriStabSta,int RandomSeed,IntegerVector Indexs,IntegerVector Values){
  mt.seed(RandomSeed);
  int Times=OriStabSta.nrow();
  Rcpp::IntegerMatrix ResMat(Times,AdjMat.nrow());// Times*nGene
  colnames(ResMat)=colnames(AdjMat);
  Rcpp::IntegerMatrix Steps(Times,2);// Each step of random vector fall in basin.
  boolnetsys attrpool;
  attrpool.Initialization(AdjMat);
  int n_control=Indexs.size();
  std::set<int> OEKDers;
  for(int ii=0;ii<n_control;ii++){
    OEKDers.insert(Indexs[ii]-1);}// [!!Note] Index: R vs C++
  for(int ii=0;ii<Times;ii++){
    attrpool.OEKD_RandomVec_PointSS(OriStabSta.row(ii),Indexs,Values);
    attrpool.candidate.clear();
    attrpool.FirstCandidate();
    Steps(ii,0)=attrpool.OEKD_GotoStable(Indexs,n_control,OEKDers);
    Steps(ii,1)=attrpool.PseudoEnergy();
    for(int jj=0;jj<attrpool.size;jj++){
      ResMat(ii,jj)=attrpool.sss[jj];}}
  Rcpp::List res_list=Rcpp::List::create(
    Named("StateMatrix")=ResMat,Named("Steps")=Steps);
  return res_list;
}

// Code is over.