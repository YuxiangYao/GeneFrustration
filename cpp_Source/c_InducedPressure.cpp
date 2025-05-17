// Test the ΔPhenotypic score and frustration
// Date: 2024-06-13

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
#include <Rcpp.h>
using namespace Rcpp;

std::mt19937 mt(12345);
std::uniform_int_distribution<int> cpp_rbinary(0,1);
std::uniform_real_distribution<double> cpp_runif(0,1);

int Sign(int x){
    if(x>0)return 1;
    else if(0==x)return 0;
    else return -1;
}
int WeightSumLogi(int *topo,int *wigt,int *vecs,int degs,int code){
    int ii,sum=0;
    for(ii=0;ii<degs;ii++){
        sum+=vecs[topo[ii]]*wigt[ii];}
    if(0==sum){// Keep old value.
        return vecs[code];}
    else {
        return (sum>0);}
}
void DeleteElementSet(std::set<int> *A,int element){
    auto it=A->find(element);
    if (it != A->end()){// Avoid element not existed in Set a. 
    A->erase(element);}
}

int DifferSets(std::vector<int> &a,std::vector<int> &b,int &num){
    int sum=0;
    for(int ii=0;ii<num;ii++){
        sum+=(a[ii]!=b[ii]);}
    return sum;
}

class boolnetsys {
public:
    int size;
    std::vector<int> TF_ind;
    std::vector<int> Ep_ind;
    std::vector<int> sss;
    std::vector<std::vector<int>> TF_topo;// TF linking
    std::vector<std::vector<int>> TF_wigt;// TF's wight 
    std::vector<std::vector<int>> TF_out;// Pointed node's TF children.
    std::vector<std::vector<int>> Ep_out;// Pointed node's EPi children.
    std::set<int> candidate;
    // Functions:
    void Initialization(IntegerMatrix AdjMat);
    void RandomVec();
    void FirstCandidate();
    int TransitionPressure_quick2(IntegerVector ControlGene, IntegerVector ControlVals,
        std::set<int> &gene_code);
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

void boolnetsys::FirstCandidate(){
    int ii,*topo=nullptr,*wigt=nullptr,*sss_vec=sss.data();
    for(ii=0;ii<size;ii++){
        if(TF_ind[ii]>0){
            topo=TF_topo[ii].data();
            wigt=TF_wigt[ii].data();
        if(sss[ii]!=WeightSumLogi(topo,wigt,sss_vec,TF_ind[ii],ii)){
            candidate.insert(ii);}}}// Node ii can change its state instead of having parents.
}

// Class-object: one transient state that has the tendency.
int boolnetsys::TransitionPressure_quick2(IntegerVector ControlGene, IntegerVector ControlVals,std::set<int> &gene_code){
    int ii, jj, sum, pressure;
    std::vector<int> news=sss;
    for(ii=0;ii<ControlGene.size();ii++){// R vs C++ index
        news[ControlGene[ii]-1]=ControlVals[ii];}
    int *topo=nullptr,*wigt=nullptr;
    pressure=0;
    for(ii=0;ii<size;ii++){// All gene checking.
      // Only check non-controller & Only check the in-deg >0 
      //if((TF_ind[ii]>0)&&(gene_code.find(ControlGene[ii]-1)==gene_code.end())){
        if((TF_ind[ii]>0)&&(gene_code.find(ii)==gene_code.end())){
            sum=0;
            topo=TF_topo[ii].data();
            wigt=TF_wigt[ii].data();
            for(jj=0;jj<TF_ind[ii];jj++){
                sum+=news[topo[jj]]*wigt[jj];}
            // Judge: whether the pointed genes should be changed!
            if(sum>0&&sss[ii]==0){
                pressure+=1;}
            else if(sum<0&&sss[ii]==1){
                pressure+=1;}
            else ;}}
    return pressure;
}


// [[Rcpp::export]]
IntegerVector Pressure2Potential_Landscape(IntegerMatrix AdjMat,IntegerMatrix OriStabSta,List MarkerList,
    IntegerVector ControlSet, IntegerVector ControlVal,int RandomSeed){
    mt.seed(RandomSeed);int ii,jj,kk;
    int Ngene=OriStabSta.ncol(),Ncell=OriStabSta.nrow();
    Rcpp::IntegerVector Pressure(Ncell);
    // Set class of Boolean network system.
    boolnetsys InducedCell;
    InducedCell.Initialization(AdjMat);
    std::vector<int> old_sss(Ngene);
    std::set<int> OEKDers;
    int n_control=ControlSet.size();
    for(int ii=0;ii<n_control;ii++){
        OEKDers.insert(ControlSet[ii]-1);}// [!!Note] Index: R vs C++
    // Index of markers ( C++'s code, start for 0 ).
    Rcpp::IntegerVector mark_p=MarkerList[0];int Np=mark_p.size();// P-marker
    Rcpp::IntegerVector mark_s=MarkerList[1];int Ns=mark_s.size();// S-marker
    Rcpp::IntegerVector mark_e=MarkerList[2];int Ne=mark_e.size();// E-marker
    Rcpp::IntegerVector mark_m=MarkerList[3];int Nm=mark_m.size();// M-marker
    Rcpp::IntegerVector mark_ps(Np+Ns);int Nps=Np+Ns;
    for(ii=0;ii<Np;ii++)mark_ps[ii]=mark_p[ii];
    for(ii=0;ii<Ns;ii++)mark_ps[ii+Np]=mark_s[ii];
    // Check each gene's case.
    for(ii=0;ii<Ncell;ii++){
        // Copy stable state.
        for(kk=0;kk<Ngene;kk++){old_sss[kk]=OriStabSta(ii,kk);}
        InducedCell.sss=old_sss;
        // Induced experiments;
        Pressure[ii]=InducedCell.TransitionPressure_quick2(ControlSet,ControlVal,OEKDers);
    }
    return Pressure;
}

// Code's over!!