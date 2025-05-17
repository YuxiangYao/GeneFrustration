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
  void FirstCandidate();
  int PseudoEnergy();
  // OE/KD reshape:
  void OEKD_RandomVec_PointSS(IntegerVector oriSS,IntegerVector gene,IntegerVector vals);
  //int OEKD_GotoStable(IntegerVector gene,int n_control,std::set<int> &gene_code,
  //    int *EachID, int *EachTempEng,IntegerVector markxx,int nx);
  Rcpp::List OEKD_GotoStable(IntegerVector gene,int n_control,std::set<int> &gene_code,
      int *EachID, int *EachTempEng,IntegerVector markxx,int nx);
  int TransitionPressure_quick(IntegerVector ControlGene, IntegerVector ControlVals,
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

void boolnetsys::FirstCandidate(){
  int ii,*topo=nullptr,*wigt=nullptr,*sss_vec=sss.data();
  for(ii=0;ii<size;ii++){
    if(TF_ind[ii]>0){
      topo=TF_topo[ii].data();
      wigt=TF_wigt[ii].data();
      if(sss[ii]!=WeightSumLogi(topo,wigt,sss_vec,TF_ind[ii],ii)){
        candidate.insert(ii);}}}// Node ii can change its state instead of having parents.
}

int boolnetsys::PseudoEnergy(){
  int ii,jj,sum,tmpr;
  std::vector<int> sss_spin=sss;
  std::vector<int> tmpw,tmpt;
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

void boolnetsys::OEKD_RandomVec_PointSS(IntegerVector oriSS,IntegerVector gene,IntegerVector vals){
  int ii,g_len=gene.size();
  for(ii=0;ii<size;ii++){
    sss[ii]=oriSS[ii];}
  for(ii=0;ii<g_len;ii++){// Set controller factors with configurated values.
    sss[gene[ii]-1]=vals[ii];}// Watch out the index between R and C++.
} 

Rcpp::List boolnetsys::OEKD_GotoStable(IntegerVector gene,int n_control,std::set<int> &gene_code, int *EachID, int *EachTempEng, IntegerVector markxx,int nx){
  int step=0,tmp_id,MAXSTEP=100000,index;
  int *topo=nullptr,*wigt=nullptr,*sss_vec=sss.data();
  // Mostly same as normal stable state finder, but remove check the pointed gene.
  auto flag=candidate.begin();
  for(const auto& xx : gene){// Remove all controllers.
    flag=candidate.find(xx-1);// !!Note: index C++ vs R
    if(flag!=candidate.end()){candidate.erase(xx-1);}}
  // Note: may reach a stable state.
  Rcpp::IntegerMatrix transistMaitrx(100,sss.size());
  Rcpp::List rreess;
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
    tmp_id=0;
    for(int pp=0;pp<nx;pp++){
      tmp_id+=(candidate.find(markxx[pp])!=candidate.end());// Find whether "xx" in the candidate Set.
    }
    EachID[step]=tmp_id;
    EachTempEng[step]=boolnetsys::PseudoEnergy();
    for(int qq=0; qq<sss.size(); qq++){
      transistMaitrx(step,qq)=sss[qq];}
    step++;}
  rreess.push_back(step,"steps");
  rreess.push_back(transistMaitrx,"TemporarySS_Mat");
  return rreess;
}

// Class-object: one transient state that has the tendency.
int boolnetsys::TransitionPressure_quick(IntegerVector ControlGene, IntegerVector ControlVals,std::set<int> &gene_code){
  int ii, jj, sum, pressure;
  std::vector<int> news=sss;
  for(ii=0;ii<ControlGene.size();ii++){// R vs C++ index
    news[ControlGene[ii]-1]=ControlVals[ii];}
  int *topo=nullptr,*wigt=nullptr;
  pressure=0;
  for(ii=0;ii<size;ii++){
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
        else ;
      }
  }
  return pressure;
}

int PhenotypicScore(std::vector<int> &ss,IntegerVector &mark_P,int mp,IntegerVector &mark_N,int mn){
  int ii,sum=0;
  for(ii=0;ii<mp;ii++)sum+=ss[mark_P[ii]];
  for(ii=0;ii<mn;ii++)sum-=ss[mark_N[ii]];
  return sum;
}

// Combined genes set 
// [[Rcpp::export]]
List Pressure2Potential(IntegerMatrix AdjMat,IntegerMatrix OriStabSta,List MarkerList,
  IntegerVector ControlSet, IntegerVector ControlVal,int RandomSeed){
  mt.seed(RandomSeed);int ii,jj,kk;
  int Ngene=OriStabSta.ncol(),Ncell=OriStabSta.nrow();
  Rcpp::IntegerVector Ori_PS(Ncell);
  Rcpp::IntegerVector Ori_EM(Ncell);
  Rcpp::IntegerVector OriEng(Ncell);
  Rcpp::IntegerVector GeneDiff(Ncell);
  Rcpp::IntegerVector NewEng(Ncell);
  Rcpp::IntegerVector New_PS(Ncell);
  Rcpp::IntegerVector New_EM(Ncell);
  Rcpp::IntegerVector Pressure(Ncell);
  Rcpp::IntegerVector NewTime(Ncell);
  Rcpp::IntegerMatrix NewState(Ncell,Ngene);
  Rcpp::IntegerMatrix Each_Update(Ncell,100);// 100 bits enough for recording the IDs.
  Rcpp::IntegerMatrix Each_TempEng(Ncell,100);// 100 bits enough for recording the TempEng.
  int TempID[100];
  int TempEng[100];
  Rcpp::List temp_List, all_process;
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
    OriEng[ii]=InducedCell.PseudoEnergy();
    Ori_PS[ii]=PhenotypicScore(old_sss,mark_p,Np,mark_s,Ns);
    Ori_EM[ii]=PhenotypicScore(old_sss,mark_e,Ne,mark_m,Nm);
    // Induced experiments;
    Pressure[ii]=InducedCell.TransitionPressure_quick(ControlSet,ControlVal,OEKDers);
    //for(int xx=0;xx<NewTime[ii];xx++){
    //   Each_Update(ii,xx)=TempID[xx];}
    InducedCell.OEKD_RandomVec_PointSS(OriStabSta.row(ii),ControlSet,ControlVal);
    InducedCell.candidate.clear();
    InducedCell.FirstCandidate();
    temp_List=InducedCell.OEKD_GotoStable(ControlSet,n_control,OEKDers,
        TempID,TempEng,mark_ps,Nps);
    NewTime[ii]=temp_List[0];
    all_process.push_back(temp_List[1]);
    for(int xx=0;xx<100;xx++){
      Each_Update(ii,xx)=-1;
      Each_TempEng(ii,xx)=1000;}
    for(int xx=0;xx<NewTime[ii];xx++){// Only record valid step;
       Each_Update(ii,xx)=TempID[xx];
       Each_TempEng(ii,xx)=TempEng[xx];}
    GeneDiff[ii]=DifferSets(InducedCell.sss,old_sss,Ngene);
    NewEng[ii]=InducedCell.PseudoEnergy();
    New_PS[ii]=PhenotypicScore(InducedCell.sss,mark_p,Np,mark_s,Ns);
    New_EM[ii]=PhenotypicScore(InducedCell.sss,mark_e,Ne,mark_m,Nm);
    for(int xx=0;xx<Ngene;xx++){
        NewState(ii,xx)=InducedCell.sss[xx];}
  }
  Rcpp::List ResList;
  ResList.push_back(Ori_PS,"Ori_PS");
  ResList.push_back(Ori_EM,"Ori_EM");
  ResList.push_back(OriEng,"OriEng");
  ResList.push_back(GeneDiff,"GeneDiff");
  ResList.push_back(NewEng,"NewEng");
  ResList.push_back(New_PS,"New_PS");
  ResList.push_back(New_EM,"New_EM");
  ResList.push_back(Pressure,"Pressure");
  ResList.push_back(NewTime,"NewTime");
  ResList.push_back(Each_Update,"EachUpdate");
  ResList.push_back(Each_TempEng,"EachTempEng");
  ResList.push_back(NewState,"NewState");
  ResList.push_back(all_process,"TransisMat");
  return ResList;
}

// Code's over!!