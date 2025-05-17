# = = = = = = = = = = = = = = = = = = = = #
# Auxiliary function: filter valid genes  # 
#                   & binarize the data   #
# = = = = = = = = = = = = = = = = = = = = #

# Filter genes.
{GeneList=c(
  "Sox2","Sall4","Esrrb","Nanog","Dppa2","Pou5f1",
  "Klf4","Jun","Jdp2","Tet2","Tgfb1","Myc",
  "miR200","miR141","Bmp4","Rif1","Tcl1","Sp1",
  "E2f1","Sumo1","Pias1","E2f4","Cebpa","Sp3",
  "Sox17","Bmp2","Lef1","Cdx2","Zic3","Smad2",
  "Smad4","Tcf7l1","Gata6","Foxd3","Prdm14","Nodal",
  "Zeb2","Cdk9","Gata4","Ezh2","Snai1","Zeb1",
  "Tfcp2l1","Fgf4","miR145","Trp53","Dnmt3a","Glis1",
  "Lif","Mapk1","Gsk3b","Ctnnb1","Mbd2","miR101",
  "miR205","miR34","Snai2","Twist1","Twist2","Gsc",
  "Ovol2","Grhl2","Cldn7","Cdh1","Trp63","Nkx3-1",
  "Akt1","Foxa2","Foxa1","Sox7","Mbd3","Ino80",
  "Zfp42","TeKlf4","Suv39h1","Suv39h2","Ep300","Sap30",
  "Setdb1","Ehmt2","Parp1","Kdm6b","Dppa4","Kdm3a",
  "Kdm4c","Zbtb7a","Runx1","Fgfr2");}
Keep.Gene<-function(GeneVec,GeneList){
  keepname=intersect(GeneVec,GeneList);
  return(keepname);
}
Fill.Gene<-function(GeneVec,Nsam,random=NA,GeneList){
  remains=setdiff(GeneList,GeneVec);
  res=matrix(sample.int(2,length(remains)*Nsam,replace=T)-1,
             Nsam,length(remains),dimnames=list(NULL,remains));
  return(res);
}
Stand.Gene<-function(GeneMatrix,Baseline){
  #Row-gene,Col-sample;length(Baseline)=sample
  res=t(GeneMatrix)/Baseline;
  return(t(res));
}
Intersect.Gene<-function(GeneMatrix,GeneList){
  res=GeneMatrix;#Row-gene,Col-sample;length(Baseline)=sample
  index=intersect(GeneList,rownames(res));
  res=res[index,];
  return(res);
}

# Binarize gene's expression.
library(mclust);
Binarization_Kmeans<-function(Vec){
  types=kmeans(Vec,centers=2)$cluster;
  if(sum(types==1)==0||sum(types==2)==0){
    types=rep(0,length(Vec));
  }else {# if the label inverse?
    if(mean(Vec[types==1])>mean(Vec[types==2])){
      types=2-types;
    }else {
      types=types-1;}}
  return(as.integer(types));
}
Binarization_GMM<-function(Vec){
  types=mclust::Mclust(Vec,G=2,verbose=FALSE)$classification
  if(sum(types==1)==0||sum(types==2)==0){
    types=rep(0,length(Vec));
  }else {# if the label inverse?
    if(mean(Vec[types==1])>mean(Vec[types==2])){
      types=2-types;
    }else {
      types=types-1;}}
  return(as.integer(types));
}
Binarization_DropOut<-function(Vec){# Cluster for zero drop-out type
  index_n0=Vec>0;# Only consider non-zero values.
  tmp_vec=Vec[index_n0];
  types=mclust::Mclust(tmp_vec,2,verbose = FALSE)$classification;
  if(sum(types==1)==0||sum(types==2)==0){
    types=rep(1,length(Vec));
  }else {
    if(mean(tmp_vec[types==1])>mean(tmp_vec[types==2])){
      types=2-types;# if the label inverse?
    }else {
      types=types-1;}}
  Res=Vec;
  Res[index_n0]=types;
  return(as.integer(Res));
}
Binarization_Bool<-function(Vec){# Cluster for zero drop-out type
  res=Vec>0
  return(as.integer(res));
}

BinaryMatrixObtain<-function(ExpExpMat,GeneList,ScaleFactor=1000){
  hoskeep=as.vector(as.matrix(ExpExpMat["Gapdh",]));
  expdata=as.matrix(ExpExpMat[Keep.Gene(rownames(ExpExpMat),GeneList),]);
  expdata=Stand.Gene(expdata,hoskeep);
  expdata=as.matrix(apply(expdata*ScaleFactor,1,Binarization_Kmeans));
  rownames(expdata)=colnames(ExpExpMat);
  return(expdata);
}

