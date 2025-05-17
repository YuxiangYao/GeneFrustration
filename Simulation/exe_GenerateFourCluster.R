# Execute simulation: Generate four clusters of significantly phenotypic cells.

# Load data & necessary function:
adjmat=readRDS("./Data_Used/AdjacentMatrix.Matrix.rds");
Rcpp::sourceCpp("./cpp_Source/c_GotoStableDistribution.cpp");

# Quickly generate four classical cell populations.
SimulatedCharacteristic<-function(adjmat,RandomSeed){
  tmp=B_ScanningAttrctors_Fast(adjmat, 1000000L, RandomSeed);
  sss=tmp$StateMatrix;
  emt=EMT_score(sss);       spt=SPT_score(sss);
  id.pe=(spt>=3 &emt>=3);   id.pm=(spt>=3 &emt<=-4);
  id.se=(spt<=-5&emt>=3);   id.sm=(spt<=-5&emt<=-4);
  pe=sss[id.pe,]; pm=sss[id.pm,];
  se=sss[id.se,]; sm=sss[id.sm,];
  ii=1;
  while(nrow(pe)<1e6||nrow(pm)<1e6||nrow(se)<1e6||nrow(sm)<1e6){
    tmp=B_ScanningAttrctors_Fast(adjmat,200000L,RandomSeed+ii);
    sss=tmp$StateMatrix;
    emt=EMT_score(sss);       spt=SPT_score(sss);
    id.pe=(spt>=3 &emt>=3);   id.pm=(spt>=3 &emt<=-4);
    id.se=(spt<=-5&emt>=3);   id.sm=(spt<=-5&emt<=-4);
    pe=rbind(pe,sss[id.pe,]);
    pm=rbind(pm,sss[id.pm,]);
    se=rbind(se,sss[id.se,]);
    sm=rbind(sm,sss[id.sm,]);#cat("Repeat times: ",ii,"\n");
    ii=ii+1;}
  pe=pe[!duplicated(pe),];pe=pe[1:min(c(1e6,nrow(pe))),];
  pm=pm[!duplicated(pm),];pm=pm[1:min(c(1e6,nrow(pm))),];
  se=se[!duplicated(se),];se=se[1:min(c(1e6,nrow(se))),];
  sm=sm[!duplicated(sm),];sm=sm[1:min(c(1e6,nrow(sm))),];
  res=list(pe=pe,pm=pm,se=se,sm=sm);
  return(res);
}

# Generate some extreme cases (~15min)
ExtremeCase=SimulatedCharacteristic(adjmat=adjmat,239921L);
saveRDS(ExtremeCase, "./Intermediate_Result/ExtremeCase_FourClusters.rds", compress=FALSE);

rm(list=ls());gc();
# Code is over.