# Prepare for 10^6 random SM-like cells.

# Load data & necessary function:
adjmat=readRDS("./Data_Used/AdjacentMatrix.Matrix.rds");
Rcpp::sourceCpp("./cpp_Source/c_GotoStableDistribution.cpp");
source("./Auxiliary_Functions_for_Figure/aux_PhenotypicScore.R");

# Quickly generate four classical cell populations.
SimulatedCharacteristic<-function(adjmat,RandomSeed){
  tmp=B_ScanningAttrctors_Fast(adjmat, 100000L, RandomSeed);
  sss=tmp$StateMatrix;
  emt=EMT_score(sss);
  spt=SPT_score(sss);
  id.sm=(spt<=-5&emt<=-4);# Set MEF-like cells.
  sm=sss[id.sm,];
  ii=1;
  while(nrow(sm)<1e5){
    tmp=B_ScanningAttrctors_Fast(adjmat,20000L,RandomSeed+ii);
    sss=tmp$StateMatrix;
    emt=EMT_score(sss);
    spt=SPT_score(sss);
    id.sm=(spt<=-5&emt<=-4);
    sm=rbind(sm,sss[id.sm,]);
    ii=ii+1;}
  sm=sm[!duplicated(sm),];
  sm=sm[1:min(c(1e5,nrow(sm))),];
  return(sm);
}

# Code is over.