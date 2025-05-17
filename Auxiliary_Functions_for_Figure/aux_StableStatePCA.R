# = = = = = = = = = = = = = = = = = = = = = = = = #
# Auxiliary function: Quick PCA of stable points  #
# = = = = = = = = = = = = = = = = = = = = = = = = #

PCA_index_StableState<-function(StablePointMatrix){
  # Row: sample, Col: gene.
  # Get PCA index.
  #obtain.pr=prcomp(sss,center=FALSE,scale.=FALSE,rank.=10);# Slow!
  obtain.pr=gmodels::fast.prcomp(StablePointMatrix,center=FALSE,scale.=FALSE);
  # PC score dimension (PC 1~5)
  pca.index=list(dims=obtain.pr$x[,1:5],rotation=obtain.pr$rotation);
  return(pca.index);
}
