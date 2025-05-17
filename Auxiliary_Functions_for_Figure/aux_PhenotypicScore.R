# = = = = = = = = = = = = = = = = = = = = = = #
# Auxiliary function: phenotypic scores of SS #
# = = = = = = = = = = = = = = = = = = = = = = #

# Sample Matrix: row/col is sample/gene
# Return sample's EMT score.
EMT_score<-function(Samples){
  emt.ips=c('Cdh1','Tcf3','Ovol2','miR200','miR34');
  emt.mef=c('Vim','Cdh2','Zeb1','Zeb2','Foxc2','Snai1','Snai2','Twist1','Twist2','miR9');
  scores=rowSums(Samples[,intersect(emt.ips,colnames(Samples))])-
    rowSums(Samples[,intersect(emt.mef,colnames(Samples))]);
  return(scores);
}
# Return sample's SPT score.
SPT_score<-function(Samples){
  spt.ips=c('Esrrb','Nanog',"Klf4","Pou5f1","Sox2","Myc","Sall4")
  spt.mef=c('Jun','Gata6','Gata4','Bmp2','Bmp4','Cebpa',"Sox17","Sox7","Cdx2",'Foxa1','Foxa2');
  scores=rowSums(Samples[,intersect(spt.ips,colnames(Samples))])-
    rowSums(Samples[,intersect(spt.mef,colnames(Samples))]);
  return(scores);
}
# Code is over.