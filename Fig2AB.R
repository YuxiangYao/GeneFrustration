# Figure: Phenotypic distribution of random initial state
#   Stable in PC2 and PC3 via two different phenotypic measures: pluripotent vs somatic,
# epithelial vs mesenchymal. They also occupy different pseudo-energy regions.

# Complie quick smooth KNN cpp-based function.
Rcpp::sourceCpp("./cpp_Source/c_SmoothKNN.cpp");
# Load PCA function.
source("./Auxiliary_Functions_for_Figure/aux_StableStatePCA.R");
# Load functions of calcualting phenotypic score.
source("./Auxiliary_Functions_for_Figure/aux_PhenotypicScore.R");
# Load our default theme of ggplot object.
source("./Auxiliary_Functions_for_Figure/aux_ggplot2Theme.R");

# Load stable states (SSs).
sss=readRDS("./Intermediate_Result/SSS_origin_nodup.rds");
# Embed SSs into a PCA-space.
ind=PCA_index_StableState(sss);
# Obtain PCA coordinates.
pca_index=cbind.data.frame(ind$dims);
colnames(pca_index)=paste0("x",c(1:5));
# Highly recommend! for subsequent using.
saveRDS(ind$rotation, file="./Intermediate_Result/PCA_index_rotation_origin.rds", compress=FALSE);
saveRDS(pca_index, file="./Intermediate_Result/PCA_index_coor.rds", compress=FALSE);
# Get KNN-based neighbor indexs of PC2 and PC3 for smoothing scores.
all_knn=RcppHNSW::hnsw_knn(as.matrix(pca_index[,2:3]), k=5, distance="l2", n_threads=64);# Please note the core number here.

# Function: Draw the PCA-phonetypic scores figures.
fig_pca_score<-function(DataFrame,LowMidHigh=c("#2156A6","white","#EE312E")){
  fig_pca_s=ggplot()+
  geom_point(data=DataFrame,aes(x=x2,y=x3,col=pheno),size=0.1)+# Should keep same as L35 code.
  labs(x=paste0('PC',2),y=paste0('PC',3))+
  scale_color_gradient2(
    limits=c(min(score_smooth),max(score_smooth)),
    low=LowMidHigh[1], mid=LowMidHigh[2], high=LowMidHigh[3],
    midpoint=0.5*(min(score_smooth)+max(score_smooth)));
  fig_pca_s=Plot_ThemeConfiguration(fig_pca_s)+theme(
  axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank(),
  axis.text.x.bottom = element_blank(),axis.text.y.left = element_blank(),
  axis.title.x.bottom = element_blank(),axis.title.y.left = element_blank(),
  aspect.ratio=1,legend.position="right");
  return(fig_pca_s);
}

# Figure of phase diagram SPT-KNN-PC23 
score_smooth=SmoothKNN(all_knn[[1]]-1, SPT_score(sss), 5);# "-1" denotes index between C++ and R.
datas=cbind.data.frame(pca_index,pheno=score_smooth);
#datas=cbind.data.frame(pca_index,SPT_score(sss));# Non-smoothed one (Actually close to the smoothed one)
fig_pca_spt23=fig_pca_score(datas,c("#2156A6","white","#EE312E"));
ggsave(filename="./Figure_OriginalVersion/pc23_spt_knn_20231127.png",fig_pca_spt23,width=7.5,height=4.5,units="in");

# Figure of phase diagram EMT-KNN-PC23 (Fig. 1b/S1c)
score_smooth=SmoothKNN(all_knn[[1]]-1, EMT_score(sss), 5);
datas=cbind.data.frame(pca_index,pheno=score_smooth);
#datas=cbind.data.frame(pca_index,EMT_score(sss));# Non-smoothed one (Actually close to the smoothed one)
fig_pca_emt23=fig_pca_score(datas,c("#3B2E7E","white","#FFB600"));
ggsave(filename="./Figure_OriginalVersion/pc23_emt_knn_20231127.png",fig_pca_emt23,width=7.5,height=4.5,units="in");

rm(list=ls());gc();
# Code is over. 