# Load original energy landscape as benchmark>>>
pc.index=readRDS("../Intermediate_Result/PCA_index_coor.rds");
ori.eng=readRDS("../Intermediate_Result/Step_Energy.rds");
nondup.id=readRDS("../Intermediate_Result/nonduplicate_index.rds")
ori.eng=ori.eng[nondup.id,2];
sss=readRDS("../Intermediate_Result/SSS_origin_nodup.rds");
  
CalculateBinnedPC_Coord_LocalFrust_Ori<-function(Specificgenes,AdjMat,TransisMat,
                                                   Group_Av=FALSE){
    #Convertas spin-like:
    sp_mat=TransisMat;
    sp_mat[sp_mat<1]=-1;
    # Set return matrix.
    ResLocFru=matrix(NA,nrow(sp_mat),length(Specificgenes));
    colnames(ResLocFru)=Specificgenes;
    ResLocFru=-sp_mat%*%t(AdjMat[Specificgenes,])*sp_mat[,Specificgenes];
    # Load rotation matrix for embedding phase.
    cat("Calculate PCs.\n");
    rotation=readRDS("../Intermediate_Result/PCA_index_rotation_origin.rds");
    transient_pc_success=TransisMat%*%rotation[,1:5];
    if(Group_Av){
      ResLocFru=apply(ResLocFru,1,mean);}
    return(list(pcs=transient_pc_success, tmp_eng=ResLocFru));
  }
# Only HVGs.
HVGs_LocalFru_Ori=CalculateBinnedPC_Coord_LocalFrust_Ori(gene_hig, adjmat, 
                                                           sss,TRUE);
# Only LVGs.
LVGs_LocalFru_Ori=CalculateBinnedPC_Coord_LocalFrust_Ori(gene_low, adjmat, 
                                                           sss, TRUE);
# All genes as benchmark.
Alls_LocalFru_Ori=CalculateBinnedPC_Coord_LocalFrust_Ori(colnames(adjmat), adjmat, 
                                                           sss, TRUE);
saveRDS(Alls_LocalFru_Ori,"../Intermediate_Result/OriginalPhoneLandLocalFru.rds")