# This code will generate figures of PCs on high/low cluster transient pseudo-energy
# on all successful/failed trajectories, and energy-landscape changed.
# Output 4 figures:
# TransientTrajectory_OriSucFail_Energy/ ..._Changed_Energy/ 
#         ..._Energy_with_PCs/ ..._SucFail_Path_Density/

# Load benchmark data.
adjmat=readRDS("./Data_Used/AdjacentMatrix.Matrix.rds");
source("./Auxiliary_Functions_for_Figure/aux_ggplot2Theme.R");
AllResult=readRDS("/dev/shm/All_Result_onlyPou5f1.rds");# This may be slow, please wait.
Rcpp::sourceCpp("./cpp_Source/c_ListMatrixCombined.cpp");
# C++ function for aftercombined_transisMat(...)
CalculateBinnedPC_Coord_Frust<-function(TempEng,TransisMatList,Index){
  # Filter pesudo-energy
  zantai_nengliang=TempEng[Index,];
  zantai_nengliang=as.vector(t(zantai_nengliang));
  index.nl=(zantai_nengliang<1000);cat("Filter.  ");
  zantai_nengliang=zantai_nengliang[index.nl];cat("Over!\n");
  # Filter transient states.
  zantai_zhuangtai=TransisMatList[Index];
  zantai_zhuangtai=aftercombined_transisMat(zantai_zhuangtai);cat("Filter.  ");
  index.zt=(rowSums(zantai_zhuangtai)>0);
  zantai_zhuangtai=zantai_zhuangtai[index.zt,];cat("Over!\n");
  # Check reasonable index, values, & stable states.
  if(sum(index.nl)!=sum(index.zt))stop("Error Data!");
  # Load rotation matrix for embedding phase.
  cat("Calculate PCs.\n");
  rotation=readRDS("./Intermediate_Result/PCA_index_rotation_origin.rds");
  transient_pc_success=zantai_zhuangtai%*%rotation[,1:5];
  return(list(pcs=transient_pc_success, tmp_eng=zantai_nengliang));
}

# Extract appropriate states, paths, & coordinates
num_id_high=c(1:1000000)[(AllResult$New_PS-AllResult$Ori_PS)>6];
num_id_lows=c(1:1000000)[(AllResult$New_PS-AllResult$Ori_PS)<6];
only.high=CalculateBinnedPC_Coord_Frust(
  AllResult$EachTempEng, AllResult$TransisMat, num_id_high);
only.lows=CalculateBinnedPC_Coord_Frust(
  AllResult$EachTempEng, AllResult$TransisMat, num_id_lows);


# Calculate PC code.
source("./Auxiliary_Functions_for_Figure/aux_Binned_PCA_Energy.R");
# Set characher of axises.
pc1_labs=as.character(c(-8.2,-7.1,-6.0,-4.8,-3.6,-2.5));
pc2_labs=as.character(c(-2.5,-1.2,-0.1,1.4,2.7,4.0));
pc3_labs=as.character(c(-2.5,-1.4,-0.3,0.8,1.9,3.0));

# Fig: Mean of upper/lower clusters of temp energy.
TempEnergy=AllResult$EachTempEng;
TempEnergy[TempEnergy==1000]=NA;# <- "1000" is the identifier of invalid records.
TempEnergy_H=TempEnergy[(AllResult$New_PS-AllResult$Ori_PS)>6,];# Temp energy of high cluster
TempEnergy_L=TempEnergy[(AllResult$New_PS-AllResult$Ori_PS)<6,];# Temp energy of low cluster
CountStep_H=rowSums(!is.na(TempEnergy_H), na.rm=TRUE);# Obtain valid reaction step.
CountStep_L=rowSums(!is.na(TempEnergy_H), na.rm=TRUE);# Obtain valid reaction step.
TempEnergy_H=TempEnergy_H[CountStep_H>0,];# more than 1 cases.
TempEnergy_L=TempEnergy_L[CountStep_L<=29,];# more than 1 cases.
Mean_SD_Type<-function(VecHasNA,mingzi,label){
  t=c(1:dim(VecHasNA)[2]);
  x=apply(VecHasNA,2,function(x){mean(x,na.rm = TRUE)});
  x_sd=apply(VecHasNA,2,function(x){sd(x,na.rm = TRUE)});
  return(data.frame(bin=t, av=x, sd=x_sd,pc=label,ty=mingzi));
}
# Dataframe of step with energy.
TempEnergy=rbind.data.frame(
  Mean_SD_Type(TempEnergy_H,"high","Step"),
  Mean_SD_Type(TempEnergy_L,"low", "Step"));
TempEnergy=TempEnergy[!is.na(TempEnergy$sd),];


# Dataframe of each PCx along with energy at high/low/step scale.
df_pc123_energy_alongwith=rbind.data.frame(
  # PC1~PC3 of high 
  data.frame(SingleBinned(only.high$pcs,only.high$tmp_eng,50L,1L), pc="PC1", ty="high"),
  data.frame(SingleBinned(only.high$pcs,only.high$tmp_eng,50L,2L), pc="PC2", ty="high"),
  data.frame(SingleBinned(only.high$pcs,only.high$tmp_eng,50L,3L), pc="PC3", ty="high"),
  # PC1~PC3 of low
  data.frame(SingleBinned(only.lows$pcs,only.lows$tmp_eng,50L,1L), pc="PC1", ty="low"),
  data.frame(SingleBinned(only.lows$pcs,only.lows$tmp_eng,50L,2L), pc="PC2", ty="low"),
  data.frame(SingleBinned(only.lows$pcs,only.lows$tmp_eng,50L,3L), pc="PC3", ty="low")

)

# Figure: change of energy along with PC1, PC2, PC3, 
fig_pc123_energy_alongwith=ggpubr::ggarrange(
  SingleBinned_Plot(df_pc123_energy_alongwith,"PC1",xLab=pc1_labs),
  SingleBinned_Plot(df_pc123_energy_alongwith,"PC2",xLab=pc2_labs),
  SingleBinned_Plot(df_pc123_energy_alongwith,"PC3",xLab=pc3_labs),
  SingleBinned_Plot(TempEnergy,"Step",xLab = seq(0,50,10)),
  ncol=2,nrow=2);
ggsave(filename="./Figure_OriginalVersion/TransientTrajectory_Energy_with_PCs.pdf",
  fig_pc123_energy_alongwith, width=8.0, height=8.0, units="in");

rm(list=ls());gc();
# Code is over.