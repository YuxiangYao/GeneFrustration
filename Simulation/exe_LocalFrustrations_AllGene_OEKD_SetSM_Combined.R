# Combine each individual file as whole file/object.

args <- commandArgs(trailingOnly=TRUE);
param1 <- as.character(args[1]);# Gene Names
part01=readRDS(paste0("/dev/shm/part_",param1,"_n0.rds"));# Note your correct temporary address.
part02=readRDS(paste0("/dev/shm/part_",param1,"_n1.rds"));
part03=readRDS(paste0("/dev/shm/part_",param1,"_n2.rds"));
part04=readRDS(paste0("/dev/shm/part_",param1,"_n3.rds"));
part05=readRDS(paste0("/dev/shm/part_",param1,"_n4.rds"));
part06=readRDS(paste0("/dev/shm/part_",param1,"_n5.rds"));
part07=readRDS(paste0("/dev/shm/part_",param1,"_n6.rds"));
part08=readRDS(paste0("/dev/shm/part_",param1,"_n7.rds"));
part09=readRDS(paste0("/dev/shm/part_",param1,"_n8.rds"));
part10=readRDS(paste0("/dev/shm/part_",param1,"_n9.rds"));
all.data=list(part01,part02,part03,part04,part05,part06,part07,part08,part09,part10);
rm(part01,part02,part03,part04,part05,part06,part07,part08,part09,part10);

mef01=readRDS(paste0("/dev/shm/part_",param1,"_n0_MEFs.rds"));
mef02=readRDS(paste0("/dev/shm/part_",param1,"_n1_MEFs.rds"));
mef03=readRDS(paste0("/dev/shm/part_",param1,"_n2_MEFs.rds"));
mef04=readRDS(paste0("/dev/shm/part_",param1,"_n3_MEFs.rds"));
mef05=readRDS(paste0("/dev/shm/part_",param1,"_n4_MEFs.rds"));
mef06=readRDS(paste0("/dev/shm/part_",param1,"_n5_MEFs.rds"));
mef07=readRDS(paste0("/dev/shm/part_",param1,"_n6_MEFs.rds"));
mef08=readRDS(paste0("/dev/shm/part_",param1,"_n7_MEFs.rds"));
mef09=readRDS(paste0("/dev/shm/part_",param1,"_n8_MEFs.rds"));
mef10=readRDS(paste0("/dev/shm/part_",param1,"_n9_MEFs.rds"));
MEFs=rbind(mef01,mef02,mef03,mef04,mef05,mef06,mef07,mef08,mef09,mef10);
rm(mef01,mef02,mef03,mef04,mef05,mef06,mef07,mef08,mef09,mef10);

# Combined all data.
combined.vector<-function(aList,Whichid){
  tmp=aList[[1]][[Whichid]];
  lens=length(tmp);
  nsam=length(aList);
  returner=rep(NA,lens*nsam);
  left.id=1;right.id=lens;
  for(ii in c(1:length(aList))){
    returner[left.id:right.id]=aList[[ii]][[Whichid]];
    left.id=left.id+lens;
    right.id=right.id+lens;}
  return(returner);
}
combined.matrix<-function(aList,Whichid){
  tmp=aList[[1]][[Whichid]];
  nsam=length(aList);
  nlen=nrow(tmp);
  nval=ncol(tmp);
  returner=matrix(NA,nlen*nsam,nval);
  left.id=1;right.id=nlen;
  for(ii in c(1:length(aList))){
    returner[left.id:right.id,]=aList[[ii]][[Whichid]];
    left.id=left.id+nlen;
    right.id=right.id+nlen;}
  return(returner);
}
combined.list<-function(aList,Whichid){
  tmp=aList[[1]][[Whichid]];
  nsam=length(aList);
  nlen=length(tmp);
  returner=as.list(rep(NA,nlen*nsam));
  left.id=1;right.id=nlen;
  for(ii in c(1:length(aList))){
    returner[left.id:right.id]=aList[[ii]][[Whichid]];
    left.id=left.id+nlen;
    right.id=right.id+nlen;}
  return(returner);
}

AllResult=all.data[[1]]; #as.slot
{
  AllResult$Ori_PS=combined.vector(all.data,1);
  AllResult$Ori_EM=combined.vector(all.data,2);
  AllResult$OriEng=combined.vector(all.data,3);
  AllResult$GeneDiff=combined.vector(all.data,4);
  AllResult$NewEng=combined.vector(all.data,5);
  AllResult$New_PS=combined.vector(all.data,6);
  AllResult$New_EM=combined.vector(all.data,7);
  AllResult$Pressure=combined.vector(all.data,8);
  AllResult$NewTime=combined.vector(all.data,9);
  AllResult$EachUpdate=combined.matrix(all.data,10);
  AllResult$EachTempEng=combined.matrix(all.data,11);
  AllResult$NewState=combined.matrix(all.data,12);
  AllResult$TransisMat=combined.list(all.data,13);
}
rm(combined.vector,combined.matrix,combined.list,all.data);
# The file size is too large. You can choose whether to compress it as required. 
# For compression, the "fastSave" library is recommended.
# saveRDS(AllResult, file="/dev/shm/All_Result_1E6_",param1,".rds", compress=FALSE);

#rm(list=ls());gc();

#To save time, we continually analyze and draw coressponding figure!
# Extract appropriate states, paths, or coordinates
num_id_high=c(1:length(AllResult$New_PS))[(AllResult$New_PS-AllResult$Ori_PS)>6];
if(0==length(num_id_high)||1==length(num_id_high)){
  stop("None of cells fits the condition.\n");
}
# Only save the successful paths.
adjmat=readRDS("../Data_Used/AdjacentMatrix.Matrix.rds");
source("../Auxiliary_Functions_for_Figure/aux_ggplot2Theme.R");

# Part1: Variable genes under pointed gene perturbations.
delta_gene=AllResult$NewState[num_id_high,]-MEFs[num_id_high,];
delta_gene=colSums(abs(delta_gene))/nrow(delta_gene);
# Set spin-like MEF-matrix and adjusted adjacent matrix.
adjmat_0=adjmat;
diag(adjmat_0)=0;# Remove the self-regulatory patterns.
mef_spin=MEFs[num_id_high,];
mef_spin[mef_spin<1]=-1;
# Calculate each gene's local-fru of each cell.
loc.fru.ori=matrix(NA,nrow(mef_spin),ncol(adjmat_0));
colnames(loc.fru.ori)=colnames(adjmat_0);
for(ii in c(1:ncol(mef_spin))){
  loc.fru.ori[,ii]=-(mef_spin%*%adjmat_0[ii,])*mef_spin[,ii];
}
# Calculate each gene's local-fru of each cell. (Have perturbation)
mef_spin[,param1]=1;# Exist external factor!
loc.fru.per=matrix(NA,nrow(mef_spin),ncol(adjmat_0));
colnames(loc.fru.per)=colnames(adjmat_0);
for(ii in c(1:ncol(mef_spin))){
  loc.fru.per[,ii]=-(mef_spin%*%adjmat_0[ii,])*mef_spin[,ii];
}

# Obtain valid varied genes.
loc.fru.del=loc.fru.per-loc.fru.ori;
id_delta_av=colMeans(loc.fru.del[,delta_gene>0]);
id_delta_sd=apply(loc.fru.del[,delta_gene>0],2,sd);

# Prepare dataframe of figure.
non0VarGene=loc.fru.ori[,delta_gene>0];
OriLocalFru=data.frame(
  x=delta_gene[delta_gene>0],# Non-zero genes, changign ratios
  y=colMeans(non0VarGene),# Mean of each gene's local frustration
  er=apply(non0VarGene,2,sd),# Standard Deviation of each gene's local frustration
  t=colMeans(mef_spin[,delta_gene>0]));# Mean of initial values of genes, [-1,+1]
rownames(OriLocalFru)=colnames(mef_spin)[delta_gene>0];

# Show the negative correalation:
tmp.cor=OriLocalFru[abs(id_delta_av)>1e-5&abs(id_delta_sd)>1e-5&rownames(OriLocalFru)!=param1,1:2];
models=lm(y~x, data=tmp.cor);
model.sum=summary(models);
# cor(tmp.cor[,1],tmp.cor[,2])
# Monotone decreasing trend of local frustration vs av.changed ratio.
figs=ggplot(OriLocalFru[abs(id_delta_av)>1e-5&abs(id_delta_sd)>1e-5&rownames(OriLocalFru)!=param1,],
  aes(x=x,y=y))+
  geom_errorbar(aes(x=x,ymin=y-er,ymax=y+er),linewidth=1.0,alpha=0.2)+
  geom_point(aes(fill=t),col="#000000",shape=21,size=2,stroke=0.5)+
  scale_fill_gradient2(low="#8359B5",mid="white",high="#008844",midpoint=0,limits=c(-1,1))+
  labs(x="Av. Changed Ratio",y="Local Frustration");
figs=Plot_ThemeConfiguration(figs)+theme(legend.position = c(0.26,0.24),aspect.ratio=1.0,)+
  annotate("text", x=0.60, y=4, label=paste0("rho=",round(cor(tmp.cor[,1],tmp.cor[,2]),4),
    ",\nR^2=",round(model.sum$adj.r.squared,4)), hjust=0, vjust=0, size=5, colour="#000000")+
  annotate("text", x=0.60, y=-8, label=paste0(param1,"_oe"), 
    hjust=0, vjust=0, size=5, colour="#000000");
ggsave(filename=paste0("/dev/shm/LocalFrustration_vs_AvChangedGene_valid_",param1,".pdf"),
  figs, width=5.5, height=4.0, units="in");

# Monotone decreasing trend of local frustration vs av.changed ratio (all nodes)
figs=ggplot(OriLocalFru,aes(x=x,y=y))+
  geom_errorbar(aes(x=x,ymin=y-er,ymax=y+er),linewidth=1.0,alpha=0.2)+
  geom_point(aes(fill=t),col="#000000",shape=21,size=2,stroke=0.5)+
  scale_fill_gradient2(low="#8359B5",mid="white",high="#008844",midpoint=0)+
  labs(x="Av. Changed Ratio",y="Local Frustration");
figs=Plot_ThemeConfiguration(figs)+theme(legend.position = c(0.26,0.24),aspect.ratio=1.0,)+
  annotate("text", x=0.60, y=-8, label=paste0(param1,"_oe"), 
  hjust=0, vjust=0, size=5, colour="#000000");
ggsave(filename=paste0("/dev/shm/LocalFrustration_vs_AvChangedGene_all_",param1,".pdf"),
  figs, width=5.5, height=4.0, units="in");


# Part2: High/low gene compared as local frustration.
# Load benchmark data.
Rcpp::sourceCpp("../cpp_Source/c_ListMatrixCombined.cpp");
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
  rotation=readRDS("../Intermediate_Result/PCA_index_rotation_origin.rds");
  transient_pc_success=zantai_zhuangtai%*%rotation[,1:5];
  return(list(pcs=transient_pc_success, tmp_eng=zantai_nengliang));
}

# Extract appropriate states, paths, or coordinates
num_id_high=c(1:length(AllResult$New_PS))[(AllResult$New_PS-AllResult$Ori_PS)>6];
# Extract appropriate states, paths, & coordinates
only.high=CalculateBinnedPC_Coord_Frust(AllResult$EachTempEng, AllResult$TransisMat, num_id_high);
# Actually, only.high$tmp_eng should be same as Alls_LocalFru$tmp_eng*88.

# Load the high/low varied genes. (PE-like versus SM-like)
var.gene=readRDS("../Intermediate_Result/var_gene_PE2SM.rds");
gene_hig=var.gene$hvgs;
gene_low=var.gene$lvgs;
# Calculate PC code.
source("../Auxiliary_Functions_for_Figure/aux_Binned_PCA_Energy.R");
# Dataframe of each PCx along with energy at high/low/step scale.
df_pc2_energy_alongwith=data.frame(# PC2 of high 
  SingleBinned(only.high$pcs,only.high$tmp_eng,50L,2L), pc="PC2", ty="high");
# Calculate local frustrations of each cluster 
CalculateBinnedPC_Coord_LocalFrust<-function(Specificgenes,AdjMat,TransisMatList,
  Index, Group_Av=FALSE){
  # Filter transient states.
  zantai_zhuangtai=TransisMatList[Index];
  zantai_zhuangtai=aftercombined_transisMat(zantai_zhuangtai);cat("Filter.  ");
  colnames(zantai_zhuangtai)=colnames(AdjMat);
  index.zt=(rowSums(zantai_zhuangtai)>0);
  zantai_zhuangtai=zantai_zhuangtai[index.zt,];cat("Over!\n");
  #Convertas spin-like:
  sp_mat=zantai_zhuangtai;
  sp_mat[sp_mat<1]=-1;
  # Set return matrix.
  ResLocFru=matrix(NA,nrow(sp_mat),length(Specificgenes));
  colnames(ResLocFru)=Specificgenes;
#   for(ii in c(1:length(Specificgenes))){
#     ResLocFru[,ii]=-sp_mat%*%AdjMat[Specificgenes[ii],]*sp_mat[,Specificgenes[ii]];
#   }
  ResLocFru=-sp_mat%*%t(AdjMat[Specificgenes,])*sp_mat[,Specificgenes];
  # Load rotation matrix for embedding phase.
  cat("Calculate PCs.\n");
  rotation=readRDS("../Intermediate_Result/PCA_index_rotation_origin.rds");
  transient_pc_success=zantai_zhuangtai%*%rotation[,1:5];
  if(Group_Av){
    ResLocFru=apply(ResLocFru,1,mean);}
  return(list(pcs=transient_pc_success, mat=zantai_zhuangtai, tmp_eng=ResLocFru));
}

# Only HVGs.
HVGs_LocalFru=CalculateBinnedPC_Coord_LocalFrust(gene_hig, adjmat, 
  AllResult$TransisMat, num_id_high, TRUE);
# Only LVGs.
LVGs_LocalFru=CalculateBinnedPC_Coord_LocalFrust(gene_low, adjmat, 
  AllResult$TransisMat, num_id_high, TRUE);
# All genes as benchmark.
Alls_LocalFru=CalculateBinnedPC_Coord_LocalFrust(colnames(adjmat), adjmat, 
  AllResult$TransisMat, num_id_high, TRUE);

# Dataframe of each PCx along with energy at high/low/step scale.
LocalFru_AlongWithPC2=rbind.data.frame(
  data.frame(SingleBinned(Alls_LocalFru$pcs, Alls_LocalFru$tmp_eng, 50L,2L), pc="PC2", ty="All"),
  data.frame(SingleBinned(LVGs_LocalFru$pcs, LVGs_LocalFru$tmp_eng, 50L,2L), pc="PC2", ty="LVGs"),
  data.frame(SingleBinned(HVGs_LocalFru$pcs, HVGs_LocalFru$tmp_eng, 50L,2L), pc="PC2", ty="HVGs"));
LocalFru_AlongWithPC2$ty=factor(LocalFru_AlongWithPC2$ty, levels=c("HVGs","LVGs","All"));

# Figure: change of energy along with PC1, PC2, PC3, 
fig_=
  SingleBinned_Plot(LocalFru_AlongWithPC2,"PC2",xLab=as.character(c(-2.5,-1.2,-0.1,1.4,2.7,4.0)))+
  guides(color = guide_legend(title = paste0(param1,"_kd")),  
    linetype= guide_legend(title = paste0(param1,"_kd")),
    fill = guide_legend(title = paste0(param1,"_kd"))) +  
  labs(y="Local Frustration")+theme(legend.position=c(0.55,0.25),);
ggsave(filename=paste0("/dev/shm/LocalFrustration_High2Low_",param1,"_inProcess.pdf"),
  fig_, width=3.5, height=3.5, units="in");

# Each single item compared between original SSS and induced processes.
EachSingleBinned_Plot<-function(DataFrame,Color,ID,xLim,xLab){
  #(TempEngs,ID,xLim,xLab,DashedPoint)
  fig_tempeng=ggplot(DataFrame[DataFrame$pc==ID,],aes(x=bin,y=av))+
    geom_ribbon(aes(ymin=av-sd, ymax=av+sd,fill=cell),alpha=0.13) +
    geom_line(aes(linetype=cell,col=cell),linewidth=1.2)+
    #scale_color_manual(values=c("#2156A6","#FFB600","#EE312E","#6c6c6c","#00ff00"))+
    #scale_fill_manual(values=c("#2156A6","#FFB600","#EE312E","#6c6c6c","#00ff00"))+
    scale_fill_manual(values=c("#2156A6","#EE312E"))+
    scale_color_manual(values=c("#2156A6","#EE312E"))+
    scale_linetype_manual(values=c("solid","dashed"))+
    labs(x=ID,y="pseudoenergy")+
    scale_x_continuous(breaks=c(0,10,20,30,40,50), labels=xLab, expand=c(0,0),);
  fig_tempeng=Plot_ThemeConfiguration(fig_tempeng)+theme(legend.position="right",aspect.ratio=1.0);
  return(fig_tempeng);
}
# Load the having been calculated Original landscape.
load("../Intermediate_Result/OriginalPhoneLandLocalFru.rds");

LocalFru_AlongWithPC2_Per_Ori=rbind.data.frame(
  data.frame(SingleBinned(Alls_LocalFru$pcs, Alls_LocalFru$tmp_eng, 50L,2L), pc="PC2", ty="All",cell="Per"),
  data.frame(SingleBinned(LVGs_LocalFru$pcs, LVGs_LocalFru$tmp_eng, 50L,2L), pc="PC2", ty="LVGs",cell="Per"),
  data.frame(SingleBinned(HVGs_LocalFru$pcs, HVGs_LocalFru$tmp_eng, 50L,2L), pc="PC2", ty="HVGs",cell="Per"),
  data.frame(SingleBinned(Alls_LocalFru_Ori$pcs, Alls_LocalFru_Ori$tmp_eng, 50L,2L), pc="PC2", ty="All",cell="Ori"),
  data.frame(SingleBinned(LVGs_LocalFru_Ori$pcs, LVGs_LocalFru_Ori$tmp_eng, 50L,2L), pc="PC2", ty="LVGs",cell="Ori"),
  data.frame(SingleBinned(HVGs_LocalFru_Ori$pcs, HVGs_LocalFru_Ori$tmp_eng, 50L,2L), pc="PC2", ty="HVGs",cell="Ori") );
LocalFru_AlongWithPC2_Per_Ori$ty=factor(LocalFru_AlongWithPC2_Per_Ori$ty, 
    levels = c("HVGs","MVGs","LVGs","All"));
LocalFru_AlongWithPC2_Per_Ori$cell=factor(LocalFru_AlongWithPC2_Per_Ori$cell, 
    levels = c("Per","Ori"));

# Each part (High/low/All compared with original landscapes)
fig_3=ggpubr::ggarrange(
  EachSingleBinned_Plot(LocalFru_AlongWithPC2_Per_Ori[LocalFru_AlongWithPC2_Per_Ori$ty=="HVGs",],
    "#2156A6","PC2",xLab=as.character(c(-2.5,-1.2,-0.1,1.4,2.7,4.0)))+
    guides(color = guide_legend(title = paste0(param1,"_")),  
    linetype= guide_legend(title = paste0(param1,"_")),
    fill = guide_legend(title = paste0(param1,"_"))) + 
    labs(y="Local Frustration")+theme(legend.position=c(0.55,0.25),),
  EachSingleBinned_Plot(LocalFru_AlongWithPC2_Per_Ori[LocalFru_AlongWithPC2_Per_Ori$ty=="LVGs",],
    "#EE312E","PC2",xLab=as.character(c(-2.5,-1.2,-0.1,1.4,2.7,4.0)))+
    guides(color = guide_legend(title = paste0(param1,"_")),  
    linetype= guide_legend(title = paste0(param1,"_kd")),
    fill = guide_legend(title = paste0(param1,"_kd"))) + 
    labs(y="Local Frustration")+theme(legend.position=c(0.55,0.25),),
  EachSingleBinned_Plot(LocalFru_AlongWithPC2_Per_Ori[LocalFru_AlongWithPC2_Per_Ori$ty=="All",],
    "#6c6c6c","PC2",xLab=as.character(c(-2.5,-1.2,-0.1,1.4,2.7,4.0)))+
    guides(color = guide_legend(title = paste0(param1,"_")),  
    linetype= guide_legend(title = paste0(param1,"_")),
    fill = guide_legend(title = paste0(param1,"_"))) + 
    labs(y="Local Frustration")+theme(legend.position=c(0.55,0.25),),
  ncol=3,nrow=1);

ggsave(filename=paste0("/dev/shm/LocalFrustration_High2Low_",param1,"_Per2Ori.pdf"),
  fig_3, width=10.5, height=3.5, units="in");


# saveRDS(only.high,"/dev/shm/Only_High.rds");# You can check them if you want.
{AllResult$Ori_PS=AllResult$Ori_PS[num_id_high];
AllResult$Ori_EM=AllResult$Ori_EM[num_id_high];
AllResult$OriEng=AllResult$OriEng[num_id_high];
AllResult$GeneDiff=AllResult$GeneDiff[num_id_high];
AllResult$NewEng=AllResult$NewEng[num_id_high];
AllResult$New_PS=AllResult$New_PS[num_id_high];
AllResult$New_EM=AllResult$New_EM[num_id_high];
AllResult$Pressure=AllResult$Pressure[num_id_high];
AllResult$NewTime=AllResult$NewTime[num_id_high];
AllResult$EachUpdate=AllResult$EachUpdate[num_id_high,];
AllResult$EachTempEng=AllResult$EachTempEng[num_id_high,];
AllResult$NewState=AllResult$NewState[num_id_high,];
AllResult$TransisMat=AllResult$TransisMat[num_id_high];}

saveRDS(AllResult,paste0("./Intermediate_Result/OnlySuccessCase_",param1,".rds"),compress=FALSE);
saveRDS(MEFs[num_id_high,],paste0("./Intermediate_Result/OnlySuccessCase_",param1,"_MEFs.rds"),compress=FALSE,)
# Code is over.
