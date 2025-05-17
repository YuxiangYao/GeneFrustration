# Compared with simulated results of Pou5f1-OE:
# Load benchmark data.
adjmat=readRDS("./Data_Used/AdjacentMatrix.Matrix.rds");
source("./Auxiliary_Functions_for_Figure/aux_ggplot2Theme.R");
MEFs=readRDS("./Intermediate_Result/ExtremeCase_FourClusters.rds");
MEFs=MEFs$sm;
AllResult=readRDS("/dev/shm/All_Result_onlyPou5f1.rds");# This may be slow, please wait.

# Obtain two clusters.
Index_High=c(1:1000000)[AllResult$New_PS-AllResult$Ori_PS>6];# The successful paths.

# Extract successful paths.
ExtractTransientCell<-function(aMatrix){# center=TRUE
  # Row: (sample) each time point, Col: gene.
  tmp.mat=aMatrix[rowSums(aMatrix)>0,];
  return(tmp.mat);
}
# Calucate the driving force during reprogramming
DrivingForce<-function(aVec,AdjMat,Removed){# Also known as induced frustration.
  CanBeChanged=aVec%*%t(AdjMat);
  CanBeChanged=sign(CanBeChanged);
  ind1=(aVec[-Removed]==1)&(CanBeChanged[-Removed]<0);
  ind2=(aVec[-Removed]==0)&(CanBeChanged[-Removed]>0);
  # CanBeChanged==0, keep same values of original aVec.
  DF=sum(ind1|ind2);
  return(DF);
}

# Output:
set.seed(20250101L);
Filtered_Succ=lapply(AllResult$TransisMat[Index_High], ExtractTransientCell);
Final_Succ=do.call(rbind, Filtered_Succ);
colnames(Final_Succ)=colnames(adjmat);
DriveForce=apply(Final_Succ,1,DrivingForce,adjmat,which(colnames(adjmat)=="Pou5f1"));
rotation=readRDS("./Intermediate_Result/PCA_index_rotation_origin.rds");
transient_pc_success=Final_Succ%*%rotation[,1:5];


SingleBinned<-function(PCAss,Values,BinWidth=40L,oneID){
  zuobiao=(PCAss[,oneID]);
  if(1==oneID){# PC1
    x_bian=c(-8.2,-2.5);
  } else if(2==oneID){# PC2
    x_bian=c(-2.5,4.0);
  } else {# PC3
    x_bian=c(-2.5,3.0);}
  seq_x=seq(x_bian[1],x_bian[2],length.out=(BinWidth+1));
  bin_x=floor((zuobiao-x_bian[1])/(seq_x[2]-seq_x[1]));
  dfs=cbind.data.frame(num=bin_x,val=Values);
  # Binning ....
  tmps=tapply(dfs$val,dfs$num,function(x){c(mean(x),sd(x))})
  tmps=matrix(unlist(tmps),ncol=2,byrow=T);
  res=matrix(NA,BinWidth,3);
  res[,1]=c(1:BinWidth);
  rownames(res)=paste0("B_",c(1:BinWidth));
  colnames(res)=c("bin","av","sd");
  res[c(min(bin_x):max(bin_x))+1,2:3]=tmps;
  return(res);
}

SingleBinned_Plot<-function(TempEngs,ID,xLim,xLab,DashedPoint){
  fig_tempeng=ggplot(TempEngs[TempEngs$pc==ID,],aes(x=bin,y=av))+
    geom_ribbon(aes(ymin=av-sd, ymax=av+sd,fill=ty),alpha=0.15) +
    geom_line(aes(color=ty,linetype=ty),linewidth=1.2)+
    # scale_color_manual(values=c("#2156A6","#FFB600","#EE312E","#6c6c6c","#00ff00"))+
    # scale_fill_manual(values=c("#2156A6","#FFB600","#EE312E","#6c6c6c","#00ff00"))+
    scale_color_manual(values=c("#000000","#2156A6","#EE312E","#FFB600","#00ff00"))+
    scale_fill_manual(values=c("#6c6c6caa","#2156A6","#EE312E","#FFB600","#00ff00"))+
    scale_linetype_manual(values = c("solid","solid","dashed"))+
    labs(x=ID,y="pseudoenergy")+
    scale_x_continuous(breaks=c(0,10,20,30,40,50), labels=xLab, expand=c(0,0),);
  fig_tempeng=Plot_ThemeConfiguration(fig_tempeng)+theme(legend.position="right",aspect.ratio=1.0);
  return(fig_tempeng);
}

# Dataframe of each PCx along with energy at high/low/step scale.
tmp_bin_df=data.frame(SingleBinned(transient_pc_success, DriveForce, 50L,2L), pc="PC2", ty="All");

# Figure: change of energy along with PC1, PC2, PC3, 
fig_=
  SingleBinned_Plot(tmp_bin_df,"PC2",xLab=as.character(c(-2.5,-1.2,-0.1,1.4,2.7,4.0)))+
  labs(y="Local Frustration")+theme(legend.position=c(0.55,0.25),);
ggsave(filename="./Figure_OriginalVersion/DrivingForce.pdf",
  fig_, width=3.5, height=3.5, units="in");





# Code is over.