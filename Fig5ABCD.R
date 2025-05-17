# Real experimental data:
# RNA-sequencing for validating accumulations in transitions.
# Real experimental data.
source("./Auxiliary_Functions_for_Figure/aux_ggplot2Theme.R");
source("./Auxiliary_Functions_for_Figure/aux_PhenotypicScore.R");
Rcpp::sourceCpp("./cpp_Source/c_SmoothKNN.cpp");
source("./Auxiliary_Functions_for_Figure/aux_Mixed2ColorBar.R");
source("./Auxiliary_Functions_for_Figure/aux_FilterValidBinarizeGene.R");
# Set intersection sub-PCA picture:
source("./Auxiliary_Functions_for_Figure/aux_StableStatePCA.R");
library(ggplot2,ggnewscale);# for multiple scale_fill_bars
# Load basic data.
ori.ss=readRDS("./Intermediate_Result/SSS_origin_nodup.rds");
adjmat=readRDS("./Data_Used/AdjacentMatrix.Matrix.rds");
# Function: Prepare the dataframe for RNA-Phenotypic-PCA.
figa_expdata_pre<-function(ori.data,exp.data){
  sss=ori.data;
  cat("Remove duplicated states ...\n");
  sss=sss[,intersect(colnames(sss),colnames(exp.data))];
  sss=sss[!duplicated(sss),];
  cat("PCA ...\n");
  ind=PCA_index_StableState(sss);# Embed SS into a PCA-space.
  pca_index=cbind.data.frame(ind$dims);# PCA coordinates.
  pca_index=as.matrix(pca_index[,2:3]);
  cat("KNN ...\n");
  all_knn=RcppHNSW::hnsw_knn(pca_index,k=5,distance="l2",verbose=TRUE,n_threads=64,progress="bar");
  cat("Set figure dataframe ...\n");
  score_1=SmoothKNN(all_knn[[1]]-1,SPT_score(sss),5);# PS
  score_2=SmoothKNN(all_knn[[1]]-1,EMT_score(sss),5);# EM
  s1=(score_1-min(score_1))/(max(score_1)-min(score_1));# Linear scaling
  s2=(score_2-min(score_2))/(max(score_2)-min(score_2));# Linear scaling
  ColorBlocks=ColMix_Vec(c('#999999',"#9999ff","#ff9999","#99ff99"),cbind(s1,s2));
  pca_bg=cbind.data.frame(x=pca_index[,1],y=pca_index[,2],c=ColorBlocks,s1=score_1,s2=score_2);
  colortype=unique(pca_bg$c);
  pca_bg$c=factor(pca_bg$c,levels=colortype);
  return(list(dt=pca_bg,cb=colortype,rota=ind$rotation));
}
# Remove texts of axis+title.
Fig_RemoveText<-function(ggplotObj){
  fig=Plot_ThemeConfiguration(ggplotObj)+
    theme(legend.position="none",aspect.ratio=1.0,
          axis.title.x.bottom=element_blank(), axis.line.x.bottom=element_blank(),
          axis.line.y.left=element_blank(),axis.title.y.left=element_blank(),
          axis.text.x.bottom = element_blank(),axis.text.y.left = element_blank(),);
  return(fig);
}
# Figure: Output figures background + scatters.
Fig_Expdata_Embed<-function(exp.dt,PreList,FileName){
  pca_bg=PreList$dt;  colortype=PreList$cb;
  # Figure: colorful background.
  fig_bg=ggplot()+geom_point(data=pca_bg,aes(x=x,y=y,col=c),size=0.5,shape=16,alpha=0.30)+
    scale_color_manual(values=colortype)+
    scale_x_continuous(breaks = c(-2,-1,0,1,2,3),limits = c(-2.7,3.8),expand = c(0,0))+
    scale_y_continuous(breaks = c(-2,-1,0,1,2),limits = c(-2.5,2.8),expand = c(0,0));
  fig_bg=Fig_RemoveText(fig_bg);
  # Figure: scater points of each time.
  fig_point=ggplot()+geom_point(data=exp.dt,aes(x=x,y=y,shape=f,col=f,size=f),stroke=1.5)+
    scale_x_continuous(breaks = c(-2,-1,0,1,2,3),limits = c(-2.7,3.8),expand = c(0,0))+
    scale_y_continuous(breaks = c(-2,-1,0,1,2),limits = c(-2.5,2.8),expand = c(0,0))+
    scale_shape_manual(values=exp.dt$s)+scale_size_manual(values=exp.dt$l)+
    scale_color_manual(values=exp.dt$c);
  fig_point=Fig_RemoveText(fig_point);
  # Output the figures.
  ggsave(paste0("./Figure_OriginalVersion/",FileName,"_bg.png"),fig_bg,height=3.5,width=3.5,units="in",dpi=300);
  ggsave(paste0("./Figure_OriginalVersion/",FileName,"_points.pdf"),fig_point,height=3.5,width=3.5,units="in",dpi=300);
  return(fig_point);
}
# Function: return instant candidate number.
InstantTr<-function(Expmat,AdjMat_Ori,Vals=1L){#Row: samples; Col: genes
  index=intersect(colnames(AdjMat_Ori),colnames(Expmat));
  adjm=AdjMat_Ori[index,index];
  expmat=Expmat[,index];
  inducer=rep(0,ncol(expmat));
  names(inducer)=colnames(expmat);
  inducer[c("Pou5f1","Sox2","Klf4","Myc")]=Vals;# Yamanaka protocol=1, non-Y protocol=0
  trigger=rep(NA,nrow(expmat));
  for(ii in c(1:nrow(expmat))){
    realss=expmat[ii,];
    pesudo=as.integer(realss|inducer);
    realssx=adjm%*%pesudo;# Candidate genes.
    trigger[ii]=sum((realssx>0&realss==0)|(realssx<0&realss==1));# Obtain number of candidates
  }
  res=cbind.data.frame(trigger=trigger);
  rownames(res)=rownames(Expmat);
  return(res);
}
# Function: Output figures for
delta_Tr_Histogram<-function(DataFrame,Filename){
  pica=ggplot(DataFrame, aes(x=x,fill=y))+
    geom_histogram(position = "identity", alpha = 0.55, binwidth=1) + 
    scale_fill_manual(values=c("#2156A6","#EE312E","#FFB600"))#+
  pica=Plot_ThemeConfiguration(pica);
  ggsave(paste0("./Figure_OriginalVersion/",Filename,".pdf"),pica,height=1.7,width=3.5,units="in",dpi=300);
  return(pica);
}
# Function: Output figures for real phenotypic scores.
Real_Pheno_Score<-function(DataFrame,Filename,ScoreType='ps'){
  if('ps'==ScoreType){
    figx=ggplot(DataFrame)+# P/S score.
      geom_line(aes(x=tt,y=spt,col=type),linewidth=1)+
      geom_point(aes(x=tt,y=spt,col=type,shape=type),size=4.5);
  } else {
    figx=ggplot(DataFrame)+# P/S score.
      geom_line(aes(x=tt,y=met,col=type),linewidth=1)+
      geom_point(aes(x=tt,y=met,col=type,shape=type),size=4.5);}
  figx=figx+scale_color_manual(values=c("#2156A6","#EE312E","#FFB600"))+
    scale_shape_manual(values=c(16,17,15))+
    scale_x_continuous(expand = c(0,0),limits = c(-1.5,9.5),
      breaks=c(-1,0,1,3,5,7,8,9),labels=c("MEF",paste0("Day",c(0,1,3,5,7)),"IPS","ESC"));
  # Output figure.
  figx=Plot_ThemeConfiguration(figx)+theme(legend.position = c(0.50,0.25),
    axis.text.x=element_text(angle=45,vjust=1.0,hjust=1.0));
  return(figx);
}


# Li2017 (GSE93027) ======================================
set.seed(93027);
g93027=read.table("./Data_Used/GSE93027.tsv",sep="\t");
g93027=g93027[,1:24];
sam.name=colnames(g93027);
g93027=g93027_ori=as.matrix(g93027[Keep.Gene(rownames(g93027),GeneList),]);
g93027=(as.matrix(apply(log1p(g93027),1,Binarization_GMM)));
rownames(g93027)=sam.name;
pre.list=figa_expdata_pre(ori.ss,g93027);# ~5min
expdata=g93027[,rownames(pre.list$rota)]%*%pre.list$rota[,2:3];
expdata=cbind.data.frame(x=expdata[,1],y=expdata[,2],
                         t=rep(c("ESCs","MEF","D0","D1","D3","D5","D7","D0-Jun","D1-Jun","D5-Jun","D7-Jun","iPS"),each=2),
                         s=rep(c("\u25CE","\u25C6",0,1,3,5,7, 0,1,3,7,"\u25C9"),each=2),
                         l=rep(c(rep(4.5,2),3.5,rep(4,4),3.5,rep(4,4)),each=2),
                         c=rep(c("black","black",rep("#2156A6",5),rep("#EE312E",4),"black"),each=2),
                         f=sprintf("%02d",c(1:nrow(expdata))) );

# Fig. 5(b): PCA of Li2017 data.
fig_93027=Fig_Expdata_Embed(expdata,pre.list,"Li2017");# ~3min


# Obtain real phenotypic scores.
g93027_p=t(g93027_ori) # Has been a realdata matrix
g93027_p=g93027_p[,intersect(colnames(adjmat),colnames(g93027_p))];
g93027_p=log(g93027_p+1,2);
score1=SPT_score(g93027_p);
score1=0.5*(score1[seq(1,23,2)]+score1[seq(2,24,2)])
score2=EMT_score(g93027_p);
score2=0.5*(score2[seq(1,23,2)]+score2[seq(2,24,2)])
rPhenoScore=data.frame(tt=c(9,-1,0,1,3,5,7,0,1,3,7,8),# Time points
  spt=score1,met=score2,type=c(rep("A",7),rep("B",4),"A"));
rPhenoScore=rbind.data.frame(rPhenoScore,rPhenoScore[2,]);
rPhenoScore$type[13]="B";# Last node is used for linking adjacent nodes.
# Fig.5(c): Curve and points figure.
li2017_rPheno1=Real_Pheno_Score(rPhenoScore,"Li2017_SPT",'ps');
li2017_rPheno2=Real_Pheno_Score(rPhenoScore,"Li2017_MET",'em');

ggsave(filename="./Figure_OriginalVersion/RealData_G093027_Pheno1.pdf",
  li2017_rPheno1, height=2.0,width=5.0,units="in",dpi=300);
ggsave(filename="./Figure_OriginalVersion/RealData_G093027_Pheno2.pdf",
  li2017_rPheno2, height=2.0,width=5.0,units="in",dpi=300);



# Energy >>>
g93027_spin=g93027;g93027_spin[g93027_spin<1]=-1;
adjmat.93027=adjmat[colnames(g93027_spin),colnames(g93027_spin)];

energys=rep(NA, nrow(g93027_spin));
for(ii in c(1:nrow(g93027_spin))){
  energys[ii]=-sum((g93027_spin[ii,]%*%adjmat.93027)*g93027_spin[ii,]);
}
names(energys)=rownames(g93027_spin)


g93027_df=cbind.data.frame(days=c(9,-1,0,1,3,5,7,0,1,3,7,8),
    type=letters[c(3,3,1,1,1,1,1,2,2,2,2,3)],
    eng_av=apply(matrix(energys,ncol=2,byrow = TRUE),1,mean),
    eng_sd=apply(matrix(energys,ncol=2,byrow = TRUE),1,sd)
  )


# Monotone decreasing trend of local frustration vs av.changed ratio (all nodes)
figs=ggplot()+
  geom_line(data=g93027_df[c(2:7,12),],aes(x=days,y=eng_av),col="#2156A6",linewidth=0.8)+
  geom_line(data=g93027_df[c(2,8:11),],aes(x=days,y=eng_av),col="#EE312E",linewidth=0.8,linetype="dashed")+
  geom_errorbar(data=g93027_df,aes(x=days,ymin=eng_av-eng_sd,ymax=eng_av+eng_sd,col=type),
    linewidth=1.0,alpha=0.3)+
  geom_point(data=g93027_df,aes(x=days,y=eng_av,col=type),shape=16,size=2)+
  scale_color_manual(values = c("#2156A6","#EE312E","#000000"),
    labels = c("a" = "Normal", "b" = "JunOE", "c" = "Other") )+
  labs(x="Time Point",y="Pseudo Energy")+
  scale_x_continuous(expand = c(0,0),limits = c(-1.5,9.5),
    breaks=c(-1,0,1,3,5,7,8,9),labels=c("MEF",paste0("Day",c(0,1,3,5,7)),"IPS","ESC"))+
  guides(color=guide_legend(nrow=1,ncol=3,byrow=TRUE));
figs=Plot_ThemeConfiguration(figs)+theme(legend.position = c(0.50,0.25),
  axis.text.x=element_text(angle=45,vjust=1.0,hjust=1.0));
ggsave(filename="./Figure_OriginalVersion/RealData_G093027_TimeEnergy.pdf",
  figs, height=2.0,width=5.0,units="in",dpi=300);


# Real experimental data.
source("./Auxiliary_Functions_for_Figure/aux_ggplot2Theme.R");
source("./Auxiliary_Functions_for_Figure/aux_PhenotypicScore.R");
Rcpp::sourceCpp("./cpp_Source/c_SmoothKNN.cpp");
source("./Auxiliary_Functions_for_Figure/aux_Mixed2ColorBar.R");
source("./Auxiliary_Functions_for_Figure/aux_FilterValidBinarizeGene.R");
# Set intersection sub-PCA picture:
source("./Auxiliary_Functions_for_Figure/aux_StableStatePCA.R");
set.seed(93027);
g93027=read.table("/mnt/yaoyuxiang/20230616/ExpData/GSE93027.tsv",sep="\t");
g93027=g93027[,1:24];
sam.name=colnames(g93027);
g93027=g93027_ori=as.matrix(g93027[Keep.Gene(rownames(g93027),GeneList),]);
g93027=(as.matrix(apply(log1p(g93027),1,Binarization_GMM)));
rownames(g93027)=sam.name;
g93027_spin=g93027;g93027_spin[g93027_spin<1]=-1;
adjmat.93027=adjmat[colnames(g93027_spin),colnames(g93027_spin)];

energys=rep(NA, nrow(g93027_spin));
for(ii in c(1:nrow(g93027_spin))){
  energys[ii]=-sum((g93027_spin[ii,]%*%t(adjmat.93027))*g93027_spin[ii,]);
}
names(energys)=rownames(g93027_spin);



# Added 1209
energys=-g93027_spin%*%t(adjmat.93027)*g93027_spin;
rownames(energys)=rownames(g93027_spin);
colnames(energys)=colnames(adjmat.93027);

#threetype=readRDS("./Intermediate_Result/ThreeTypeGene.rds")
twotypes=readRDS("./Intermediate_Result/var_gene_PE2SM.rds");
highs=rowSums(energys[,intersect(colnames(adjmat.93027),twotypes$hvgs)]);
lowss=rowSums(energys[,intersect(colnames(adjmat.93027),twotypes$lvgs)]);

g93027_df_high=cbind.data.frame(days=c(9,-1,0,1,3,5,7,0,1,3,7,8),
    type=letters[c(3,3,1,1,1,1,1,2,2,2,2,3)],
    eng_av=apply(matrix(highs,ncol=2,byrow = TRUE),1,mean),
    eng_sd=apply(matrix(highs,ncol=2,byrow = TRUE),1,sd))

g93027_df_lows=cbind.data.frame(days=c(9,-1,0,1,3,5,7,0,1,3,7,8),
    type=letters[c(3,3,1,1,1,1,1,2,2,2,2,3)],
    eng_av=apply(matrix(lowss,ncol=2,byrow = TRUE),1,mean),
    eng_sd=apply(matrix(lowss,ncol=2,byrow = TRUE),1,sd))

# Monotone decreasing trend of local frustration vs av.changed ratio (all nodes summation)
figs_LF_h=ggplot()+
  geom_line(data=g93027_df_high[c(2:7,12),],aes(x=days,y=eng_av),col="#2156A6",linewidth=0.8)+
  geom_line(data=g93027_df_high[c(2,8:11),],aes(x=days,y=eng_av),col="#EE312E",linewidth=0.8,linetype="dashed")+
  geom_errorbar(data=g93027_df_high,aes(x=days,ymin=eng_av-eng_sd,ymax=eng_av+eng_sd,col=type),
    linewidth=1.0,alpha=0.3)+
  geom_point(data=g93027_df_high,aes(x=days,y=eng_av,col=type),shape=16,size=2)+
  scale_color_manual(values = c("#2156A6","#EE312E","#000000"),
    labels = c("a" = "Normal", "b" = "JunOE", "c" = "Other") )+
  labs(x="Time Point",y="Pseudo Energy")+
  scale_x_continuous(expand = c(0,0),limits = c(-1.5,9.5),
    breaks=c(-1,0,1,3,5,7,8,9),labels=c("MEF",paste0("Day",c(0,1,3,5,7)),"IPS","ESC"))+
  guides(color=guide_legend(nrow=1,ncol=3,byrow=TRUE));
figs_LF_h=Plot_ThemeConfiguration(figs_LF_h)+theme(legend.position = c(0.50,0.25),
  axis.text.x=element_text(angle=45,vjust=1.0,hjust=1.0));
ggsave(filename="./Figure_OriginalVersion/RealData_G093027_TimeLocalFru_H_sum.pdf",
  figs_LF_h, height=2.0,width=5.0,units="in",dpi=300);

figs_LF_l=ggplot()+
  geom_line(data=g93027_df_lows[c(2:7,12),],aes(x=days,y=eng_av),col="#2156A6",linewidth=0.8)+
  geom_line(data=g93027_df_lows[c(2,8:11),],aes(x=days,y=eng_av),col="#EE312E",linewidth=0.8,linetype="dashed")+
  geom_errorbar(data=g93027_df_lows,aes(x=days,ymin=eng_av-eng_sd,ymax=eng_av+eng_sd,col=type),
    linewidth=1.0,alpha=0.3)+
  geom_point(data=g93027_df_lows,aes(x=days,y=eng_av,col=type),shape=16,size=2)+
  scale_color_manual(values = c("#2156A6","#EE312E","#000000"),
    labels = c("a" = "Normal", "b" = "JunOE", "c" = "Other") )+
  labs(x="Time Point",y="Pseudo Energy")+
  scale_x_continuous(expand = c(0,0),limits = c(-1.5,9.5),
    breaks=c(-1,0,1,3,5,7,8,9),labels=c("MEF",paste0("Day",c(0,1,3,5,7)),"IPS","ESC"))+
  guides(color=guide_legend(nrow=1,ncol=3,byrow=TRUE));
figs_LF_l=Plot_ThemeConfiguration(figs_LF_l)+theme(legend.position = c(0.50,0.25),
  axis.text.x=element_text(angle=45,vjust=1.0,hjust=1.0));
ggsave(filename="./Figure_OriginalVersion/RealData_G093027_TimeLocalFru_L_sum.pdf",
  figs_LF_l, height=2.0,width=5.0,units="in",dpi=300);


# The high-succ genes are summarized as one values, should be independent.

twotypes=readRDS("./Intermediate_Result/var_gene_PE2SM.rds");
highs=(energys[,intersect(colnames(adjmat.93027),twotypes$hvgs)])
lowss=(energys[,intersect(colnames(adjmat.93027),twotypes$lvgs)])
highs=cbind(highs[2*c(0:11)+1,],highs[2*c(1:12),]);
lowss=cbind(lowss[2*c(0:11)+1,],lowss[2*c(1:12),]);
g93027_df_high=cbind.data.frame(days=c(9,-1,0,1,3,5,7,0,1,3,7,8),
  type=letters[c(3,3,1,1,1,1,1,2,2,2,2,3)],
  eng_av=apply(highs, 1, mean),
  eng_sd=apply(highs, 1, sd))

g93027_df_lows=cbind.data.frame(days=c(9,-1,0,1,3,5,7,0,1,3,7,8),
  type=letters[c(3,3,1,1,1,1,1,2,2,2,2,3)],
  eng_av=apply(lowss, 1, mean),
  eng_sd=apply(lowss, 1, sd))

# Monotone decreasing trend of local frustration vs av.changed ratio (all nodes average)
figs_LF_h=ggplot()+
  geom_line(data=g93027_df_high[c(2:7,12),],aes(x=days,y=eng_av),col="#2156A6",linewidth=0.8)+
  geom_line(data=g93027_df_high[c(2,8:11),],aes(x=days,y=eng_av),col="#EE312E",linewidth=0.8,linetype="dashed")+
  geom_errorbar(data=g93027_df_high,aes(x=days,ymin=eng_av-eng_sd,ymax=eng_av+eng_sd,col=type),
                linewidth=1.0,alpha=0.3)+
  geom_point(data=g93027_df_high,aes(x=days,y=eng_av,col=type),shape=16,size=2)+
  scale_color_manual(values = c("#2156A6","#EE312E","#000000"),
                     labels = c("a" = "Normal", "b" = "JunOE", "c" = "Other") )+
  labs(x="Time Point",y="Pseudo Energy")+
  scale_x_continuous(expand = c(0,0),limits = c(-1.5,9.5),
                     breaks=c(-1,0,1,3,5,7,8,9),labels=c("MEF",paste0("Day",c(0,1,3,5,7)),"IPS","ESC"))+
  guides(color=guide_legend(nrow=1,ncol=3,byrow=TRUE));
figs_LF_h=Plot_ThemeConfiguration(figs_LF_h)+theme(legend.position = c(0.50,0.25),
                                                   axis.text.x=element_text(angle=45,vjust=1.0,hjust=1.0));

ggsave(filename="./Figure_OriginalVersion/RealData_G093027_TimeLocalFru_H_av.pdf",
       figs_LF_h, height=2.0,width=5.0,units="in",dpi=300);
