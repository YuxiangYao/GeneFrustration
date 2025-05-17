# Figure: OE/KD unique/combined factors to perturb SM-like cells.

source("./Simulation/exe_GenerateFourCluster.R");# Generate data.
# Load data, function.
Rcpp::sourceCpp("./cpp_Source/c_OEKD_StateTransition.cpp");
source("./Auxiliary_Functions_for_Figure/aux_PhenotypicScore.R");
source("./Auxiliary_Functions_for_Figure/aux_ggplot2Theme.R");
adjmat=readRDS("./Data_Used/AdjacentMatrix.Matrix.rds");
# Load four clusters of significant cells.
extreme=readRDS("./Intermediate_Result/ExtremeCase_FourClusters.rds");

# Set the marker's list.
Marker_List=list(
  mark_p=match(intersect(c('Esrrb','Nanog',"Klf4","Pou5f1","Sox2","Myc","Sall4"),
    colnames(adjmat)),colnames(adjmat))-1,
  mark_s=match(intersect(c('Jun','Gata6','Gata4','Bmp2','Bmp4','Cebpa',"Sox17","Sox7","Cdx2",'Foxa1','Foxa2'),
    colnames(adjmat)),colnames(adjmat))-1,
  mark_e=match(intersect(c('Cdh1','Tcf3','Ovol2','miR200','miR34'),
    colnames(adjmat)),colnames(adjmat))-1,
  mark_m=match(intersect(c('Vim','Cdh2','Zeb1','Zeb2','Foxc2','Snai1','Snai2','Twist1','Twist2','miR9'),
    colnames(adjmat)),colnames(adjmat))-1); 

# Function: Output log10-heatmap of P/S_δ vs Tr(x,ε)
fig_PressurePhone<-function(oekdss,Type='s2p',Show.Legend=FALSE){# oekdss,Filename,Type='s2p'
  # Obtain the delta of two phenotypic scores.
  tmps=cbind(oekdss$New_PS-oekdss$Ori_PS,oekdss$New_EM-oekdss$Ori_EM);
  if('s2p'==Type){
    ratio=sum(tmps[,1]>=6)/length(oekdss$Pressure);# Here we only cosider the PS score.
    xx=factor(tmps[,1],levels=seq(-5,15,1));
    yy=factor(tmps[,2],levels=seq(-5,15,1));
  }else {
    ratio=sum(tmps[,1]<=-5)/length(oekdss$Pressure);
    xx=factor(tmps[,1],levels=seq(-15,5,1));
    yy=factor(tmps[,2],levels=seq(-15,5,1));}
  # Generate the dataframe.
  tmpDT=cbind.data.frame(xx,yy);
  colnames(tmpDT)=c("xx","yy");
  tmpDT=table(tmpDT$xx,tmpDT$yy)
  tmpDT=as.data.frame(tmpDT);
  tmpDT$Var1=as.numeric(as.vector(tmpDT$Var1));
  tmpDT$Var2=as.numeric(as.vector(tmpDT$Var2));
  tmpDT$Freq=log((1+tmpDT$Freq),base = 10);
  # Output the figure 
  figs=ggplot(tmpDT)+
    geom_tile(aes(x=Var1,y=Var2,fill=Freq),color=NA)+
    scale_fill_viridis_c(limits=c(0,6))+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0)); 
  if(Show.Legend){# Have legend
    figs=Plot_ThemeConfiguration(figs)+theme(legend.position = "right");
  } else {# no legend
    figs=Plot_ThemeConfiguration(figs)+theme(legend.position = "none");
  }
  if('s2p'==Type){
    figs=figs+annotate("text",x=0,y=13,label=as.character(ratio))+
      coord_fixed(ratio = 1.0);
  } else {
    figs=figs+annotate("text",x=-5,y=-10,label=as.character(ratio))+
      coord_fixed(ratio = 1.0);
  }
  #ggsave(paste0("/dev/shm/",Filename,".pdf"), figs,width=3.0, height=2.5,units="in");
  return(figs);
}


# Perturbation of some critical genes
oekd.pressure.01=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Pou5f1'),colnames(adjmat)),c(1),503321L);
oekd.pressure.02=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Sox2'),colnames(adjmat)),c(1),503322L);
oekd.pressure.03=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Klf4'),colnames(adjmat)),c(1),503323L);
oekd.pressure.04=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Myc'),colnames(adjmat)),c(1),503324L);
oekd.pressure.05=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Esrrb'),colnames(adjmat)),c(1),503325L);
oekd.pressure.06=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Sall4'),colnames(adjmat)),c(1),503326L);
oekd.pressure.07=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Glis1'),colnames(adjmat)),c(1),503327L);
oekd.pressure.08=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Nanog'),colnames(adjmat)),c(1),503328L);

ips.c1=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Pou5f1','Sox2'),colnames(adjmat)),c(1,1),771901L);
ips.c2=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Klf4','Pou5f1'),colnames(adjmat)),c(1,1),771902L);
ips.c3=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Sox2','Klf4'),colnames(adjmat)),c(1,1),771903L);
ips.c4=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Pou5f1','Sox2','Klf4'),colnames(adjmat)),c(1,1,1),771904L);
ips.c5=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Pou5f1','Esrrb'),colnames(adjmat)),c(1,1),771905L);
ips.c6=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Pou5f1','Nanog'),colnames(adjmat)),c(1,1),771906L);
ips.c7=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Esrrb','Nanog'),colnames(adjmat)),c(1,1),771907L);
ips.c8=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Esrrb','Nanog','Sall4'),colnames(adjmat)),c(1,1,1),771908L);
  
fig_all_perturbated=gridExtra::grid.arrange(
    fig_PressurePhone(oekd.pressure.01),# Pou5f1
    fig_PressurePhone(oekd.pressure.02),# Sox2
    fig_PressurePhone(oekd.pressure.03),# Klf4
    fig_PressurePhone(oekd.pressure.04),# Myc
    fig_PressurePhone(oekd.pressure.05),# Esrrb
    fig_PressurePhone(oekd.pressure.06),# Sall4
    fig_PressurePhone(oekd.pressure.07),# Glis1
    fig_PressurePhone(oekd.pressure.08),# Nanog
    fig_PressurePhone(ips.c1),# Pou5f1 + Sox2
    fig_PressurePhone(ips.c2),# Pou5f1 + Klf4
    fig_PressurePhone(ips.c3),# Sox2 + Klf4
    fig_PressurePhone(ips.c4),# Pou5f1 + Sox2 + Klf4
    fig_PressurePhone(ips.c5),# Pou5f1 + Esrrb
    fig_PressurePhone(ips.c6),# Pou5f1 + Nanog
    fig_PressurePhone(ips.c7),# Esrrb + Nanog
    fig_PressurePhone(ips.c8),# Esrrb + Nanog + Sall4
  nrow=4,ncol=4);

ggsave("./Figure_OriginalVersion/PerturbatedCriticalGene_All_2pheno.pdf",
    fig_all_perturbated, width=15.0, height=12.0, units="in");

ggsave("./Figure_OriginalVersion/PerturbatedCriticalGene_One.pdf",
    fig_PressurePhone(oekd.pressure.01,Show.Legend=TRUE), width=4.0, height=4.0, units="in");
# Recommend!
saveRDS(oekd.pressure.01,file="./Intermediate_Result/oekd_Pouf51_sample.rds",compress=FALSE);

rm(list=ls());gc();
# Code is over.