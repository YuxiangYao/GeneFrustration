
library("ggplot2");# Need reshape2, Rcpp;
source("../Auxiliary_Functions_for_Figure/aux_ggplot2Theme.R");
Rcpp::sourceCpp("../cpp_Source/c_DoubleBinCountAv.cpp");# C++: DoubleBins(...)
BinnedPCA_DataFrame<-function(PCAss,Values,IDs=c(2,3),Paths=FALSE){
  zuobiao=cbind.data.frame(PCAss[,IDs]);
  bin=40L;
  if(identical(IDs,c(2,3))){# PC2 + PC3
    x_bian=c(-2.5,4.0);
    y_bian=c(-2.5,3.0);}
  else {# PC1 + PC2
    x_bian=c(-8.2,-2.5);
    y_bian=c(-2.5,4.0);}
  seq_x=seq(x_bian[1],x_bian[2],length.out=(bin+1));
  seq_y=seq(y_bian[1],y_bian[2],length.out=(bin+1));
  bin_x=floor((zuobiao[,1]-x_bian[1])/(seq_x[2]-seq_x[1]));# Gongcha-X
  bin_y=floor((zuobiao[,2]-y_bian[1])/(seq_y[2]-seq_y[1]));# Gongcha-Y
  # Binning ....
  if(Paths){
    res=DoubleCount(bin_x,bin_y,Values,bin,bin);
  }else {
    res=DoubleBins(bin_x,bin_y,Values,bin,bin);}
  colnames(res)=rownames(res)=paste0("x_",c(1:bin));
  res=reshape2::melt(res);
  if(Paths){
    ;# No need operation
  }else {
    res$value[res$value>-1]=NA;# NULL set as NA.
  }
  return(res);
}

BinnedPCA_Plot<-function(PCAss,Values,IDs=c(2,3),xLabs,yLabs){
  tmp_df=BinnedPCA_DataFrame(PCAss,Values,IDs);
  figs=ggplot(tmp_df)+
    geom_tile(aes(x=Var1,y=Var2,fill=value))+
    labs(x=paste0("PC",IDs[1]),y=paste0("PC",IDs[2]))+
    scale_fill_viridis_c(limits=c(-225,-35),# <- this value is original limited values.
      name="Pseudo-\nenergy", option = "magma", na.value = "#ffffff")+
    scale_x_discrete(expand=c(0,0), breaks=paste0("x_",c(1,9,17,25,33,40)),
      labels=xLabs)+
    scale_y_discrete(expand=c(0,0), breaks=paste0("x_",c(1,9,17,25,33,40)),
      labels=yLabs);
  figs=Plot_ThemeConfiguration(figs)+theme(legend.position="right",aspect.ratio = 1.0);
  return(figs);
}

BinnedPCA_Diff_Plot<-function(PCAss1,Values1,PCAss2,Values2,IDs=c(2,3),xLabs,yLabs){
  df_1=BinnedPCA_DataFrame(PCAss1,Values1,IDs);
  df_2=BinnedPCA_DataFrame(PCAss2,Values2,IDs);
  dfss=cbind.data.frame(df_1[,1:2],value=df_2$value-df_1$value);# <- obtain changed energy.
  figs=ggplot(dfss)+
    geom_tile(aes(x=Var1,y=Var2,fill=value))+
    labs(x=paste0("PC",IDs[1]),y=paste0("PC",IDs[2]))+
    scale_fill_gradient2(low="#2156A6", high="#EE312E",mid="#dadada",midpoint =0,
      limits=c(-75,75),# <- this value is original limited values.
      na.value="#ffffff")+
    scale_x_discrete(expand=c(0,0), breaks=paste0("x_",c(1,9,17,25,33,40)),
      labels=xLabs)+
    scale_y_discrete(expand=c(0,0), breaks=paste0("x_",c(1,9,17,25,33,40)),
      labels=yLabs);
  figs=Plot_ThemeConfiguration(figs)+theme(legend.position="right",aspect.ratio = 1.0);
  return(figs);
}

BinnedPCA_PathDensity<-function(PCAss,Values,IDs=c(2,3),xLabs,yLabs){
  tmp_df=BinnedPCA_DataFrame(PCAss,Values,IDs,TRUE);
  figs=ggplot(tmp_df)+
    geom_tile(aes(x=Var1,y=Var2,fill=log1p(value)))+
    labs(x=paste0("PC",IDs[1]),y=paste0("PC",IDs[2]))+
    scale_fill_viridis_c(limits=c(0,13))+
    scale_x_discrete(expand=c(0,0), breaks=paste0("x_",c(1,9,17,25,33,40)),
      labels=xLabs)+
    scale_y_discrete(expand=c(0,0), breaks=paste0("x_",c(1,9,17,25,33,40)),
      labels=yLabs);
  figs=Plot_ThemeConfiguration(figs)+theme(legend.position="right",aspect.ratio = 1.0);
  return(figs);
}

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
    scale_color_manual(values=c("#2156A6","#EE312E","#6c6c6c","#FFB600","#00ff00"))+
    scale_fill_manual(values=c("#2156A6","#EE312E","#6c6c6c","#FFB600","#00ff00"))+
    scale_linetype_manual(values = c("solid","solid","dashed"))+
    labs(x=ID,y="pseudoenergy")+
    scale_x_continuous(breaks=c(0,10,20,30,40,50), labels=xLab, expand=c(0,0),);
  fig_tempeng=Plot_ThemeConfiguration(fig_tempeng)+theme(legend.position="right",aspect.ratio=1.0);
  return(fig_tempeng);
}
