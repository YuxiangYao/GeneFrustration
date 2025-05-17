# Figures of initial pseudo-energy of success/failure cells.
# Output: 2 figures.
# InitialEnergy_UpperLower InitialEnergy_Redo

# Load data, function.
Rcpp::sourceCpp("./cpp_Source/c_OEKD_StateTransition.cpp");
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

Case_Pou5f1=Pressure2Potential(adjmat,extreme$sm,Marker_List,
  match(c('Pou5f1'),colnames(adjmat)),c(1),503321L);
# Before execute this file, [Fig_OEKD_transition] has been done.
# Thus the temporary file [oekd_Pouf51_sample.rds] exists. 
# or do this: -->  Case_Pou5f1=readRDS("./Intermediate_Result/oekd_Pouf51_sample.rds");
Index_High=Case_Pou5f1$New_PS-Case_Pou5f1$Ori_PS>6;
Index_Lows=Case_Pou5f1$New_PS-Case_Pou5f1$Ori_PS<6;
# Redo each part.
Case_Pou5f1_upP=Pressure2Potential(adjmat,extreme$sm[Index_High,], Marker_List,
  match(c('Pou5f1'),colnames(adjmat)),c(1),779843L);
Case_Pou5f1_lowP=Pressure2Potential(adjmat,extreme$sm[Index_Lows,], Marker_List,
  match(c('Pou5f1'),colnames(adjmat)),c(1),661029L);


# Funciton: Count and fit energy distribution.
OriEngery.count.fitting<-function(SamMat,Label){
  sam.diff=SamMat;
  tmp.fit=MASS::fitdistr(SamMat,"normal");
  muss=tmp.fit$estimate[1];
  tmp.fit=cbind.data.frame(seq(min(sam.diff),max(sam.diff),length.out=500),
    dnorm(seq(min(sam.diff),max(sam.diff),length.out=500),log=FALSE,
    mean=tmp.fit$estimate[1],sd=tmp.fit$estimate[2]),Label);
  colnames(tmp.fit)=c("x","y","c");
  sam.diff=table(sam.diff);
  sam.diff=data.frame(
    x=as.integer(rownames(as.matrix(sam.diff))),
    y=as.matrix(sam.diff));
  sam.diff$y=sam.diff$y/sum(sam.diff$y);
  return(list(sam.diff,tmp.fit,muss));# 1-bar;2-fitted line
}
# Function: Generate barplot and fitting line.
source("./Auxiliary_Functions_for_Figure/aux_ggplot2Theme.R");
BarplotFittingLine<-function(CaseList, yLims, yLabel,Bins){
  High_ID=CaseList$New_PS-CaseList$Ori_PS>6;
  Low_ID=CaseList$New_PS-CaseList$Ori_PS<6;
  hihi1=OriEngery.count.fitting(CaseList$OriEng[High_ID],"High");
  hihi2=OriEngery.count.fitting(CaseList$OriEng[Low_ID],"Low" );
  diff.all=rbind.data.frame(hihi2[[2]],hihi1[[2]]);
  diff.all$c=factor(diff.all$c,levels = c('Low','High'));
  oriEng1=cbind.data.frame(x=CaseList$OriEng[High_ID],type="High");
  oriEng2=cbind.data.frame(x=CaseList$OriEng[Low_ID],type="Low");
  fig_s=ggplot()+ 
    geom_histogram(data=oriEng2, 
      aes(x = x,y = after_stat(count)/sum(after_stat(count)) ), 
      position = "identity", alpha = 0.55, fill="#EE312E",binwidth = Bins) + 
    geom_histogram(data=oriEng1, 
      aes(x = x,y = after_stat(count)/sum(after_stat(count)) ), 
      position = "identity", alpha = 0.55, fill="#2156A6",binwidth = Bins) + 
    geom_line(data=diff.all,aes(x=x,y=y,col=c),linewidth=1.2)+
    scale_color_manual(values = c("#EE312E","#2156A6"))+
    labs(x="Initial E(x)",y="Density")+
    scale_y_continuous(expand=c(0,0),limits=yLims,
      breaks=yLabel,labels=yLabel)+
    scale_x_continuous(expand=c(0,0),limits=c(-225,0),
      breaks=seq(-200,0,50),labels=c('-200','-150','-100','-50','0'))+
    annotate("text", x=-200, y=0.015, label=hihi2[[3]],
      hjust=0, vjust=0, size=5, colour="#000000")+
    annotate("text", x=-25, y=0.015, label=hihi1[[3]],
      hjust=0, vjust=0, size=5, colour="#000000");
  fig_s=Plot_ThemeConfiguration(fig_s)+
    theme(legend.position=c(.1,.65), aspect.ratio=1,);
  return(fig_s);
}

# Figure: Original Energy of Upper and Lower clusters.
ggsave(filename="./Figure_OriginalVersion/InitialEnergy_UpperLower.pdf",
  plot=BarplotFittingLine(Case_Pou5f1, c(0,0.03), seq(0,0.03,0.01),2), width=5.0, height=4.0);
fig_redo=ggpubr::ggarrange(
    BarplotFittingLine(Case_Pou5f1_lowP, c(0,0.08), seq(0,0.08,0.02),4),
    BarplotFittingLine(Case_Pou5f1_upP, c(0,0.08), seq(0,0.08,0.02),4));
ggsave(filename="./Figure_OriginalVersion/InitialEnergy_Redo.pdf",
  plot=fig_redo, width=12, height=4.0);

# Check similarity of two clusters.
Case_Pou5f1_Redo=Pressure2Potential(adjmat, extreme$sm, Marker_List,
  match(c('Pou5f1'), colnames(adjmat)), c(1), 613321L);

New_Index_High=Case_Pou5f1_Redo$New_PS-Case_Pou5f1_Redo$Ori_PS>6;# Obtain index.
New_Index_Lows=Case_Pou5f1_Redo$New_PS-Case_Pou5f1_Redo$Ori_PS<6;

ind101=c(1:1000000)[Index_High];
ind102=c(1:1000000)[Index_Lows];
ind201=c(1:1000000)[New_Index_High];
ind202=c(1:1000000)[New_Index_Lows];

jaccard_similarity <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(intersection / union)
}
cat("Jaccard Index:\n");
matrix(
  c(jaccard_similarity(ind101,ind201),
    jaccard_similarity(ind101,ind202),
    jaccard_similarity(ind102,ind201),
    jaccard_similarity(ind102,ind202)),
  2,2,byrow=TRUE);

cat("Cell Number:\n");
matrix(
  c(NA, length(ind201),length(ind202),
    length(ind101),length(intersect(ind101,ind201)),length(intersect(ind101,ind202)),
    length(ind102),length(intersect(ind102,ind201)),length(intersect(ind102,ind202))),
  3,3,byrow=TRUE);

# Code is over.