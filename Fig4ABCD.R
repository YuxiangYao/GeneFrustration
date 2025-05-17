# Figures: OE/KD reshape cell's phenotypic score in various scaling law patterns.

# Load and complie data, functions. 
Rcpp::sourceCpp("./cpp_Source/c_OEKD_Testing.cpp");
source("./Auxiliary_Functions_for_Figure/aux_ggplot2Theme.R");
# source("./Auxiliary_Functions_for_Figure/aux_LogLogBin.R");
sss=readRDS("./Intermediate_Result/SSS_origin_nodup.rds");
adjmat=readRDS("./Data_Used/AdjacentMatrix.Matrix.rds");
keepname=intersect(colnames(adjmat),c(
  c('Cdh1','Tcf3','Ovol2','miR200','miR34'),#EM-Positive
  c('Vim','Cdh2','Zeb1','Zeb2','Foxc2','Snai1','Snai2','Twist1','Twist2','miR9'),#EM-Negative
  c('Esrrb','Nanog',"Klf4","Pou5f1","Sox2","Myc","Sall4"),#SP-Positive
  c('Jun','Gata6','Gata4','Bmp2','Bmp4','Cebpa',"Sox17","Sox7","Cdx2",'Foxa1','Foxa2')));#SP-Negative


# Simulate all genes perturbations, return OE_Δ cases, KD_Δ cases, OE_δ, KD_δ marker-states.
# Get nearly one hour. (Calculate, too large, too slow, Recommended background operation)
oe.kd=OEKD_Testing(adjmat,sss,match(keepname,colnames(adjmat)),821992L);
colnames(oe.kd$OE_Diff)=colnames(oe.kd$KD_Diff)=colnames(adjmat);

# Figure 4A: OE/KD changed state distribution (Δs). 
oe.kd$OE_Diff->tmp;tmp=as.vector(tmp);tmp=tmp[tmp>0];tab.oe=table(tmp);# Only non-zero cases.
oe.kd$KD_Diff->tmp;tmp=as.vector(tmp);tmp=tmp[tmp>0];tab.kd=table(tmp);
OEs=data.frame(x=as.numeric(names(tab.oe)), y=as.matrix(tab.oe), c="OE");
KDs=data.frame(x=as.numeric(names(tab.kd)), y=as.matrix(tab.kd), c="KD");
alldata=rbind.data.frame(KDs,OEs);
alldata$x=log2(alldata$x);
alldata$y=log2(alldata$y);
alldata$c=factor(alldata$c,levels=c('OE','KD'));# Assign the order.
coeffs=cbind.data.frame(x=log2(c(1:16)),y=rowMeans(cbind(log2(OEs[1:16,2]),log2(KDs[1:16,2]))))
coeffs=lm(coeffs$y ~ coeffs$x)$coefficients;# Instead of bin-log, here we use fitting line to guide.
fitline=cbind.data.frame(seq(0.5,3.7,by=0.1),seq(0.5,3.7,by=0.1)*coeffs[2]+coeffs[1]-1.5);
colnames(fitline)=c('x','y');
fig_oekd=ggplot()+
  geom_point(data=alldata,aes(x=x,y=y,col=c,shape=c),size=2.0)+
  scale_color_manual(values = c('#EE312E','#2156A6','#EE312E','#2156A6'))+
  scale_shape_manual(values = c(16,16,1,1))+
  geom_line(data=fitline,aes(x=x,y=y),linewidth=1.5,col="#222222");
fig_oekd=Plot_ThemeConfiguration(fig_oekd)+
  scale_y_continuous(limits=c(-0.5,30),breaks=seq(0,30,5), 
    labels=c(latex2exp::TeX("${2^0}$"),latex2exp::TeX("${2^5}$"),
      latex2exp::TeX("$2^{10}$"),latex2exp::TeX("$2^{15}$"),latex2exp::TeX("$2^{20}$"),
      latex2exp::TeX("$2^{25}$"),latex2exp::TeX("$2^{30}$")))+
  scale_x_continuous(limits=c(-0.5,6),breaks=c(0:6), labels=c(1,2,4,8,16,32,64))+
  annotate("text",x=1.7,y=20,label="-2.72",size=3)+ # Write slope
  theme(legend.position =c(.25,.30),aspect.ratio=0.6,
    legend.background = element_rect(fill="transparent"));
ggsave("./Figure_OriginalVersion/OEKD_ScalingLaw_GeneticChange.pdf",
  fig_oekd,width=4.5,height=4.5,units="in");

# Show gene' names. (for Table Gene distribution)
# OE case >=32
oe.kd$OE_Diff->tmp.oe;
tmp.oe=(tmp.oe>=32);
oe_scenarios=colSums(tmp.oe);
oe_scenarios=oe_scenarios[oe_scenarios>0];
oe_scenarios=oe_scenarios[order(oe_scenarios, decreasing=T)];
# OE case >=32
oe.kd$KD_Diff->tmp.kd;
tmp.kd=(tmp.kd>=32);
kd_scenarios=colSums(tmp.kd);
kd_scenarios=kd_scenarios[kd_scenarios>0];
kd_scenarios=kd_scenarios[order(kd_scenarios, decreasing=T)];
saveRDS(list(oe_scenarios=oe_scenarios, kd_scenarios=kd_scenarios),
  file="./Intermediate_Result/Tab_LargerThan32Cases.rds",compress = FALSE);
#rm(tmp.oe, tmp.kd, oe_scenarios, kd_scenarios);




# Figure：top utility genes of Poisson Distribution & fitting show.
PoissonFittingNonZero<-function(vec){
  vec=vec[vec>0&vec<20];# avoid extreme value >=20
  poi.fit=MASS::fitdistr(vec,"poisson");
  return(c(poi.fit$estimate, poi.fit$sd, max(vec)));
}
lambdaset=rbind(apply(oe.kd$OE_Diff,2,PoissonFittingNonZero),
                apply(oe.kd$KD_Diff,2,PoissonFittingNonZero));# ~10min
lambdaset=t(lambdaset);
rownames(lambdaset)=colnames(adjmat);

lambda.p=cbind.data.frame(x=rep(colnames(adjmat),time=2),
  y=c(lambdaset[,1],lambdaset[,4]),c=rep(c('OE','KD'),each=88));
lambda.p$x=factor(lambda.p$x,levels=colnames(adjmat)[order((lambdaset[,1]+lambdaset[,4]),decreasing = T)]);
lambda.p$c=factor(lambda.p$c,levels = c('OE','KD'));
lambda.p=rbind.data.frame(lambda.p[89:176,],lambda.p[1:88,]);

# sub-insert figures. (Figure 4CD inserts)
oekd_subinsert<-function(OE.KD,targets,Color){
  Vec2TabMat<-function(vec){
    res=vec;
    res=table(res[res>0]);
    res=cbind(as.numeric(names(res)),res);
    colnames(res)=c('x','y');
    return(res);
  }
  poireal=rbind.data.frame(
    cbind.data.frame(Vec2TabMat(OE.KD$OE_Diff[,targets]),c='aOE'),
    cbind.data.frame(Vec2TabMat(OE.KD$KD_Diff[,targets]),c='bKD') )
  poireal=poireal[poireal[,1]<=15,];
  index=(poireal$c=='bKD');poireal[index,2]=poireal[index,2]/sum(poireal[index,2]);
  index=(poireal$c=='aOE');poireal[index,2]=poireal[index,2]/sum(poireal[index,2]);
  c.name=rep(colnames(adjmat),2);
  index1=lambda.p[c.name==targets&lambda.p[,3]=='OE',2];
  index2=lambda.p[c.name==targets&lambda.p[,3]=='KD',2];
  x.data=seq(0,15,0.1);
  poifit=rbind.data.frame(
    cbind.data.frame(x=x.data,c='aOE',y=index1^x.data*exp(-index1)/factorial(x.data)),
    cbind.data.frame(x=x.data,c='bKD',y=index2^x.data*exp(-index2)/factorial(x.data)));
  fig_tmp=ggplot()+
    geom_point(data=poireal,aes(x=x,y=y,shape=c),col=Color,size=2.5,stroke=0.25)+
    geom_line(data=poifit,aes(x=x,y=y,linetype=c),col=Color,linewidth=1.5)+
    scale_shape_manual(values=c(16,2))+
    scale_y_continuous(expand=c(0,0.05))+
    scale_x_continuous(expand=c(0.05,0),limits=c(0,15),breaks=seq(0,15,5),labels=seq(0,15,5));
  fig_tmp=Plot_ThemeConfiguration(fig_tmp)+theme(legend.position=c(.66,.70),aspect.ratio=0.45,);
  return(fig_tmp);
}
# Output: 
ggsave(filename="./Figure_OriginalVersion/OEKD_LambdaFit1_Pou5f1_in.pdf", width=4.5, height=4.5,
  plot=oekd_subinsert(oe.kd,"Pou5f1","#2156A6"), units="in");
ggsave(filename="./Figure_OriginalVersion/OEKD_LambdaFit2_Sox17_in.pdf", width=4.5, height=4.5,
  plot=oekd_subinsert(oe.kd,"Sox17","#EE312E"), units="in");


# Figure 4B
tmpdata=readRDS("./Intermediate_Result/Tab_LargerThan32Cases.rds")
oe_scenarios=tmpdata[[1]];
kd_scenarios=tmpdata[[2]];

genenames1=c("Nanog","Sox2","Pou5f1","Esrrb","Cebpa",
             "Tgfb1","Klf4","Smad4","Sall4","Prdm14",
             "Smad2","Nodal","Tcl1","miR101","Bmp2",
             "Zeb1","Snai1","Jdp2");
genenames2=c("Tcf7l1","Sox17","Mbd3","Mapk1","Jun",
             "miR145","Mbd2","Trp53","Twist2","Gata6");
# case1=c(8221,4126,1609,323,222,
#         147,123,98,96,51, 47,45,32,23,20,18,16,12);
case1=oe_scenarios[oe_scenarios>10]

# case2=c(1970,865,571,307,271,
#         211,42,33,24,11);
case2=kd_scenarios[kd_scenarios>10]

data_cases=cbind.data.frame(
  x=c(genenames1,genenames2),
  y=c(case1,case2),
  t=c(rep("OE",18),rep("KD",10))
)
data_cases$x=factor(data_cases$x,levels=c(genenames1,genenames2));

# OE/KD lambdas of each gene. (Figure 4B)
fig_lambda3=
  ggplot()+geom_bar(data=data_cases,aes(x=x,y=log2(y),fill=t),
                    stat="identity",position="stack",alpha=0.75, width=1.0,linewidth=0)+
  scale_fill_manual(values=c("#2156A6","#EE312E"))+
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  scale_x_discrete(expand=c(0.01,0),)+
  scale_y_continuous(limits=c(0,14),breaks=seq(1,13,2),expand=c(0,0),labels=seq(1,13,2));
fig_lambda3=Plot_ThemeConfiguration(fig_lambda3)+
  theme(legend.position =c(.85,.75),aspect.ratio=0.60,
        axis.text.x=element_text(size=12,vjust=0.5,hjust=1,angle=90,face="bold",));
ggsave("./Figure_OriginalVersion/OEKD_LambdaFit0_All_Genes_4B.pdf",
       fig_lambda3, width=7.5, height=4.5, units="in");

# Figure:top utility genes of Poisson Distribution & fitting show.
# Figure 4C
# PoissonFittingNonZero<-function(vec){
#   vec=vec[vec>0&vec<20];# avoid extreme value >=20
#   poi.fit=MASS::fitdistr(vec,"poisson");
#   return(c(poi.fit$estimate, poi.fit$sd, max(vec)));
# }
# lambdaset=rbind(apply(oe.kd$OE_Diff,2,PoissonFittingNonZero),
#                 apply(oe.kd$KD_Diff,2,PoissonFittingNonZero));# ~10min
# lambdaset=t(lambdaset);
# rownames(lambdaset)=colnames(adjmat);
lambda.p=cbind.data.frame(x=rep(colnames(adjmat),time=2),
                          y=c(lambdaset[,1],lambdaset[,4]),c=rep(c('OE','KD'),each=88));
lambda.p$x=factor(lambda.p$x,levels=colnames(adjmat)[order((lambdaset[,1]),decreasing = T)]);
lambda.p$c=factor(lambda.p$c,levels = c('OE','KD'));
lambda.p=rbind.data.frame(lambda.p[89:176,],lambda.p[1:88,]);
keepname=colnames(adjmat)[(lambdaset[,1])>1.0001];# Least one the effects.
# OE/KD lambdas of each gene.
fig_lambda=
  ggplot()+geom_bar(data=lambda.p[lambda.p$x%in%keepname&lambda.p$c=="OE",],aes(x=x,y=y,fill=c),
                    stat="identity",position="stack",alpha=0.75, width=1.0,linewidth=0)+
  scale_fill_manual(values=c("#aaaaaa","#666666"))+
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  scale_x_discrete(expand=c(0.01,0),)+
  scale_y_continuous(limits=c(0,8.5),breaks=seq(1,8,1),expand=c(0,0),labels=seq(1,8,1));
colorlink=rep("#000000",88);# Set x-axis's label color.
colorlink[1]='#2156A6';colorlink[4] ='#EE312E';
fig_lambda=Plot_ThemeConfiguration(fig_lambda)+
  theme(legend.position =c(.85,.75),aspect.ratio=0.35,
        axis.text.x=element_text(size=12,vjust=0.5,hjust=1,angle=90,face="bold",colour=colorlink,));
ggsave("./Figure_OriginalVersion/OEKD_LambdaFit0_All_Genes_4C.pdf",
       fig_lambda, width=11.0, height=4.5, units="in");

# Figure 4D
lambda.p$x=factor(lambda.p$x,levels=colnames(adjmat)[order((lambdaset[,4]),decreasing = T)]);
lambda.p$c=factor(lambda.p$c,levels = c('OE','KD'));
lambda.p=rbind.data.frame(lambda.p[89:176,],lambda.p[1:88,]);
keepname=colnames(adjmat)[(lambdaset[,4])>1.0001];# Least one the effects.
# OE/KD lambdas of each gene.
fig_lambda2=
  ggplot()+geom_bar(data=lambda.p[lambda.p$x%in%keepname&lambda.p$c=="KD",],aes(x=x,y=y,fill=c),
                    stat="identity",position="stack",alpha=0.75, width=1.0,linewidth=0)+
  scale_fill_manual(values=c("#aaaaaa","#666666"))+
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  scale_x_discrete(expand=c(0.01,0),)+
  scale_y_continuous(limits=c(0,3.5),breaks=seq(1,3,1),expand=c(0,0),labels=seq(1,3,1));
colorlink=rep("#000000",88);# Set x-axis's label color.
colorlink[1]='#2156A6';colorlink[4] ='#EE312E';
fig_lambda2=Plot_ThemeConfiguration(fig_lambda2)+
  theme(legend.position =c(.85,.75),aspect.ratio=0.35,
        axis.text.x=element_text(size=12,vjust=0.5,hjust=1,angle=90,face="bold",colour=colorlink,));
ggsave("./Figure_OriginalVersion/OEKD_LambdaFit0_All_Genes_4D.pdf",
       fig_lambda2, width=11.0, height=4.5, units="in");

# Code is over.