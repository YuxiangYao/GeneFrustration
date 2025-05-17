# Figure: Signficantly phenotypic cell's location and proportions.

# Loaded data & Necessary functions.
sss=readRDS("./Intermediate_Result/SSS_origin_nodup.rds");
index.nonrep=readRDS("./Intermediate_Result/nonduplicate_index.rds");
Steps=readRDS("./Intermediate_Result/Step_Energy.rds");
source("./Auxiliary_Functions_for_Figure/aux_PhenotypicScore.R");
source("./Auxiliary_Functions_for_Figure/aux_ggplot2Theme.R");

# Check phenotypic feature cells.
top.20=as.integer(dim(sss)[1]*0.20);# Top 20% cell as 
x.sp=SPT_score(sss);# Score of P/S
x.me=EMT_score(sss);# Score of E/M
rank.sp=order(x.sp);# from S to P. (negative -> positive)
rank.me=order(x.me);# from M to E. (negative -> positive)
index.pe=intersect(rev(rank.sp)[1:top.20],rev(rank.me)[1:top.20]);# P+E type
index.pm=intersect(rev(rank.sp)[1:top.20],rank.me[1:top.20]);# P+M type
index.se=intersect(rank.sp[1:top.20],rev(rank.me)[1:top.20]);# S+E type
index.sm=intersect(rank.sp[1:top.20],rank.me[1:top.20]);# S+M type
# Obtain the corresponding cells.
pe.sss=sss[index.pe,];
sm.sss=sss[index.sm,];
pe.activated.ratio=colMeans(pe.sss);
sm.activated.ratio=colMeans(sm.sss);
# Prepare for data-frame.
library(ggplot2);
delta.ratio=pe.activated.ratio-sm.activated.ratio;
Orders=order(abs(delta.ratio));
av.gene.change=cbind.data.frame(x=c(1:length(delta.ratio)),
  n=names(delta.ratio)[Orders],y=delta.ratio[Orders],c="up");
rownames(av.gene.change)=names(delta.ratio)[Orders];
av.gene.change$c[av.gene.change$y<0]="down";
av.gene.change$c=factor(av.gene.change$c,levels=c("up","down"));
av.gene.change$n=factor(av.gene.change$n,levels=rownames(av.gene.change));

x.lab.col=rep("#000000",length(delta.ratio));
x.lab.col[av.gene.change$y<(-0.7)]="#EE312E";
x.lab.col[av.gene.change$y>(0.7)]="#2156A6";
x.line.col=x.lab.col[x.lab.col!="#000000"]
x.line.col[x.line.col=="#EE312E"]="#EE312E30";
x.line.col[x.line.col=="#2156A6"]="#2156A630";
library(ggplot2);

figs=ggplot()+
  geom_hline(yintercept = 0.7, color = "#acacac", linetype = "dashed")+
  geom_hline(yintercept =-0.7, color = "#acacac", linetype = "dashed")+
  geom_vline(xintercept = c((length(x.lab.col)-length(x.line.col)+1):
                            length(x.lab.col)),
             color=x.line.col,linetype = "solid")+
  geom_point(data=av.gene.change,aes(x=x,y=y,col=c,shape=c),size=2.5)+
  scale_color_manual(values=c("#2156A6","#EE312E"))+#"white",
  scale_shape_manual(values=c(15,16))+
  scale_x_continuous(expand = c(0,0),breaks =seq(1,88,1),
                     labels=av.gene.change$n,limits=c(0,89))+
  scale_y_continuous(expand = c(0,0),breaks =c(-1.0,-0.5,0,0.5,1.0),
                     labels=c(-1.0,-0.5,0,0.5,1.0),limits=c(-1.1,1.1))+
  labs(x="Gene",y="Expression Ratio of PE-SM");
figs=Plot_ThemeConfiguration(figs)+
  theme(axis.text.x.bottom=element_text(angle = 90,hjust = 1.0,colour = x.lab.col));
#figs
ggsave("./Figure_OriginalVersion/Fig_PhenoPhases_Observe_SM2PE.pdf",figs,height = 4.0,width=12.0,units="in")

hvgs=as.character(av.gene.change$n[abs(av.gene.change$y)>0.7]);
lvgs=as.character(av.gene.change$n[abs(av.gene.change$y)<0.7]);
var.gene=list(hvgs=hvgs,lvgs=lvgs);
saveRDS(var.gene,"./Intermediate_Result/var_gene_PE2SM.rds");

rm(list=ls());gc();
# Code is over.