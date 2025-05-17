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
# mid.sp=0.5*(max(x.sp)+min(x.sp)); Middle S-P part.
# mid.me=0.5*(max(x.me)+min(x.me)); Middle M-E part.
# mid.mix=intersect(order(abs(mid.me-x.me))[1:top.20],order(abs(mid.sp-x.sp))[1:top.20]); Middle part.

# Figure: Signficantly phenotypic cell in the landscape 
pca_index=readRDS("./Intermediate_Result/PCA_index_coor.rds");# Previous processing
datas2=cbind.data.frame(pca_index[,2:3],"CC");
colnames(datas2)=c("x","y","group");
datas2[index.pe,3]="PE";
datas2[index.se,3]="SE";
datas2[index.pm,3]="PM";
datas2[index.sm,3]="SM";
datas2=datas2[order(datas2[,3]),];# Display trivial cells on the bottom.
fig_extre_dist=ggplot()+
  geom_point(data=datas2,aes(x=x, y=y, col=group),size=0.1,alpha=1.00)+
  scale_color_manual(values=c('#c3c3c3','#d64e3d','#4087af','#eec67f','#9a91c8'));
fig_extre_dist=Plot_ThemeConfiguration(fig_extre_dist)+
  theme(axis.line.x.bottom=element_blank(), axis.line.y.left=element_blank(),
    axis.text.x.bottom=element_blank(), axis.text.y.left=element_blank(),
    axis.title.x.bottom=element_blank(), axis.title.y.left=element_blank(),
    legend.position="none", aspect.ratio=1,
    legend.background=element_rect(fill="transparent"));
# Many dots, use png format. Other information added manually.
ggsave("./Figure_OriginalVersion/SignificantPhenotypic_Location.png",fig_extre_dist,width=4.5,height=4.5,units="in");
# Show detail numbers:
cat("PE:", length(index.pe),"\n", "PM:", length(index.pm),"\n",
    "SE:", length(index.se),"\n", "SM:", length(index.sm),"\n");

rm(list=ls());gc();
# Code is over.
