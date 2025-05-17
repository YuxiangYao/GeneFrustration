# Figure: Analyze the large/middle/small basin's average energy.
# Output: BasinSize_ThreeAvEnergy.pdf

# Load data & necessary function:
source("./Auxiliary_Functions_for_Figure/aux_ggplot2Theme.R");
mc.attr=readRDS("./Intermediate_Result/mc_attr_origin.rds");
# Get stable states (including duplicate cases).
sss=mc.attr$StateMatrix;
# Turn to character vectors. (To compare and count the basin)
sss.nameid=apply(sss,1,paste,collapse="");

ttt=cbind.data.frame(eng=mc.attr$Steps[,2], names=sss.nameid);
cat("here lines:",nrow(ttt),"\n");

# Each isoform stable state has how much attracting area.
sss.nameid=as.matrix(table(sss.nameid));
# Obtain index of large basin sizes.
index_Large=rownames(sss.nameid)[sss.nameid[,1]>=8];
index_Small=rownames(sss.nameid)[sss.nameid[,1]==1];
index_Middle=rownames(sss.nameid)[1<sss.nameid[,1]&sss.nameid[,1]<8];
# Obtain stable states of large basin sizes (including duplicate cases).

ttt=ttt[!duplicated(ttt$names),];
cat("here lines:",nrow(ttt),"\n");

index_Large_Ori=(ttt$names)%in%index_Large;
index_Small_Ori=(ttt$names)%in%index_Small;
index_Middle_Ori=(ttt$names)%in%index_Middle;

cat("L+M+S: ",sum(index_Large_Ori),
    sum(index_Middle_Ori),sum(index_Small_Ori),"\n");

# Get the energy vector and each part of large/small 
#energy=mc.attr$Steps[,2];
energy_Large_BS=ttt$eng[index_Large_Ori];
energy_Small_BS=ttt$eng[index_Small_Ori];
energy_Middle_BS=ttt$eng[index_Middle_Ori];
bs_three_part=list(large=index_Large_Ori, middle=index_Middle_Ori, small=index_Small_Ori);

# Highly recommend! for subsequent using.
saveRDS(bs_three_part, file="./Intermediate_Result/BasinSizeIndex_LargeMiddleSmall.rds", compress=FALSE);

# Prepare for dataframe of ggplot;
datas=data.frame(parts=NA, value=ttt$eng);
datas$parts[bs_three_part$small]="Small";
datas$parts[bs_three_part$middle]="Middle";
datas$parts[bs_three_part$large]="Large";
datas$parts=factor(datas$parts,levels=c("Small","Middle","Large"));

# Generate and output figures.
library(ggplot2);
library(ggnewscale); # for more one scale slots.
figs=ggplot(data=datas,aes(x=parts, y=value, fill=parts, color=parts))+
  geom_violin(trim = FALSE)+
  scale_fill_manual(values = c("#EE312E88","#FFB60088","#2156A688"))+
  scale_color_manual(values = c("#EE312E","#FFB600","#2156A6"))+
  geom_boxplot(width = 0.1,  color = "black", outlier.shape = NA)+
  new_scale_fill() +
  scale_fill_manual(values = c("#EE312E","#FFB600","#2156A6"))+
  new_scale_color()+
  scale_color_manual(values = c("#000000","#000000","#000000"));
figs=Plot_ThemeConfiguration(figs)+
  theme(aspect.ratio=0.60,legend.position="right");

# Output the figure.
ggsave("./Figure_OriginalVersion/BasinSize_ThreeAvEnergy.pdf",figs,width=6.0,height=4.5,units="in");

# Show averages and standard deviations of three sizes of basin.
cat("Large av:", mean(energy_Large_BS),"sd:",sd(energy_Large_BS),"\n");
cat("Middle av:", mean(energy_Middle_BS),"sd:",sd(energy_Middle_BS),"\n");
cat("Small av:", mean(energy_Small_BS),"sd:",sd(energy_Small_BS),"\n");

rm(list=ls());gc();
# Code is over.