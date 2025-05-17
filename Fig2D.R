# Figure: Energy landscape of phenotypic phase (PCA landscape).

# Load data & function: (Data from previous processing)
index.nonrep=readRDS("./Intermediate_Result/nonduplicate_index.rds");
Steps=readRDS("./Intermediate_Result/Step_Energy.rds");
pca_index=readRDS("./Intermediate_Result/PCA_index_coor.rds");
pca_index=cbind.data.frame(pca_index[,2:3]);

# Bin the location at X/Y axis. 
bin=40;# Here we use the bin of 40.
x_edge=c(-2.5,4.0)#c(min(pca_index[,1]),max(pca_index[,1]));# Set -2.5 -> 4.0
y_edge=c(-2.5,3.0)#c(min(pca_index[,2]),max(pca_index[,2]));# Set -2.5 -> 3.0
seq_x=seq(x_edge[1],x_edge[2],length.out=(bin+1));
seq_y=seq(y_edge[1],y_edge[2],length.out=(bin+1));
bin_x=floor((pca_index[,1]+2.5)/(seq_x[2]-seq_x[1]));# Gongcha-X
bin_y=floor((pca_index[,2]+2.5)/(seq_y[2]-seq_y[1]));# Gongcha-Y
energys=Steps[index.nonrep,2];# [-233, +29], in fact, but interval shrinking via binned avg.

# cpp-based quick count average energy in each bin.
Rcpp::sourceCpp("./cpp_Source/c_DoubleBinCountAv.cpp");
count_bin=DoubleBins(bin_x,bin_y,energys,40L,40L);
colnames(count_bin)=rownames(count_bin)=paste0("x_",c(1:40));

# Prepare for dataframe of ggplot-geom_tile.
library(ggplot2);
source("./Auxiliary_Functions_for_Figure/aux_ggplot2Theme.R");
count_bin=reshape2::melt(count_bin);
count_bin$value[count_bin$value>-1]=NA;
figs=ggplot(count_bin)+geom_tile(aes(x=Var1,y=Var2,fill=value),)+
  scale_fill_viridis_c(limits=c(-225,-35),option = "magma",
    name="Pseudo-energy", na.value = "#ffffff")+
  scale_x_discrete(expand=c(0,0), breaks=paste0("x_",c(1,9,17,25,33,40)),
    labels=as.character(c(-2.5,-1.2,-0.1,1.4,2.7,4.0)))+
  scale_y_discrete(expand=c(0,0), breaks=paste0("x_",c(1,9,17,25,33,40)),
    labels=as.character(c(-2.5,-1.4,-0.3,0.8,1.9,3.0)))+
  labs(x="PC2", y="PC3");
figs=Plot_ThemeConfiguration(figs)+
  theme(aspect.ratio=1, legend.position="right");
ggsave(filename="./Figure_OriginalVersion/BinCount_Energy.pdf",
  figs, width=7.5, height=4.5, units="in");

rm(list=ls());gc();
# Code is over.