# Figure: Obtain the basin size of 10^7 step. We can obtian a log-log curve.

# Load data & necessary function:
source("./Auxiliary_Functions_for_Figure/aux_ggplot2Theme.R");
mc.attr=readRDS("./Intermediate_Result/mc_attr_origin.rds");
# Get stable states (including duplicate cases).
sss=mc.attr$StateMatrix;
# Turn to character vectors. (To compare and count the basin)
sss=apply(sss,1,paste,collapse="");
# Each isoform stable state has how much attracting area.
sss=as.matrix(table(sss));
# Counter the area-size distribution.
sss=as.matrix(table(sss));

# Prepare for data.frame of log-log of basin size and fitting lines.
part1=cbind.data.frame(log2(as.numeric(rownames(sss))),log2(sss),"single");
colnames(part1)=c('x','y','c');
coeffs=lm(part1$y ~ part1$x)$coefficients;# Instead of bin-log, here we use fitting line to guide.
part2=cbind.data.frame(seq(1,6,by=0.1),seq(1,6,by=0.1)*coeffs[2]+coeffs[1]+1.0,"fitting");
colnames(part2)=colnames(part1);
cat(coeffs[2],"\n");# This is the slope.

# Figure for distributions of attractor basins.
library(ggplot2);
fig_attr=ggplot()+
  geom_point(data=part1,aes(x=x,y=y,col=c,shape=c),size=2.2,alpha=0.5)+
  scale_color_manual(values=c("#878787","#cfcfcf"))+
  scale_shape_manual(values = c(16,15))+
  geom_line(data=part2,aes(x=x,y=y),col="#222222",linewidth=1.5)+
  scale_x_continuous(limits=c(-0.5,7.3),breaks=c(0:7), 
    labels=c(latex2exp::TeX("${2^0}$",bold=TRUE),latex2exp::TeX("${2^1}$"),
      latex2exp::TeX("${2^2}$"),latex2exp::TeX("${2^3}$"),latex2exp::TeX("${2^4}$"),
      latex2exp::TeX("${2^5}$"),latex2exp::TeX("${2^6}$"),latex2exp::TeX("${2^7}$")))+
  scale_y_continuous(limits=c(-0.5,24.5),breaks=c(0,4,8,12,16,20,24), 
    labels=c(latex2exp::TeX("${2^0}$"),latex2exp::TeX("$2^4$"),
      latex2exp::TeX("$2^8$"),latex2exp::TeX("$2^{12}$"),
      latex2exp::TeX("$2^{16}$"),latex2exp::TeX("$2^{20}$"),latex2exp::TeX("$2^{24}$")));
fig_attr=Plot_ThemeConfiguration(fig_attr)+theme(aspect.ratio=0.60,legend.position=c(.75,.75));
ggsave(filename="./Figure_OriginalVersion/BasinSize_LogLogCurve.pdf",fig_attr,width=4.0,height=4.0,units="in");

rm(list=ls());gc();
# Code is over. 