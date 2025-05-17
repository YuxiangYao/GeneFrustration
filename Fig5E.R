
# Functions: Filter Cells.
FilterSuitbaleCell<-function(Cell.Matrix){
  cellmat=Cell.Matrix;
  # Filter some invalid samples and genes
  sam_c=colSums(cellmat);# Counter number.
  sam_g=colSums(cellmat>0);# Non-zero gene's number.
  # Total-counter of samples fitted as normal distribution.
  mean_sd=MASS::fitdistr(sam_c,"normal")$estimate;
  # Get up-down boundary of valid sample.
  upper_count=qnorm(0.95,mean_sd[1],mean_sd[2]);
  down_count =qnorm(0.05,mean_sd[1],mean_sd[2]);
  # Total-counter of samples fitted as normal distribution.
  mean_sd=MASS::fitdistr(sam_g,"normal")$estimate;
  # Get up-down boundary of valid sample.
  upper_genenum=qnorm(0.95,mean_sd[1],mean_sd[2]);
  down_genenum =qnorm(0.05,mean_sd[1],mean_sd[2]);
  # Select appropriate counting & gene-expression samples.
  cellmat=cellmat[,sam_c<upper_count&sam_c>down_count&sam_g<upper_genenum&sam_g>down_genenum];
  cellmat=as.matrix(cellmat);
  return(cellmat);
}
# Load data: GSE103221(OSK-10X)
# Please download from [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103221] and extract them to your directory
osk=list(
  mef=as.data.frame(data.table::fread("./your/directory/GSM3629847_10x_osk_mef.csv",sep=",")),
  d00=as.data.frame(data.table::fread("./your/directory/GSM3629850_10x_osk_d0.csv",sep=",")),
  d01=as.data.frame(data.table::fread("./your/directory/GSM3629849_10x_osk_d1.csv",sep=",")),
  d03=as.data.frame(data.table::fread("./your/directory/GSM3629851_10x_osk_d3.csv",sep=",")),
  d05=as.data.frame(data.table::fread("./your/directory/GSM3629852_10x_osk_d5.csv",sep=",")),
  d06=as.data.frame(data.table::fread("./your/directory/GSM3629853_10x_osk_d6.csv",sep=",")),
  d07=as.data.frame(data.table::fread("./your/directory/GSM3629854_10x_osk_d7.csv",sep=",")),
  d08=as.data.frame(data.table::fread("./your/directory/GSM3629855_10x_osk_d8.csv",sep=",")),
  esc=as.data.frame(data.table::fread("./your/directory/GSM3629848_10x_osk_esc.csv",sep=",")))
cell_name=osk[[1]][,1];
for(ii in c(1:9)){
  osk[[ii]]=osk[[ii]][,-1];
  rownames(osk[[ii]])=cell_name;
  colnames(osk[[ii]])=paste0(names(osk)[ii],"_",sprintf("%04d",c(1:ncol(osk[[ii]]))))}
osk_mat=osk[[1]];
for(ii in c(2:9)){
  osk_mat=cbind.data.frame(osk_mat,osk[[ii]]);}

# Filter valid cells.
osk_mat=FilterSuitbaleCell(osk_mat);
osk_mat=t(osk_mat);
# Set time-point labels.
mef=grep("^mef",rownames(osk_mat)); d00=grep("^d00",rownames(osk_mat));
d01=grep("^d01",rownames(osk_mat)); d03=grep("^d03",rownames(osk_mat));
d05=grep("^d05",rownames(osk_mat)); d06=grep("^d06",rownames(osk_mat));
d07=grep("^d07",rownames(osk_mat)); d08=grep("^d08",rownames(osk_mat));
esc=grep("^esc",rownames(osk_mat));

# Function: return instant candidate number.
InstantTr<-function(Expmat,AdjMat_Ori){#Row: samples; Col: genes
  index=intersect(colnames(AdjMat_Ori),colnames(Expmat));
  adjm=AdjMat_Ori[index,index];
  expmat=Expmat[,index];
  inducer=rep(0,ncol(expmat));
  names(inducer)=colnames(expmat);
  inducer[c("Pou5f1","Sox2","Klf4","Myc")]=1;# Yamanaka factors. 
  # Although not shown, in fact, the exogenous effects of these factors persist
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

# Load data & obtain instant Tr(x,ε) in each time point.
adjmat=readRDS("./Data_Used/AdjacentMatrix.Matrix.rds");
source("./Auxiliary_Functions_for_Figure/aux_ggplot2Theme.R");
inst_tr_xe=InstantTr(osk_mat,adjmat);

set.seed(103221);
# Set sub.dataframe for integrated figure.
Histo_Trigger2FitGaussion<-function(a_vec,celltype="word"){
  FittedData=cbind.data.frame(a_vec,celltype);
  colnames(FittedData)=c("x","t")
  return(FittedData);
}
# Prepare the dataframe.
datas=rbind.data.frame(Histo_Trigger2FitGaussion(inst_tr_xe[mef,1],"aMEF"),
  Histo_Trigger2FitGaussion(inst_tr_xe[d00,1],"Day0"),Histo_Trigger2FitGaussion(inst_tr_xe[d01,1],"Day1"),
  Histo_Trigger2FitGaussion(inst_tr_xe[d03,1],"Day3"),Histo_Trigger2FitGaussion(inst_tr_xe[d05,1],"Day5"),
  Histo_Trigger2FitGaussion(inst_tr_xe[d06,1],"Day6"),Histo_Trigger2FitGaussion(inst_tr_xe[d07,1],"Day7"),
  Histo_Trigger2FitGaussion(inst_tr_xe[d08,1],"Day8") );
datas$t<-factor(datas$t, levels=c(paste0("Day",c(8:5,3,1,0)),"aMEF"));# Set label order!

t.test(datas[datas$t=="aMEF",1],datas[datas$t=="Day0",1],alternative = "greater")$p.value
t.test(datas[datas$t=="Day0",1],datas[datas$t=="Day1",1],alternative = "greater")$p.value
t.test(datas[datas$t=="Day1",1],datas[datas$t=="Day3",1],alternative = "greater")$p.value
t.test(datas[datas$t=="Day3",1],datas[datas$t=="Day5",1],alternative = "greater")$p.value
t.test(datas[datas$t=="Day5",1],datas[datas$t=="Day6",1],alternative = "greater")$p.value
t.test(datas[datas$t=="Day6",1],datas[datas$t=="Day7",1],alternative = "greater")$p.value
t.test(datas[datas$t=="Day7",1],datas[datas$t=="Day8",1],alternative = "greater")$p.value

# Combinded figure for Possion curves.

FitCurves<-function(a_vec,chrx){
  tmp=cbind.data.frame(round(mean(a_vec),2),round(sd(a_vec),2),chrx);
  colnames(tmp)=c("av","sd","name");
  return(tmp);}

datas_scatter=rbind.data.frame(FitCurves(inst_tr_xe[mef,1],"aMEF"),
  FitCurves(inst_tr_xe[d00,1],"Day0"),FitCurves(inst_tr_xe[d01,1],"Day1"),
  FitCurves(inst_tr_xe[d03,1],"Day3"),FitCurves(inst_tr_xe[d05,1],"Day5"),
  FitCurves(inst_tr_xe[d06,1],"Day6"),FitCurves(inst_tr_xe[d07,1],"Day7"),
  FitCurves(inst_tr_xe[d08,1],"Day8") );
datas_scatter$name<-factor(datas_scatter$name, levels=c(paste0("Day",c(8:5,3,1,0)),"aMEF"));# Set label order!

# Fig. 5(g): Boxplot for each time point's Tr(x,ε).
figs=ggplot()+
  geom_hline(yintercept = c(21,24,27,30),linewidth=0.50,color="#aaaaaa")+
  geom_hline(yintercept = c(22.5,25.5,28.5),linewidth=0.25,color="#cfcfcf")+
  geom_pointrange(data=datas_scatter,aes(x = name, y = av, ymin = av-sd, ymax =av+sd,col=name),linewidth=1.5) +
  scale_color_manual(values=c('#343434','#d53e4f','#f46d43',
    '#fdae61','#abd3a4','#66c2a5','#3288bd','#5e4fa2')[8:1])+
  guides(color=guide_legend(nrow=1));#+coord_flip();
figs=Plot_ThemeConfiguration(figs)+
  theme(aspect.ratio=0.333,axis.line=element_blank(),legend.position="top",);
ggsave(filename="./Figure_OriginalVersion/RNAseq_sc_boxplot.pdf",figs,width=7.0,height=3.0,units="in");

rm(list=ls());gc();
# Code is over.