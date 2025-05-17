# Execute the pressure vs Frustration vs Potential values.
# Related to Local Frustration series, Cracking PCA series.
# Due to 10^6 cases and recording many detail, we calculate them in parallel (See end of this file)
args <- commandArgs(trailingOnly=TRUE);
param1 <- as.integer(as.numeric(args[1]));# Random Seed
param2 <- as.character(args[2]);# Gene Names
param3 <- as.character(args[3])# File Names 

# Load data and complie function
Rcpp::sourceCpp("../cpp_Source/c_OEKD_StateTransition_onlyPou5f1.cpp");
adjmat=readRDS("../Data_Used/AdjacentMatrix.Matrix.rds");
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


source("../Simulation/exe_LocalFrustrations_AllGene_OEKD_SetSM.R");
MEFs=SimulatedCharacteristic(adjmat,param1);# <-This generate 10^5 MEF-like cells.
PartialResult=Pressure2Potential(adjmat, MEFs, Marker_List,
  match(param2,colnames(adjmat)), c(1), param1+1);
# Note that you can change the temporary directory. 
saveRDS(PartialResult, file=paste0("/dev/shm/part_",param2,"_",param3,".rds") ,compress=FALSE);
saveRDS(MEFs, file=paste0("/dev/shm/part_",param2,"_",param3,"_MEFs.rds") ,compress=FALSE);
# Code is over.

# Parallel analysis:

# [code.order]:
# Rscript	exe_LocalFrustrations_AllGene_OEKD_SetSM_Parallel.R 12940 Sox2 n0
# Rscript	exe_LocalFrustrations_AllGene_OEKD_SetSM_Parallel.R	12941 Sox2 n1
# Rscript	exe_LocalFrustrations_AllGene_OEKD_SetSM_Parallel.R	12942 Sox2 n2
# Rscript	exe_LocalFrustrations_AllGene_OEKD_SetSM_Parallel.R	12943 Sox2 n3
# Rscript	exe_LocalFrustrations_AllGene_OEKD_SetSM_Parallel.R	12944 Sox2 n4
# Rscript	exe_LocalFrustrations_AllGene_OEKD_SetSM_Parallel.R	12945 Sox2 n5
# Rscript	exe_LocalFrustrations_AllGene_OEKD_SetSM_Parallel.R	12946 Sox2 n6
# Rscript	exe_LocalFrustrations_AllGene_OEKD_SetSM_Parallel.R	12947 Sox2 n7
# Rscript	exe_LocalFrustrations_AllGene_OEKD_SetSM_Parallel.R	12948 Sox2 n8
# Rscript	exe_LocalFrustrations_AllGene_OEKD_SetSM_Parallel.R	12949 Sox2 n9

# Do parallel (in Linux)
# parallel --jobs 10 --joblog OE_Pou5f1.log < code.order &

# One also can only use single threaded analysis by: 
# AllResult=Pressure2Potential(adjmat, MEFs, Marker_List,
#   match(c('Pou5f1'),colnames(adjmat)), c(1), 1234567L);