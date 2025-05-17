# Execute the pressure vs Frustration vs Potential values.
# Related to Local Frustration series, Cracking PCA series.
# Due to 10^6 cases and recording many detail, we calculate them in parallel (See end of this file)
args <- commandArgs(trailingOnly=TRUE);
param1 <- as.integer(as.numeric(args[1]));# Line_Start
param2 <- as.integer(as.numeric(args[2]));# Line_End
param3 <- as.integer(as.numeric(args[3]));# Random Seed
param4 <-as.character(args[4]);# File Names

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

extreme=readRDS("../Intermediate_Result/ExtremeCase_FourClusters.rds");
MEFs=extreme$sm;
PartialResult=Pressure2Potential(adjmat, MEFs[param1:param2,], Marker_List,
  match(c('Pou5f1'),colnames(adjmat)), c(1), param3);
# Note that you can change the temporary directory. 
saveRDS(PartialResult,file =paste0("/dev/shm/part_",param4,".rds") ,compress=FALSE);

# Code is over.

# Parallel analysis:

# [code.order]:
# Rscript	exe_OverExpressSinglePou5f1_ClassicalCases_Parallel.R	1	100000	241113	A1
# Rscript	exe_OverExpressSinglePou5f1_ClassicalCases_Parallel.R	100001	200000	241114	A2
# Rscript	exe_OverExpressSinglePou5f1_ClassicalCases_Parallel.R	200001	300000	241115	A3
# Rscript	exe_OverExpressSinglePou5f1_ClassicalCases_Parallel.R	300001	400000	241116	A4
# Rscript	exe_OverExpressSinglePou5f1_ClassicalCases_Parallel.R	400001	500000	241117	A5
# Rscript	exe_OverExpressSinglePou5f1_ClassicalCases_Parallel.R	500001	600000	241118	A6
# Rscript	exe_OverExpressSinglePou5f1_ClassicalCases_Parallel.R	600001	700000	241119	A7
# Rscript	exe_OverExpressSinglePou5f1_ClassicalCases_Parallel.R	700001	800000	241120	A8
# Rscript	exe_OverExpressSinglePou5f1_ClassicalCases_Parallel.R	800001	900000	241121	A9
# Rscript	exe_OverExpressSinglePou5f1_ClassicalCases_Parallel.R	900001	1000000	241122	Aa

# Do parallel (in Linux)
# parallel --jobs 10 --joblog OE_Pou5f1.log < code.order &

# One also can only use single threaded analysis by: 
# AllResult=Pressure2Potential(adjmat, MEFs, Marker_List,
#   match(c('Pou5f1'),colnames(adjmat)), c(1), 1234567L);