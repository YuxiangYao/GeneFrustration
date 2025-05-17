# Execute 1E+7 random initial cells into stable states.

# Compile simulation function: B_ScanningAttrctors_Fast()
Rcpp::sourceCpp("./cpp_Source/c_GotoStableDistribution.cpp");

# Load genetic regulatory network.
adjmat=readRDS("./Data_Used/AdjacentMatrix.Matrix.rds");
# Simulation (mc.attr: random Monte Carlo [[1]]Stable State; [[2]] step/energy)
mc.attr=B_ScanningAttrctors_Fast(adjmat,10000000L,343780L);
# Extract stable states of all 10^7 initials
sss=mc.attr$StateMatrix;
# Remove the duplicated cases.
index.nonrep=!duplicated(sss);
sss=sss[index.nonrep,];
# Highly recommend! for subsequent using.
saveRDS(mc.attr, file="./Intermediate_Result/mc_attr_origin.rds", compress=FALSE);
saveRDS(index.nonrep, file="./Intermediate_Result/nonduplicate_index.rds", compress=FALSE);
saveRDS(mc.attr$Steps, file="./Intermediate_Result/Step_Energy.rds", compress=FALSE);
saveRDS(sss, file="./Intermediate_Result/SSS_origin_nodup.rds", compress=FALSE);
# Code is over.