# Combine each individual file as whole file/object.
part01=readRDS("/dev/shm/part_A1.rds");# Note your correct temporary address.
part02=readRDS("/dev/shm/part_A2.rds");
part03=readRDS("/dev/shm/part_A3.rds");
part04=readRDS("/dev/shm/part_A4.rds");
part05=readRDS("/dev/shm/part_A5.rds");
part06=readRDS("/dev/shm/part_A6.rds");
part07=readRDS("/dev/shm/part_A7.rds");
part08=readRDS("/dev/shm/part_A8.rds");
part09=readRDS("/dev/shm/part_A9.rds");
part10=readRDS("/dev/shm/part_Aa.rds");
all.data=list(part01,part02,part03,part04,part05,part06,part07,part08,part09,part10);
rm(part01,part02,part03,part04,part05,part06,part07,part08,part09,part10);

# Combined all data.
combined.vector<-function(aList,Whichid){
  tmp=aList[[1]][[Whichid]];
  lens=length(tmp);
  nsam=length(aList);
  returner=rep(NA,lens*nsam);
  left.id=1;right.id=lens;
  for(ii in c(1:length(aList))){
    returner[left.id:right.id]=aList[[ii]][[Whichid]];
    left.id=left.id+lens;
    right.id=right.id+lens;}
  return(returner);
}
combined.matrix<-function(aList,Whichid){
  tmp=aList[[1]][[Whichid]];
  nsam=length(aList);
  nlen=nrow(tmp);
  nval=ncol(tmp);
  returner=matrix(NA,nlen*nsam,nval);
  left.id=1;right.id=nlen;
  for(ii in c(1:length(aList))){
    returner[left.id:right.id,]=aList[[ii]][[Whichid]];
    left.id=left.id+nlen;
    right.id=right.id+nlen;}
  return(returner);
}
combined.list<-function(aList,Whichid){
  tmp=aList[[1]][[Whichid]];
  nsam=length(aList);
  nlen=length(tmp);
  returner=as.list(rep(NA,nlen*nsam));
  left.id=1;right.id=nlen;
  for(ii in c(1:length(aList))){
    returner[left.id:right.id]=aList[[ii]][[Whichid]];
    left.id=left.id+nlen;
    right.id=right.id+nlen;}
  return(returner);
}

AllResult=all.data[[1]]; #as.slot
{
  AllResult$Ori_PS=combined.vector(all.data,1);
  AllResult$Ori_EM=combined.vector(all.data,2);
  AllResult$OriEng=combined.vector(all.data,3);
  AllResult$GeneDiff=combined.vector(all.data,4);
  AllResult$NewEng=combined.vector(all.data,5);
  AllResult$New_PS=combined.vector(all.data,6);
  AllResult$New_EM=combined.vector(all.data,7);
  AllResult$Pressure=combined.vector(all.data,8);
  AllResult$NewTime=combined.vector(all.data,9);
  AllResult$EachUpdate=combined.matrix(all.data,10);
  AllResult$EachTempEng=combined.matrix(all.data,11);
  AllResult$NewState=combined.matrix(all.data,12);
  AllResult$TransisMat=combined.list(all.data,13);
}
rm(combined.vector,combined.matrix,combined.list,all.data);
# The file size is too large. You can choose whether to compress it as required. 
# For compression, the "fastSave" library is recommended.
saveRDS(AllResult, file="../Intermediate_Result/All_Result_onlyPou5f1.rds", compress=FALSE);

rm(list=ls());gc();
# Code is over.