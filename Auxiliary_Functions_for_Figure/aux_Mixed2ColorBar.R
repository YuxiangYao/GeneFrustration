# = = = = = = = = = = = = = = = = = = = = #
# Auxiliary function: double-mixed colbar #
# = = = = = = = = = = = = = = = = = = = = #
ColMix<-function(ColVec,FracVec){# ColVals: four colour,
  col01=as.vector(col2rgb(ColVec[3])); col11=as.vector(col2rgb(ColVec[4]));
  col00=as.vector(col2rgb(ColVec[1])); col10=as.vector(col2rgb(ColVec[2]));
  x=FracVec[1];y=FracVec[2];
  MixCol=(1-x)*(1-y)*col00+x*(1-y)*col10+
    (1-x)*y*col01+x*y*col11;
  MixCol=rgb(MixCol[1],MixCol[2],MixCol[3],maxColorValue=255);
  return(MixCol);
}
ColMix_Vec<-function(ColVec,FracVec){# ColVals: four colour,
  col01=as.vector(col2rgb(ColVec[3])); col11=as.vector(col2rgb(ColVec[4]));
  col00=as.vector(col2rgb(ColVec[1])); col10=as.vector(col2rgb(ColVec[2]));
  x=FracVec[,1];y=FracVec[,2];
  MixCol_r=(1-x)*(1-y)*col00[1]+x*(1-y)*col10[1]+(1-x)*y*col01[1]+x*y*col11[1];
  MixCol_g=(1-x)*(1-y)*col00[2]+x*(1-y)*col10[2]+(1-x)*y*col01[2]+x*y*col11[2];
  MixCol_b=(1-x)*(1-y)*col00[3]+x*(1-y)*col10[3]+(1-x)*y*col01[3]+x*y*col11[3];
  MixCol=rgb(MixCol_r,MixCol_g,MixCol_b,maxColorValue=255);
  return(MixCol);
}
TestColMix<-function(ColVec=c("red","green","blue","white"),
  xyScale=c(10,10),Ratios=1.0){
  x.s=seq(0,1,length.out=xyScale[1]);
  y.s=seq(0,1,length.out=xyScale[2]);
  tmp=cbind(
    x=rep(x.s,each=xyScale[2]),
    y=rep(y.s,time=xyScale[1]))
  colss=rep(NA,nrow(tmp));
  for(ii in c(1:nrow(tmp))){
    colss[ii]=ColMix(ColVec,tmp[ii,]);}
  tmp=cbind.data.frame(tmp,colss);
  colnames(tmp)=c("x","y","c");
  colorbars=unique(tmp$c);
  tmp$c=factor(tmp$c,levels = colorbars);
  fig=ggplot()+
    geom_tile(data=tmp,aes(x=x,y=y,fill=c))+
    scale_fill_manual(values=colorbars)+
    theme(legend.position = "none",aspect.ratio =Ratios)
  return(fig);
}

#TestColMix(c('#cdcdcd',"#2156A6","#EE312E","#eedd33"),c(25,25))
#TestColMix(c('#000000',"#ff0000","#00ff00","#0000ff"),c(25,25))
