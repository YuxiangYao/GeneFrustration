# = = = = = = = = = = = = = = = = = = = = = = = #
# Auxiliary function: ggplot self-defined theme #
# = = = = = = = = = = = = = = = = = = = = = = = #

# color1: [two phenotypic score.]"#EE312E","#2156A6","#FFB600", "#3B2E7E"
# color2: [four cluster] '#909090','#d64e3d','#4087af','#eec67f','#9a91c8'
# color3: [-1,+1] low="#8359B5", high="#008844",mid = "white",

library(ggplot2);
Plot_ThemeConfiguration<-function(GGplotObject){
  # Only adjust theme:
  pics=GGplotObject;
  pics=pics+theme(
    axis.line=element_blank(),
    axis.ticks.x.top=element_line(linewidth=0.75,colour="#000000"),
    axis.ticks.x.bottom=element_line(linewidth=0.75,colour="#000000"),
    axis.ticks.y.left=element_line(linewidth=0.75,colour="#000000"),
    axis.ticks.y.right=element_line(linewidth=0.75,colour="#000000"),
    panel.grid=element_blank(), # Delete grid
    panel.background=element_blank(), # Delete bakcground
    panel.border = element_rect(color = "#000000",fill = NA,linewidth = 1.5),
    legend.key=element_rect(fill="#FFFFFF"),
    legend.title=element_text(face="bold",size=14),
    legend.text=element_text(face="bold",size=12),
    title=element_text(size=20),
    legend.background = element_rect(fill = "transparent"),
    legend.position=c(.70,.50),
    plot.title=element_text(hjust=0.5),
    axis.text.x=element_text(size=12,colour="#000000",face="bold",vjust=0.5,hjust=0.5,angle=0),
    axis.text.y=element_text(size=12,colour="#000000",face="bold",hjust=1.0,vjust=0.5,angle=0),
    axis.title.x.bottom=element_text(size=15,face="bold",colour="black",angle=0),
    axis.title.y.left=element_text(size=15,face="bold",colour="black"));
  return(pics);
}
