
library(ggplot2)
t<-read.table("/Volumes/Temp3/filip/trap_model/clusterscore/fracqual/output/summary/gapless_th0")
names(t)<-c("cl","cov","sc","id")
t$assembly<-gsub("_.*","",t$id)
t$condition<-gsub(".*_ml","",t$id)
t$condition<-gsub("_th0","",t$condition)


c<-subset(t,select=c("cl","cov","assembly","condition"))
names(c)<-c("cl","value","assembly","condition")
c$type<-c("coverage")
s<-subset(t,select=c("cl","sc","assembly","condition"))
names(s)<-c("cl","value","assembly","condition")
s$type<-c("softclip")
t<-rbind(c,s)
t$cl<-as.factor(t$cl)

t<-subset(t,type=="coverage")

g<-ggplot(t,aes(x=cl,y=value))+geom_bar(stat="identity",position="dodge",color="black")+ylab("fraction of outlier signals")+xlab("piRNA cluster ID")+ylim(0.00,1.00)+theme_bw()+facet_wrap(~assembly,ncol = 5)+theme(text=element_text(size=8),axis.text.x = element_text(angle = 45,hjust=1,size=5),legend.title = element_blank(),legend.position = "None")
plot(g)
#,alpha=type


