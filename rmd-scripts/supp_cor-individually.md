correlation\_supp
================
Filip Wierzbicki
11/16/2022

This scripts contains the pipeline for the supp correlation figure
(inidiviual assemblies) of cluster and non-cluster insertions. Please,
see main-histo.Rmd for preprocessing of required input files.

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)
library(ggpubr)

#######
###For population frequency Info based on Kofler et al. 2015 PLOS Genetics
info1<-read.table("/Users/filipwierzbicki/Desktop/evolution_cluster/temp/TEfamInfo_correct")
names(info1)<-c("name","TE","order","AF","popins")

###exclude somatically regulated TEs based on Malone et al. 2009 Cell
info1<-subset(info1,name!="gypsy10"&name!="gypsy"&name!="ZAM"&name!="gtwin"&name!="gypsy5"&name!="Tabor")

info<-subset(info1,select=c("TE","AF"))
info$AF<-round(info$AF,digits = 1)

infoP<-subset(info1,select=c("name","AF"))
infoP$AF<-round(infoP$AF,digits = 1)
names(infoP)<-c("TE","AF")




####
#cusco and TAS:

t1<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Canton-S_gapped_cusco_tas_summary.forR")
t2<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/DGRP-732_gapped_cusco_tas_summary.forR")
t3<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Iso1_gapped_cusco_tas_summary.forR")
t4<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Oregon-R_gapped_cusco_tas_summary.forR")
t5<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Pi2_gapped_cusco_tas_summary.forR")

tq<-rbind(t1,t2,t3,t4,t5)

names(tq)<-c("count","id","TE","region")

t<-tq
t$id2<-paste(t$id,t$TE,sep="+")

cl<-subset(t,region=="cluster")
cl<-subset(cl,select=c("count","id2"))

rest<-subset(t,region!="cluster") #merges non-cluster and ref regions
rest<-subset(rest,select=c("count","id2"))

for (sid in unique(rest$id2)) { 
  i <- rest$id2 == sid
  a = sum(rest$count[i])
  rest$sum[i] = a
}
rest<-subset(rest,select=c("sum","id2"))
names(rest)<-c("count","id2")
rest<-unique(rest)

cr<-full_join(cl,rest,by="id2")
names(cr)<-c("cluster","id2","noncluster")

cr$id<-gsub("\\+.*","",cr$id2)
cr$TE<-gsub(".*\\+","",cr$id2)

t<-cr


t[is.na(t)] <- 0



t<-left_join(t,info,by="TE")

#including AF threshold
t<-subset(t,AF!="NA")##remove missing AFs
t<-subset(t,AF<=0.2)



cC<-t

cC$logcluster<-log10(cC$cluster+1)
cC$lognoncluster<-log10(cC$noncluster+1)


gC<-ggplot(cC,aes(x=lognoncluster,y=logcluster))+geom_point()+stat_cor(method = "kendall", label.x = 0, label.y = 2.5,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-cluster insertions")+ylab("cluster insertions")+scale_x_continuous(breaks=c(0,1,2,3),labels=c("0","9","99","999"))+scale_y_continuous(breaks=c(0,1,2,3),labels=c("0","9","99","999"))+facet_wrap(~id,nrow=1)


plot(gC)
```

![](supp_cor-individually_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/correlation_supp.pdf",width=8,height=3)
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/correlation_supp.png",width=8,height=3)
```
