supp\_butterfly-cor
================
Filip Wierzbicki
10/10/2022

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

dsl<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/butterfly/signatures/output-V3/TE/filtered/forR/TE_w500.forR")
names(dsl)<-c("TE","chr","start","end","ls","la","rs","ra","strain")

a1<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Canton-S_gapped_cusco_tas_summary.forR")
a2<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/DGRP-732_gapped_cusco_tas_summary.forR")
a3<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Iso1_gapped_cusco_tas_summary.forR")
a4<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Oregon-R_gapped_cusco_tas_summary.forR")
a5<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Pi2_gapped_cusco_tas_summary.forR")
A<-rbind(a1,a2,a3,a4,a5)
names(A)<-c("count","strain","TE","region")
A$id<-paste(A$strain,A$TE,sep="+")

###
ic<-dsl
ic$id<-paste(ic$strain,ic$TE,sep="+")
for (sid in unique(ic$id)) { 
  i <- ic$id == sid
  x = nrow(subset(ic,ic$id==sid))
  ic$butterfly[i] = x
}
bfs<-subset(ic,select = c("butterfly","id"))
bfs<-unique(bfs)
ac<-subset(A,region=="cluster")
ac<-subset(ac,select=c("count","id"))
TE<-full_join(bfs,ac,by="id")
anc<-subset(A,region!="cluster")
anc<-subset(anc,select=c("count","id"))

for (sid in unique(anc$id)) { 
  i <- anc$id == sid
  a = sum(anc$count[i])
  anc$sum[i] = a
}
anc<-subset(anc,select=c("sum","id"))
anc<-unique(anc)

names(anc)<-c("rest","id")
TE<-full_join(TE,anc,by="id")

TE$strain<-gsub("\\+.*","",TE$id)
TE$TE<-gsub(".*\\+","",TE$id)

t<-TE

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


t[is.na(t)] <- 0

for (sid in unique(t$TE)) { 
  i <- t$TE == sid
  a = mean(t$count[i])
  b = mean(t$butterfly[i])
  t$avrcl[i] = a
  t$avrbut[i] = b
}

t<-subset(t,select = c("TE","avrcl","avrbut"))
t<-unique(t)

t<-left_join(t,info,by="TE")

#including AF threshold
t<-subset(t,AF!="NA")##remove missing AFs
t<-subset(t,AF<=0.2)
t$clusterlog<-log10(t$avrcl+1)
t$butterflylog<-log10(t$avrbut+1)


gcor<-ggplot(t,aes(x=clusterlog,y=butterflylog))+geom_point()+stat_cor(cor.coef.name="tau", method = "kendall", label.x = 0.25, label.y = 1.35,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("cluster insertions")+ylab("DSL insertions")+scale_x_continuous(breaks=c(0,1.041393,2.004321,3.000434),labels=c("0","10","100","1000"))+scale_y_continuous(breaks=c(0,1.041393,2.004321,3.000434),labels=c("0","10","100","1000"))

plot(gcor)
```

![](supp_butterfly-cor_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/supp_butterfly-cor.pdf",width=7,height=6)
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/supp_butterfly-cor.png",width=7,height=6)


ar<-subset(A,region!="cluster")

ar$id<-paste(ar$strain,ar$TE,sep="+")
for (sid in unique(ar$id)) { 
  i <- ar$id == sid
  a = sum(ar$count[i])
  ar$sum[i] = a
}
ar<-subset(ar,select=c("sum","id"))
ar<-unique(ar)

ta<-full_join(ar,TE,by="id")
ta[is.na(ta)] <- 0
ta<-left_join(ta,info,by="TE")

#including AF threshold
ta<-subset(ta,AF!="NA")##remove missing AFs
ta<-subset(ta,AF<=0.2)



ta0<-subset(ta,count==0)

hier<-read.table("/Users/filipwierzbicki/Desktop/trap_model/data/other/new_dmel_132cons_hier",header =TRUE)
hier<-subset(hier,select=c("id","family"))
names(hier)<-c("TE","family")
ta0<-left_join(ta0,hier,by="TE")
ta0<-subset(ta0,select=c("family","sum","butterfly","count","strain"))
names(ta0)<-c("family","non-cluster","other-SL","cluster","strain")

write.table(ta0,"/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/butterfly/output/missing-clusterinsertions.txt",quote=FALSE)
```
