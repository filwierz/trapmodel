correlation\_main
================
Filip Wierzbicki
6/22/2022

This scripts contains the pipeline for the main correlation figure of
cluster and non-cluster insertions. Please, see main-histo.Rmd for
preprocessing of required input files.

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
library(ggtext)
####
#simulation
output<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/constant-u/run3-seed/combined/output-constant_u")


output<-output[,-28]
names(output)<-c("replicate","generation","delim1",
                 "fwt", "w",    "tes",  "popfreq",  "fixed","delim2",
                 "fwcli",   "cluins",   "cluins_popfreq","cluins_fixed",    "phase","delim3",
                 "fwrefi",  "refins",   "refins_popfreq", "refins_fixed","delim4",
                 "novel",   "sites",    "clusites", "tes_stdev" ,"cluins_stdev" ,"fw0", "w_min","popsize")



ts<-subset(output,generation==2000)

ts$tesc<-ts$tes-ts$cluins
#ts$tesr<-ts$tes-ts$refins#-ts$cluins
t<-subset(ts,select = c("tesc","cluins"))
names(t)<-c("global","local")

t$global<-t$global/2 #for halpoid abundance
t$local<-t$local/2 #for halpoid abundance


gS<-ggplot(t, aes(x=global, y=local)) + geom_point()+stat_cor(cor.coef.name="tau",method = "kendall", label.x = 50.1, label.y = 6.5,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-cluster insertions")+ylab("cluster insertions")+xlim(50,200)+ggtitle("expected - simulations")+theme(plot.title = element_markdown(size=11))





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



###popTE2:

#t<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/popTE2/sep.mf30_brenecluster.forR")
#names(t)<-c("count","id","TE","region")

tseg<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/popTE2/haploid_counts/sep.mf30_brenecluster-seg")
tfix<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/popTE2/haploid_counts/sep.mf30_brenecluster-fixed")
t<-rbind(tseg,tfix)
names(t)<-c("count","id","TE","region","id2")
for (sid in unique(t$id2)) { 
  i <- t$id2 == sid
  a = sum(t$count[i])
  t$sum[i] = a
}
t<-subset(t,select=c("sum","id","TE","region"))
t<-unique(t)
names(t)<-c("count","id","TE","region")
t$count<-round(t$count+0.00000001)

t$id2<-paste(t$id,t$TE,sep="_")

cl<-subset(t,region=="cl")
cl<-subset(cl,select=c("count","id2"))

rest<-subset(t,region!="cl")
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

cr$id<-gsub("_.*","",cr$id2)
cr$TE<-gsub(".*_","",cr$id2)
t<-cr


t[is.na(t)] <- 0


for (sid in unique(t$TE)) { 
  i <- t$TE == sid
  a = mean(t$cluster[i])
  b = mean(t$noncluster[i])
  t$avrcl[i] = a
  t$avrrest[i] = b
}

t<-subset(t,select = c("TE","avrcl","avrrest"))
t<-unique(t)



t<-left_join(t,infoP,by="TE")
#including AF threshold
t<-subset(t,AF!="NA")##remove missing AFs
t<-subset(t,AF<=0.2)

cC<-t

cC$cluster<-log10(cC$avrcl+1)
cC$noncluster<-log10(cC$avrrest+1)


gPop<-ggplot(cC,aes(x=noncluster,y=cluster))+geom_point()+stat_cor(cor.coef.name="tau",method = "kendall", label.x = 0, label.y = 1.42,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-cluster insertions")+ylab("cluster insertions")+scale_x_continuous(breaks=c(0,1.041393,2.004321,3.000434),labels=c("0","10","100","1000"))+scale_y_continuous(breaks=c(0,1.041393,2.004321,3.000434),labels=c("0","10","100","1000"))+ggtitle("observed - reference assembly and reference annotations")+theme(plot.title = element_markdown(size=11))







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

for (sid in unique(t$TE)) { 
  i <- t$TE == sid
  a = mean(t$cluster[i])
  b = mean(t$noncluster[i])
  t$avrcl[i] = a
  t$avrrest[i] = b
}

t<-subset(t,select = c("TE","avrcl","avrrest"))
t<-unique(t)

t<-left_join(t,info,by="TE")

#including AF threshold
t<-subset(t,AF!="NA")##remove missing AFs
t<-subset(t,AF<=0.2)



cC<-t

cC$cluster<-log10(cC$avrcl+1)
cC$noncluster<-log10(cC$avrrest+1)


gC<-ggplot(cC,aes(x=noncluster,y=cluster))+geom_point()+stat_cor(cor.coef.name="tau",method = "kendall", label.x = 0.5, label.y = 2.35,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-cluster insertions")+ylab("cluster insertions")+scale_x_continuous(breaks=c(0,1.041393,2.004321,3.000434),labels=c("0","10","100","1000"))+scale_y_continuous(breaks=c(0,1.041393,2.004321,3.000434),labels=c("0","10","100","1000"))+ggtitle("observed -*de novo* assemblies and reference annotations")+theme(plot.title = element_markdown(size=11))



###proTRAC_p0.05 only 

t1<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/protrac/new-cluster_bed/Canton-S_p0.05_gapped_protrac_summary.forR")
t2<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/protrac/new-cluster_bed/DGRP-732_p0.05_gapped_protrac_summary.forR")
t3<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/protrac/new-cluster_bed/Iso1_p0.05_gapped_protrac_summary.forR")
t4<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/protrac/new-cluster_bed/Oregon-R_p0.05_gapped_protrac_summary.forR")
t5<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/protrac/new-cluster_bed/Pi2_p0.05_gapped_protrac_summary.forR")

t<-rbind(t1,t2,t3,t4,t5)

names(t)<-c("count","id","TE","region")

t$id2<-paste(t$id,t$TE,sep="+")

cl<-subset(t,region=="cluster")
cl<-subset(cl,select=c("count","id2"))

rest<-subset(t,region!="cluster")
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

for (sid in unique(t$TE)) { 
  i <- t$TE == sid
  a = mean(t$cluster[i])
  b = mean(t$noncluster[i])
  t$avrcl[i] = a
  t$avrrest[i] = b
}

t<-subset(t,select = c("TE","avrcl","avrrest"))
t<-unique(t)

t<-left_join(t,info,by="TE")

#including AF threshold
t<-subset(t,AF!="NA")##remove missing AFs
t<-subset(t,AF<=0.2)



cC<-t

cC$cluster<-log10(cC$avrcl+1)
cC$noncluster<-log10(cC$avrrest+1)


gP<-ggplot(cC,aes(x=noncluster,y=cluster))+geom_point()+stat_cor(cor.coef.name="tau",method = "kendall", label.x = 1, label.y = 2,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-cluster insertions")+ylab("cluster insertions")+scale_x_continuous(breaks=c(0,1.041393,2.004321,3.000434),labels=c("0","10","100","1000"))+scale_y_continuous(breaks=c(0,1.041393,2.004321,3.000434),labels=c("0","10","100","1000"))+ggtitle("observed -*de novo* assemblies and annotations")+theme(plot.title = element_markdown(size=11))







gcor<-ggarrange(gS, gPop, gC, gP,
                labels = c("A", "B", "C","D"),
                ncol = 2, nrow = 2)
```

    ## Warning: Removed 10 rows containing non-finite values (stat_cor).

    ## Warning: Removed 10 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 10 rows containing missing values (geom_point).

``` r
plot(gcor)
```

![](main-cor_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/correlation_main.pdf",width=9,height=6)
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/correlation_main.png",width=9,height=6)
```
