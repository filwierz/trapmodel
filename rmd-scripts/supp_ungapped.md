supp\_ungapped-CUSCO+TAS
================
Filip Wierzbicki
7/19/2022

Additional analysis of TE abundance in gapless CUSCO clusters and TAS.

``` bash

cd /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/combined-distinct
for i in *_cluster.bed;do n=${i%_cluster.bed};mkdir ${n};python /Users/filipwierzbicki/Desktop/trap_model/scripts/assembly_TE-abundance_3types.py --clu $i --ref ../ref_recover/ref_bed/${n}_ref.bed --rm /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/whole-genome/repeatmasker/${n}.fasta.out --output ${n}/ --sample ${n} --approach gapless_cusco_tas --minlen 100 --maxdiv 10.0 ;done
```

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


t1<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/combined-distinct/Canton-S_gapless_cusco_tas_summary.forR")
t2<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/combined-distinct/DGRP-732_gapless_cusco_tas_summary.forR")
t3<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/combined-distinct/Iso1_gapless_cusco_tas_summary.forR")
t4<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/combined-distinct/Oregon-R_gapless_cusco_tas_summary.forR")
t5<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/combined-distinct/Pi2_gapless_cusco_tas_summary.forR")

tq<-rbind(t1,t2,t3,t4,t5)
names(tq)<-c("count","id","TE","region")

ht<-tq
ht$id2<-paste(ht$id,ht$TE,sep="+")

cl<-subset(ht,region=="cluster")
cl<-subset(cl,select=c("count","id2"))

rest<-subset(ht,region!="cluster")
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
ht<-cr
ht<-left_join(ht,info,by="TE")

#including AF threshold
ht<-subset(ht,AF!="NA")##remove missing AFs
ht<-subset(ht,AF<=0.2)

ht[is.na(ht)] <- 0

ht$id3<-paste(ht$cluster,ht$TE,sep="_")

for (sid in unique(ht$id3)) { 
  i <- ht$id3 == sid
  a = nrow(subset(ht,id3==sid))
  ht$sum[i] = a
}

ht<-subset(ht,select=c("cluster","TE","sum"))
ht<-unique(ht)
for (sid in unique(ht$cluster)) { 
  i <- ht$cluster == sid
  b = sum(ht$sum[i])
  ht$inds[i] = b
}

ht<-subset(ht,select=c("cluster","inds"))
ht<-unique(ht)

ht$indsrel<-ht$inds/sum(ht$inds)


real<-ggplot(ht, aes(x=cluster, y=indsrel)) + geom_bar(stat="identity")+ylab("frequency of individuals")+xlab("number of cluster insertions")#+xlim(0,30)#+ geom_vline( xintercept =alq,col="red") + geom_vline( xintercept =auq,col="red")


t<-tq
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



gC<-ggplot(cC,aes(x=noncluster,y=cluster))+geom_point()+stat_cor(method = "kendall", label.x = 1, label.y = 2,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-cluster insertions")+ylab("cluster insertions")+scale_x_continuous(breaks=c(0,1,2,3),labels=c("0","9","99","999"))+scale_y_continuous(breaks=c(0,1,2,3),labels=c("0","9","99","999"))


####combine:
gGAP<-ggarrange(gC,real,
                labels = c("A","B"),
                ncol = 2, nrow = 1)
plot(gGAP)
```

![](supp_ungapped_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/supp_ungapped.png",width=8,height=3)
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/supp_ungapped.pdf",width=8,height=3)
```
