supp\_butterfly-onChr
================
Filip Wierzbicki
9/30/2022

This script generates a plot of butterfly signatures along the genome in
each strains.

The code to generate the input files can be found in the
butterfly\_cor.Rmd

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

th=15


##the filter with this Info table is not used here since its rather a demonstration than the trap model test:
##(Note in the analysis we filter for germline TEs with low population frequencies)
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

b2ch<-read.table("/Users/filipwierzbicki/Desktop/trap_model/data/core/chromosome-names/DGRP-732.tsv")
b4ch<-read.table("/Users/filipwierzbicki/Desktop/trap_model/data/core/chromosome-names/Oregon-R.tsv")

#butterfly signatures:
b1<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/other-parameter/TE/Canton-S_w500.txt")
names(b1)<-c("TE","chr","start","end","ls","la","rs","ra","strain")
b1$pos<-b1$start+((b1$end-b1$start)/2)
b1$chromo<-gsub("_RaGOO","",b1$chr)
b1$signal<-b1$ls+b1$la+b1$rs+b1$ra

b2<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/other-parameter/TE/DGRP-732_w500.txt")
b2<-left_join(b2,b2ch,by="V2")
names(b2)<-c("TE","accession","start","end","ls","la","rs","ra","strain","chr")
b2<-subset(b2,select=c("TE","chr","start","end","ls","la","rs","ra","strain"))
b2$pos<-b2$start+((b2$end-b2$start)/2)
b2$chromo<-b2$chr
b2$signal<-b2$ls+b2$la+b2$rs+b2$ra

b3<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/other-parameter/TE/Iso1_w500.txt")
names(b3)<-c("TE","chr","start","end","ls","la","rs","ra","strain")
#chromosomes are good to go!
b3$pos<-b3$start+((b3$end-b3$start)/2)
b3$chromo<-b3$chr
b3$signal<-b3$ls+b3$la+b3$rs+b3$ra

b4<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/other-parameter/TE/Oregon-R_w500.txt")
b4<-left_join(b4,b4ch,by="V2")
names(b4)<-c("TE","accession","start","end","ls","la","rs","ra","strain","chr")
b4<-subset(b4,select=c("TE","chr","start","end","ls","la","rs","ra","strain"))
b4$pos<-b4$start+((b4$end-b4$start)/2)
b4$chromo<-b4$chr
b4$signal<-b4$ls+b4$la+b4$rs+b4$ra

b5<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/other-parameter/TE/Pi2_w500.txt")
names(b5)<-c("TE","chr","start","end","ls","la","rs","ra","strain")
b5$pos<-b5$start+((b5$end-b5$start)/2)
b5$chromo<-gsub("_RaGOO","",b5$chr)
b5$signal<-b5$ls+b5$la+b5$rs+b5$ra

#cluster annotations to exclude from butterfly signatures:
c1<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Canton-S_cluster.bed")
names(c1)<-c("chr","start","end","cluster","ig","ir")
c1$chromo<-gsub("_RaGOO","",c1$chr)

c2<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/DGRP-732_cluster.bed")
c2<-left_join(c2,b2ch, by = c("V1" = "V2"))
names(c2)<-c("chr","start","end","cluster","ig","ir","chromo")

c3<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Iso1_cluster.bed")
names(c3)<-c("chr","start","end","cluster","ig","ir")
#chromosomes are good to go!
c3$chromo<-c3$chr

c4<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Oregon-R_cluster.bed")
c4<-left_join(c4,b4ch, by = c("V1" = "V2"))
names(c4)<-c("chr","start","end","cluster","ig","ir","chromo")


c5<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Pi2_cluster.bed")
names(c5)<-c("chr","start","end","cluster","ig","ir")
c5$chromo<-gsub("_RaGOO","",c5$chr)

c1$start<-c1$start+1
c1$end<-c1$end+1

c2$start<-c2$start+1
c2$end<-c2$end+1

c3$start<-c3$start+1
c3$end<-c3$end+1

c4$start<-c4$start+1
c4$end<-c4$end+1

c5$start<-c5$start+1
c5$end<-c5$end+1


b1$id<-paste(b1$TE,b1$chr,b1$start,b1$end,sep = "_")
b2$id<-paste(b2$TE,b2$chr,b2$start,b2$end,sep = "_")
b3$id<-paste(b3$TE,b3$chr,b3$start,b3$end,sep = "_")
b4$id<-paste(b4$TE,b4$chr,b4$start,b4$end,sep = "_")
b5$id<-paste(b5$TE,b5$chr,b5$start,b5$end,sep = "_")
###
#excluding signatures from clusters: 
ic=b1[FALSE,]

for(i in 1:nrow(b1)) {
  row <- b1[i,]
  
  for(k in 1:nrow(c1)){
    line <- c1[k,]
    if ((line$chr==row$chr && line$start<=row$start&&row$start<=line$end)||(line$chr==row$chr && line$start<=row$end&&row$end<=line$end)){
      ic[nrow(ic) + 1,] <- row
      break
    }
    
  }
}

nc<-anti_join(b1,ic,by="id")

ic=nc[FALSE,]

for(i in 1:nrow(nc)) {
  row <- nc[i,]
  if (row$la>=th && row$rs>=th){
    ic[nrow(ic) + 1,] <- row
  }
}

tb<-subset(ic,chromo=="2L"|chromo=="2R"|chromo=="3L"|chromo=="3R"|chromo=="X")
tb$chr_f = factor(tb$chromo, levels=c('X','2L','2R','3L','3R'))

tb<-subset(tb,select=c("chr_f","pos","signal"))
tb$type<-c("stand-alone")
c1$pos<-c1$start+((c1$end-c1$start)/2)
c1x<-subset(c1,chromo=="2L"|chromo=="2R"|chromo=="3L"|chromo=="3R"|chromo=="X")
c1x$chr_f = factor(c1x$chromo, levels=c('X','2L','2R','3L','3R'))
c1x<-subset(c1x,select=c("chr_f","pos"))
c1x$signal<-500
c1x$type<-c("cluster")
tb1<-rbind(tb,c1x)
#

ic=b2[FALSE,]

for(i in 1:nrow(b2)) {
  row <- b2[i,]
  
  for(k in 1:nrow(c2)){
    line <- c2[k,]
    if ((line$chr==row$chr && line$start<=row$start&&row$start<=line$end)||(line$chr==row$chr && line$start<=row$end&&row$end<=line$end)){
      ic[nrow(ic) + 1,] <- row
      break
    }
    
  }
}

nc<-anti_join(b2,ic,by="id")

ic=nc[FALSE,]

for(i in 1:nrow(nc)) {
  row <- nc[i,]
  if (row$la>=th && row$rs>=th){
    ic[nrow(ic) + 1,] <- row
  }
}

tb<-subset(ic,chromo=="2L"|chromo=="2R"|chromo=="3L"|chromo=="3R"|chromo=="X")
tb$chr_f = factor(tb$chromo, levels=c('X','2L','2R','3L','3R'))

tb<-subset(tb,select=c("chr_f","pos","signal"))
tb$type<-c("stand-alone")
c2$pos<-c2$start+((c2$end-c2$start)/2)
c1x<-subset(c2,chromo=="2L"|chromo=="2R"|chromo=="3L"|chromo=="3R"|chromo=="X")
c1x$chr_f = factor(c1x$chromo, levels=c('X','2L','2R','3L','3R'))
c1x<-subset(c1x,select=c("chr_f","pos"))
c1x$signal<-500
c1x$type<-c("cluster")
tb2<-rbind(tb,c1x)
#

ic=b3[FALSE,]

for(i in 1:nrow(b3)) {
  row <- b3[i,]
  
  for(k in 1:nrow(c3)){
    line <- c3[k,]
    if ((line$chr==row$chr && line$start<=row$start&&row$start<=line$end)||(line$chr==row$chr && line$start<=row$end&&row$end<=line$end)){
      ic[nrow(ic) + 1,] <- row
      break
    }
    
  }
}

nc<-anti_join(b3,ic,by="id")

ic=nc[FALSE,]

for(i in 1:nrow(nc)) {
  row <- nc[i,]
  if (row$la>=th && row$rs>=th){
    ic[nrow(ic) + 1,] <- row
  }
}

tb<-subset(ic,chromo=="2L"|chromo=="2R"|chromo=="3L"|chromo=="3R"|chromo=="X")
tb$chr_f = factor(tb$chromo, levels=c('X','2L','2R','3L','3R'))

tb<-subset(tb,select=c("chr_f","pos","signal"))
tb$type<-c("stand-alone")
c3$pos<-c3$start+((c3$end-c3$start)/2)
c1x<-subset(c3,chromo=="2L"|chromo=="2R"|chromo=="3L"|chromo=="3R"|chromo=="X")
c1x$chr_f = factor(c1x$chromo, levels=c('X','2L','2R','3L','3R'))
c1x<-subset(c1x,select=c("chr_f","pos"))
c1x$signal<-500
c1x$type<-c("cluster")
tb3<-rbind(tb,c1x)
#

ic=b4[FALSE,]

for(i in 1:nrow(b4)) {
  row <- b4[i,]
  
  for(k in 1:nrow(c4)){
    line <- c4[k,]
    if ((line$chr==row$chr && line$start<=row$start&&row$start<=line$end)||(line$chr==row$chr && line$start<=row$end&&row$end<=line$end)){
      ic[nrow(ic) + 1,] <- row
      break
    }
    
  }
}

nc<-anti_join(b4,ic,by="id")

ic=nc[FALSE,]

for(i in 1:nrow(nc)) {
  row <- nc[i,]
  if (row$la>=th && row$rs>=th){
    ic[nrow(ic) + 1,] <- row
  }
}

tb<-subset(ic,chromo=="2L"|chromo=="2R"|chromo=="3L"|chromo=="3R"|chromo=="X")
tb$chr_f = factor(tb$chromo, levels=c('X','2L','2R','3L','3R'))

tb<-subset(tb,select=c("chr_f","pos","signal"))
tb$type<-c("stand-alone")
c4$pos<-c4$start+((c4$end-c4$start)/2)
c1x<-subset(c4,chromo=="2L"|chromo=="2R"|chromo=="3L"|chromo=="3R"|chromo=="X")
c1x$chr_f = factor(c1x$chromo, levels=c('X','2L','2R','3L','3R'))
c1x<-subset(c1x,select=c("chr_f","pos"))
c1x$signal<-500
c1x$type<-c("cluster")
tb4<-rbind(tb,c1x)
#

ic=b5[FALSE,]

for(i in 1:nrow(b5)) {
  row <- b5[i,]
  
  for(k in 1:nrow(c5)){
    line <- c5[k,]
    if ((line$chr==row$chr && line$start<=row$start&&row$start<=line$end)||(line$chr==row$chr && line$start<=row$end&&row$end<=line$end)){
      ic[nrow(ic) + 1,] <- row
      break
    }
    
  }
}

nc<-anti_join(b5,ic,by="id")

ic=nc[FALSE,]

for(i in 1:nrow(nc)) {
  row <- nc[i,]
  if (row$la>=th && row$rs>=th){
    ic[nrow(ic) + 1,] <- row
  }
}

tb<-subset(ic,chromo=="2L"|chromo=="2R"|chromo=="3L"|chromo=="3R"|chromo=="X")
tb$chr_f = factor(tb$chromo, levels=c('X','2L','2R','3L','3R'))

tb<-subset(tb,select=c("chr_f","pos","signal"))
tb$type<-c("stand-alone")
c5$pos<-c5$start+((c5$end-c5$start)/2)
c1x<-subset(c5,chromo=="2L"|chromo=="2R"|chromo=="3L"|chromo=="3R"|chromo=="X")
c1x$chr_f = factor(c1x$chromo, levels=c('X','2L','2R','3L','3R'))
c1x<-subset(c1x,select=c("chr_f","pos"))
c1x$signal<-500
c1x$type<-c("cluster")
tb5<-rbind(tb,c1x)
tb1$strain<-c("Canton-S")
tb2$strain<-c("DGRP-732")
tb3$strain<-c("Iso1")
tb4$strain<-c("Oregon-R")
tb5$strain<-c("Pi2")
tb<-rbind(tb1,tb2,tb3,tb4,tb5)

pt<-ggplot(tb,aes(x=pos,y=log(signal),color=type))+geom_segment(aes(xend=pos),yend=0)+facet_grid(strain~chr_f, scales="free_x", space = "free_x")+xlab("position")+ylab("signal")+
  scale_x_continuous(breaks=c(0,5000000,10000000,15000000,20000000,25000000),labels=c("0","5m","10m","15m","20m","25m"))+
  theme(axis.text=element_text(size=6),axis.title = element_text(size=9),, strip.background = element_rect(colour="black",size=0.5),strip.text = element_text(size=8,margin = margin(0.05, 0.05, 0.05, 0.05, "cm")))
#+#ylim(0,6000)+
plot(pt)
```

![](supp_butterfly-onChr_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/butterfly-onChr_supp.png",width=7,height=6)
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/butterfly-onChr_supp.pdf",width=7,height=6)
```