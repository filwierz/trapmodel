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

th=15
#th=30

b1<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/other-parameter/TE/Canton-S_w500.txt")
#b1<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/output-V2/Canton-S.txt")
names(b1)<-c("TE","chr","start","end","ls","la","rs","ra","strain")
b2<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/other-parameter/TE/DGRP-732_w500.txt")
#b2<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/output-V2/DGRP-732.txt")
names(b2)<-c("TE","chr","start","end","ls","la","rs","ra","strain")
b3<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/other-parameter/TE/Iso1_w500.txt")
#b3<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/output-V2/Iso1.txt")
names(b3)<-c("TE","chr","start","end","ls","la","rs","ra","strain")
b4<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/other-parameter/TE/Oregon-R_w500.txt")
#b4<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/output-V2/Oregon-R.txt")
names(b4)<-c("TE","chr","start","end","ls","la","rs","ra","strain")
b5<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/other-parameter/TE/Pi2_w500.txt")
#b5<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/output-V2/Pi2.txt")
names(b5)<-c("TE","chr","start","end","ls","la","rs","ra","strain")


a1<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Canton-S_gapped_cusco_tas_summary.forR")
a2<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/DGRP-732_gapped_cusco_tas_summary.forR")
a3<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Iso1_gapped_cusco_tas_summary.forR")
a4<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Oregon-R_gapped_cusco_tas_summary.forR")
a5<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Pi2_gapped_cusco_tas_summary.forR")
names(a1)<-c("count","strain","TE","region")
names(a2)<-c("count","strain","TE","region")
names(a3)<-c("count","strain","TE","region")
names(a4)<-c("count","strain","TE","region")
names(a5)<-c("count","strain","TE","region")

c1<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Canton-S_cluster.bed")
names(c1)<-c("chr","start","end","cluster","ig","ir")
c2<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/DGRP-732_cluster.bed")
names(c2)<-c("chr","start","end","cluster","ig","ir")
c3<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Iso1_cluster.bed")
names(c3)<-c("chr","start","end","cluster","ig","ir")
c4<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Oregon-R_cluster.bed")
names(c4)<-c("chr","start","end","cluster","ig","ir")
c5<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Pi2_cluster.bed")
names(c5)<-c("chr","start","end","cluster","ig","ir")

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


for (sid in unique(ic$TE)) { 
  i <- ic$TE == sid
  x = nrow(subset(ic,ic$TE==sid))
  ic$butterfly[i] = x
}
bfs<-subset(ic,select = c("TE","butterfly","strain"))
bfs<-unique(bfs)
ac<-subset(a1,region=="cluster")
ac<-subset(ac,select=c("TE","count"))
TE1<-full_join(bfs,ac,by="TE")

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

for (sid in unique(ic$TE)) { 
  i <- ic$TE == sid
  x = nrow(subset(ic,ic$TE==sid))
  ic$butterfly[i] = x
}
bfs<-subset(ic,select = c("TE","butterfly","strain"))
bfs<-unique(bfs)
ac<-subset(a2,region=="cluster")
ac<-subset(ac,select=c("TE","count"))
TE2<-full_join(bfs,ac,by="TE")

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

for (sid in unique(ic$TE)) { 
  i <- ic$TE == sid
  x = nrow(subset(ic,ic$TE==sid))
  ic$butterfly[i] = x
}
bfs<-subset(ic,select = c("TE","butterfly","strain"))
bfs<-unique(bfs)
ac<-subset(a3,region=="cluster")
ac<-subset(ac,select=c("TE","count"))
TE3<-full_join(bfs,ac,by="TE")

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

for (sid in unique(ic$TE)) { 
  i <- ic$TE == sid
  x = nrow(subset(ic,ic$TE==sid))
  ic$butterfly[i] = x
}
bfs<-subset(ic,select = c("TE","butterfly","strain"))
bfs<-unique(bfs)
ac<-subset(a4,region=="cluster")
ac<-subset(ac,select=c("TE","count"))
TE4<-full_join(bfs,ac,by="TE")

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

for (sid in unique(ic$TE)) { 
  i <- ic$TE == sid
  x = nrow(subset(ic,ic$TE==sid))
  ic$butterfly[i] = x
}
bfs<-subset(ic,select = c("TE","butterfly","strain"))
bfs<-unique(bfs)
ac<-subset(a5,region=="cluster")
ac<-subset(ac,select=c("TE","count"))
TE5<-full_join(bfs,ac,by="TE")



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



#summary that needs to be kept:
t<-rbind(TE1,TE2,TE3,TE4,TE5)

t[is.na(t)] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
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


gcor<-ggplot(t,aes(x=clusterlog,y=butterflylog))+geom_point()+stat_cor(method = "pearson", label.x = 0.25, label.y = 1.0,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("cluster insertions")+ylab("butterfly insertions")+scale_x_continuous(breaks=c(0,1,2),labels=c("0","9","99"))+scale_y_continuous(breaks=c(0,1,2),labels=c("0","9","99"))

plot(gcor)
```

![](supp_butterfly-cor_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/supp_butterfly-cor.pdf",width=7,height=6)
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/supp_butterfly-cor.png",width=7,height=6)

a1r<-subset(a1,region!="cluster")

a1r$id2<-paste(a1r$strain,a1r$TE,sep="+")
for (sid in unique(a1r$id2)) { 
  i <- a1r$id2 == sid
  a = sum(a1r$count[i])
  a1r$sum[i] = a
}
a1r<-subset(a1r,select=c("sum","TE","strain"))
a1r<-unique(a1r)

t1a<-full_join(a1r,TE1,by="TE")
t1a[is.na(t1a)] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated
    
    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
t1a<-left_join(t1a,info,by="TE")

#including AF threshold
t1a<-subset(t1a,AF!="NA")##remove missing AFs
t1a<-subset(t1a,AF<=0.2)



a2r<-subset(a2,region!="cluster")

a2r$id2<-paste(a2r$strain,a2r$TE,sep="+")
for (sid in unique(a2r$id2)) { 
  i <- a2r$id2 == sid
  a = sum(a2r$count[i])
  a2r$sum[i] = a
}
a2r<-subset(a2r,select=c("sum","TE","strain"))
a2r<-unique(a2r)

t2a<-full_join(a2r,TE2,by="TE")
t2a[is.na(t2a)] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
t2a<-left_join(t2a,info,by="TE")

#including AF threshold
t2a<-subset(t2a,AF!="NA")##remove missing AFs
t2a<-subset(t2a,AF<=0.2)




a3r<-subset(a3,region!="cluster")

a3r$id2<-paste(a3r$strain,a3r$TE,sep="+")
for (sid in unique(a3r$id2)) { 
  i <- a3r$id2 == sid
  a = sum(a3r$count[i])
  a3r$sum[i] = a
}
a3r<-subset(a3r,select=c("sum","TE","strain"))
a3r<-unique(a3r)

t3a<-full_join(a3r,TE3,by="TE")
t3a[is.na(t3a)] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
t3a<-left_join(t3a,info,by="TE")

#including AF threshold
t3a<-subset(t3a,AF!="NA")##remove missing AFs
t3a<-subset(t3a,AF<=0.2)




a4r<-subset(a4,region!="cluster")

a4r$id2<-paste(a4r$strain,a4r$TE,sep="+")
for (sid in unique(a4r$id2)) { 
  i <- a4r$id2 == sid
  a = sum(a4r$count[i])
  a4r$sum[i] = a
}
a4r<-subset(a4r,select=c("sum","TE","strain"))
a4r<-unique(a4r)

t4a<-full_join(a4r,TE4,by="TE")
t4a[is.na(t4a)] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated
    
    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
t4a<-left_join(t4a,info,by="TE")

#including AF threshold
t4a<-subset(t4a,AF!="NA")##remove missing AFs
t4a<-subset(t4a,AF<=0.2)



a5r<-subset(a5,region!="cluster")

a5r$id2<-paste(a5r$strain,a5r$TE,sep="+")
for (sid in unique(a5r$id2)) { 
  i <- a5r$id2 == sid
  a = sum(a5r$count[i])
  a5r$sum[i] = a
}
a5r<-subset(a5r,select=c("sum","TE","strain"))
a5r<-unique(a5r)

t5a<-full_join(a5r,TE5,by="TE")
t5a[is.na(t5a)] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated
    
    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
t5a<-left_join(t5a,info,by="TE")

#including AF threshold
t5a<-subset(t5a,AF!="NA")##remove missing AFs
t5a<-subset(t5a,AF<=0.2)



ta<-rbind(t1a,t2a,t3a,t4a,t5a)
ta0<-subset(ta,count==0)

hier<-read.table("/Users/filipwierzbicki/Desktop/trap_model/data/other/new_dmel_132cons_hier",header =TRUE)
hier<-subset(hier,select=c("id","family"))
names(hier)<-c("TE","family")
ta0<-left_join(ta0,hier,by="TE")
ta0<-subset(ta0,select=c("family","sum","butterfly","count","strain.x"))
names(ta0)<-c("family","non-cluster","other-SL","cluster","strain")

write.table(ta0,"/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/butterfly/output/missing-clusterinsertions.txt",quote=FALSE)
```