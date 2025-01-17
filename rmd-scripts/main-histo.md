main\_histo
================
Filip Wierzbicki
3/1/2022

This scripts contains the whole pipeline for the main histogram figure
of TE distribution in clusters.

Note that simulated TE copy numbers are per diploid, while counts based
populationTE2 are per haploid and assemblies represent a single
chromosome (however could be chimeric).

``` bash
cd /Volumes/Temp3/filip/trap_model/popTE2/pipeline

java -Xmx16g -jar /Volumes/Temp3/filip/programs/popte2-v1.10.04.jar ppileup --bam Canton-S.sort.bam --bam DGRP-732.sort.bam --bam Iso1.sort.bam --bam Oregon-R.sort.bam --bam Pi2.sort.bam --map-qual 15 --hier /Volumes/Temp3/filip/cluster20/popte2/release5/rel5/ref/new_dmel_132cons_hier --output ../pipeline/ppileup.gz

java -Xmx6g -jar /Volumes/Temp3/filip/programs/popte2-v1.10.04.jar identifySignatures --ppileup ppileup.gz --mode separate --output sep_signatures.msm --min-count 2 --signature-window minimumSampleMedian --min-valley minimumSampleMedian

java -jar /Volumes/Temp3/filip/programs/popte2-v1.10.04.jar frequency --ppileup ppileup.gz --signature sep_signatures.msm --output sep_freqsignatures

java -jar /Volumes/Temp3/filip/programs/popte2-v1.10.04.jar filterSignatures --input sep_freqsignatures --output sep_filtered.freqsignatures --min-coverage 2 --max-otherte-count 2 --max-structvar-count 2

java -jar /Volumes/Temp3/filip/programs/popte2-v1.10.04.jar pairupSignatures --signature sep_filtered.freqsignatures --ref-genome /Volumes/Temp3/filip/cluster20/popte2/release5/rel5/ref/rel5_masked_132-dmel-cons_combined.fasta --hier /Volumes/Temp3/filip/cluster20/popte2/release5/rel5/ref/new_dmel_132cons_hier --output sep_teinsertions

#move to local machine

python ../../../scripts/popTE2-cluster_VS_rest.py --popte2 sep_teinsertions --pic resources/rel5_brenecluster.gff --minfreq 0.3 > sep.mf30_brenecluster

cat sep.mf30_brenecluster|awk '{print $1"_"$5"_"$10}'|sort|uniq -c|awk -F "_" '{print $1,$2,$3}' > sep.mf30_brenecluster.forR
#haploid:
cat sep.mf30_brenecluster|awk '$9>=0.6'|awk '{print $1"_"$5"_"$10}'|sort|uniq -c|awk -F "_" '{print $1,$2,$3}'|awk '{print $1,$2,$3,$4,$2"_"$3"_"$4}' > haploid_counts/sep.mf30_brenecluster-fixed
cat sep.mf30_brenecluster|awk '$9<0.6'|awk '{print $1"_"$5"_"$10}'|sort|uniq -c|awk -F "_" '{print $1,$2,$3}'|awk '{print $1/2,$2,$3,$4,$2"_"$3"_"$4}' > haploid_counts/sep.mf30_brenecluster-seg
```

``` bash
cd /home/vetlinux05/Filip/trap_model/simulations/storm2/constant-u/run3-seed
nohup sh -c 'python ../../../scripts/simstorm-constant-u-seed.py --number 300 --threads 60 --output output-constant_u- --invade ../../../invade-v0808.jar --silent' &
#move the folder 'run3-seed' to local machine
cd run3-seed
mkdir combined 
cat tally* > combined/tally-constant_u
```

``` bash
cd /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct
for i in *_cluster.bed;do n=${i%_cluster.bed};mkdir ${n};python /Users/filipwierzbicki/Desktop/trap_model/github/trapmodel/helper-scripts/assembly_TE-abundance_3types.py --clu $i --ref ../ref_recover/ref_bed/${n}_ref.bed --rm /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/whole-genome/repeatmasker/${n}.fasta.out --output ${n}/ --sample ${n} --approach gapped_cusco_tas --minlen 100 --maxdiv 10.0 ;done
```

``` bash
cd /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/protrac/new-cluster_bed
for i in *_cluster.bed;do n=${i%_cluster.bed};mkdir ${n};python /Users/filipwierzbicki/Desktop/trap_model/github/trapmodel/helper-scripts/assembly_TE-abundance_3types.py --clu $i --rm /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/whole-genome/repeatmasker/${n}.fasta.out --output ${n}/ --sample ${n}_p0.05 --approach gapped_protrac --minlen 100 --maxdiv 10.0 ;done
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
library(ggtext)

#simulations:
t<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/constant-u/run3-seed/combined/tally-constant_u")

names(t)<-c("replicate","generation","type","abundance","number")
t$abundance<-t$abundance/2 #for halpoid abundance

t<-subset(t,generation==2000)
t<-subset(t,type=="cluster")

for (sid in unique(t$abundance)) { 
  i <- t$abundance == sid
  a = sum(t$number[i])
  t$sum[i] = a
} 

ts<-unique(subset(t,select=c("abundance","sum")))
ts$rel<-ts$sum/sum(ts$sum)

lq=0
uq=0
tso<-ts[order(ts$abundance),]
for (row in 1:nrow(tso)) { 
  lq=lq+tso$rel[row]
  
  if(lq>0.01){
    alq=(tso$abundance[row]+tso$abundance[row-1])/2
    break
  }
}
for (row in 1:nrow(tso)) { 
  uq=uq+tso$rel[row]
  
  if(uq>0.99){
    auq=(tso$abundance[row]+tso$abundance[row+1])/2
    break
  }
}

simu<-ggplot(ts, aes(x=abundance, y=rel)) + geom_vline( xintercept =alq,col="red")+ geom_vline( xintercept =auq,col="red") + geom_bar(stat = "identity")+ylab("frequency")+xlab("cluster insertions")+ggtitle("expected - simulations")+theme(plot.title = element_markdown(size=11))


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

#new supplementary table for overview of used TEs:
towrite<-info1
towrite$AFr<-round(towrite$AF,digits = 1)
towrite<-subset(towrite,AFr<=0.2)
towrite<-subset(towrite,select=c("name","TE","order","AF"))
towrite$AF<-round(towrite$AF*100,digits = 1)
names(towrite)<-c("TE-family","seqID","order","population-frequency")
write.table(towrite,file="/Users/filipwierzbicki/Desktop/trap_model/data/other/suptab_usedTEs.txt",quote = FALSE, row.names = FALSE)


##cusco+TAS
t1<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Canton-S_gapped_cusco_tas_summary.forR")
t2<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/DGRP-732_gapped_cusco_tas_summary.forR")
t3<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Iso1_gapped_cusco_tas_summary.forR")
t4<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Oregon-R_gapped_cusco_tas_summary.forR")
t5<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/Pi2_gapped_cusco_tas_summary.forR")


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


real<-ggplot(ht, aes(x=cluster, y=indsrel))+ geom_vline( xintercept =alq,col="red")+ geom_vline( xintercept =auq,col="red") + geom_bar(stat="identity")+ylab("frequency")+xlab("cluster insertions")+xlim(-1,31)+ggtitle("observed -*de novo* assemblies and reference annotations")+theme(plot.title = element_markdown(size=11))


##ref (not used currently)

reftq<-tq
reftq$id2<-paste(reftq$id,reftq$TE,sep="+")

cl<-subset(reftq,region=="ref")
cl<-subset(cl,select=c("count","id2"))

rest<-subset(reftq,region=="non-cluster")
rest<-subset(rest,select=c("count","id2"))

cr<-full_join(cl,rest,by="id2")
names(cr)<-c("cluster","id2","noncluster")

cr$id<-gsub("\\+.*","",cr$id2)
cr$TE<-gsub(".*\\+","",cr$id2)

t<-cr
t<-left_join(t,info,by="TE")

#including AF threshold
t<-subset(t,AF!="NA")##remove missing AFs
t<-subset(t,AF<=0.2)


t[is.na(t)] <- 0

cR<-t
cR$id3<-paste(cR$cluster,cR$TE,sep="_")

for (sid in unique(cR$id3)) { 
  i <- cR$id3 == sid
  a = nrow(subset(cR,id3==sid))
  cR$sum[i] = a
}

cR<-subset(cR,select=c("cluster","TE","sum"))
cR<-unique(cR)
for (sid in unique(cR$cluster)) { 
  i <- cR$cluster == sid
  b = sum(cR$sum[i])
  cR$inds[i] = b
}

cR<-subset(cR,select=c("cluster","inds"))
cR<-unique(cR)

cR$indsrel<-cR$inds/sum(cR$inds)


ref<-ggplot(cR, aes(x=cluster, y=indsrel))+ geom_vline( xintercept =alq,col="red")+ geom_vline( xintercept =auq,col="red") + geom_bar(stat="identity")+ylab("frequency")+xlab("number of reference insertions")#+xlim(-1,160)+ylim(0,0.35)


###whole genome (not used currently)

tq$id2<-paste(tq$id,tq$TE,sep="_")
for (sid in unique(tq$id2)) { 
  i <- tq$id2 == sid
  b = sum(tq$count[i])
  tq$insertions[i] = b
}
tqu<-subset(tq,select=c("id","TE","insertions"))
tqu<-unique(tqu)

tqu<-left_join(tqu,info,by="TE")

#including AF threshold
tqu<-subset(tqu,AF!="NA")##remove missing AFs
tqu<-subset(tqu,AF<=0.2)

tqu[is.na(tqu)] <- 0

tqu$id3<-paste(tqu$insertions,tqu$TE,sep="_")

for (sid in unique(tqu$id3)) { 
  i <- tqu$id3 == sid
  a = nrow(subset(tqu,id3==sid))
  tqu$sum[i] = a
}

tqu<-subset(tqu,select=c("insertions","TE","sum"))
tqu<-unique(tqu)
for (sid in unique(tqu$insertions)) { 
  i <- tqu$insertions == sid
  b = sum(tqu$sum[i])
  tqu$inds[i] = b
}

tqu<-subset(tqu,select=c("insertions","inds"))
tqu<-unique(tqu)

tqu$indsrel<-tqu$inds/sum(tqu$inds)


whole<-ggplot(tqu, aes(x=insertions, y=indsrel)) + geom_bar(stat="identity")+ylab("frequency")+xlab("number of genomic insertions")#+xlim(-1,160)+ylim(0,0.35)


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

rest<-subset(t,region!="cl") ##### currently no third region class (e.g. ref) considered
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
pt<-cr

pt<-left_join(pt,infoP,by="TE")
#including AF threshold
pt<-subset(pt,AF!="NA")##remove missing AFs
pt<-subset(pt,AF<=0.2)

pt[is.na(pt)] <- 0

id<-c(1,2,3,4,5)
strain<-c("Canton-S","DGRP-732","Iso1","Oregon-R","Pi2")
ids<-data.frame(id,strain)

pt<-merge(pt,ids,by="id")
pt$id<-pt$strain

pt$id3<-paste(pt$cluster,pt$TE,sep="_")

for (sid in unique(pt$id3)) { 
  i <- pt$id3 == sid
  a = nrow(subset(pt,id3==sid))
  pt$sum[i] = a
}

pt<-subset(pt,select=c("cluster","TE","sum"))
pt<-unique(pt)
for (sid in unique(pt$cluster)) { 
  i <- pt$cluster == sid
  b = sum(pt$sum[i])
  pt$inds[i] = b
}

pt<-subset(pt,select=c("cluster","inds"))
pt<-unique(pt)

pt$indsrel<-pt$inds/sum(pt$inds)


popTE2<-ggplot(pt, aes(x=cluster, y=indsrel))+ geom_vline( xintercept =alq,col="red")+ geom_vline( xintercept =auq,col="red") + geom_bar(stat="identity")+ylab("frequency")+xlab("cluster insertions")+ggtitle("observed - reference assembly and reference annotations")+theme(plot.title = element_markdown(size=11))




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

rest<-subset(t,region!="cluster") ##### currently no third region class (e.g. ref) considered
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

prt<-cr
prt<-left_join(prt,info,by="TE")

#including AF threshold
prt<-subset(prt,AF!="NA")##remove missing AFs
prt<-subset(prt,AF<=0.2)

prt[is.na(prt)] <- 0

prt$id3<-paste(prt$cluster,prt$TE,sep="_")

for (sid in unique(prt$id3)) { 
  i <- prt$id3 == sid
  a = nrow(subset(prt,id3==sid))
  prt$sum[i] = a
}

prt<-subset(prt,select=c("cluster","TE","sum"))
prt<-unique(prt)
for (sid in unique(prt$cluster)) { 
  i <- prt$cluster == sid
  b = sum(prt$sum[i])
  prt$inds[i] = b
}

prt<-subset(prt,select=c("cluster","inds"))
prt<-unique(prt)

prt$indsrel<-prt$inds/sum(prt$inds)


proT<-ggplot(prt, aes(x=cluster, y=indsrel))+ geom_vline( xintercept =alq,col="red")+ geom_vline( xintercept =auq,col="red") + geom_bar(stat="identity")+ylab("frequency")+xlab("cluster insertions")+xlim(-1,31)+ggtitle("observed -*de novo* assemblies and annotations")+theme(plot.title = element_markdown(size=11))




g<-ggarrange(simu, popTE2, real, proT,
             labels = c("A", "B","C","D"),
             ncol = 2, nrow = 2)
```

    ## Warning: Removed 62 rows containing missing values (position_stack).

    ## Warning: Removed 1 rows containing missing values (geom_bar).

    ## Warning: Removed 13 rows containing missing values (position_stack).

    ## Warning: Removed 1 rows containing missing values (geom_bar).

``` r
plot(g)
```

![](main-histo_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/histogram_main.pdf",width=9,height=6)
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/histogram_main.png",width=9,height=6)
```
