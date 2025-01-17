supp-Onecode
================
Filip Wierzbicki
10/14/2022

This script contains the pipeline to merge fragmented repeatmasker
annotations using onecodetofindthemall.pl: Bailly-Bechet M., Haudry A. &
Lerat E. (2014) “One code to find them all”: a Perl tool to conveniently
parse RepeatMasker output files . Mobile DNA 5:13.

The purpose of this supplementary analysis is to what extend the merging
would change ob observed patterns in the correlation between cluster and
non-cluster insertions and the distribution of cluster insertions.

``` bash

for i in *.out;do n=${i%.fasta.out};/Volumes/Temp3/filip/programs/Onecodetofindthemall/build_dictionary.pl --rm /Volumes/Temp3/filip/trap_model/whole-genome/repeatmasker/${i} --unknown > /Volumes/Temp3/filip/trap_model/onecode/file_dictionary/${n};done



for i in *.out;do n=${i%.fasta.out};/Volumes/Temp3/filip/programs/Onecodetofindthemall/one_code_to_find_them_all.pl --rm $i --ltr file_dictionary/${n} --length /Volumes/Temp3/filip/cluster20/resources/new_dmel_consensus_TEs.lengths --unknown;done


for i in *.out;do n=${i%.fasta.out};cat ${i}_*.elem_sorted.csv|grep '^###' > process/${n}_elem-compiled.rm;done


cp to local:
cd /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/onecode

for i in *_cluster.bed;do n=${i%_cluster.bed};mkdir ${n};python /Users/filipwierzbicki/Desktop/trap_model/github/trapmodel/helper-scripts/assembly_TE-abundance_3types.py --clu $i --rm /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/onecode-rm/${n}_elem-compiled.rm --output ${n}/ --sample ${n} --approach gapped_cusco_tas-onecode --minlen 100 --maxdiv 10.0 ;done
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


#######
###For population frequency Info based on Kofler et al. 2015 PLOS Genetics
info1<-read.table("/Users/filipwierzbicki/Desktop/evolution_cluster/temp/TEfamInfo_correct")
names(info1)<-c("name","TE","order","AF","popins")
###exclude somatically regulated TEs based on Malone et al. 2009 Cell
info1<-subset(info1,name!="gypsy10"&name!="gypsy"&name!="ZAM"&name!="gtwin"&name!="gypsy5"&name!="Tabor")
info<-subset(info1,select=c("TE","AF"))
info$AF<-round(info$AF,digits = 1)


t1<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/onecode/Canton-S_gapped_cusco_tas-onecode_summary.forR")
t2<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/onecode/DGRP-732_gapped_cusco_tas-onecode_summary.forR")
t3<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/onecode/Iso1_gapped_cusco_tas-onecode_summary.forR")
t4<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/onecode/Oregon-R_gapped_cusco_tas-onecode_summary.forR")
t5<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct/onecode/Pi2_gapped_cusco_tas-onecode_summary.forR")


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


real<-ggplot(ht, aes(x=cluster, y=indsrel))+ geom_vline( xintercept =alq,col="red")+ geom_vline( xintercept =auq,col="red") + geom_bar(stat="identity")+ylab("frequency")+xlab("cluster insertions")#+ geom_vline( xintercept =alq,col="red") + geom_vline( xintercept =auq,col="red")#+xlim(-1,160)+ylim(0,0.35)


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



gC<-ggplot(cC,aes(x=noncluster,y=cluster))+geom_point()+stat_cor(cor.coef.name="tau", method = "kendall", label.x = 0.5, label.y = 2,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-cluster insertions")+ylab("cluster insertions")+scale_x_continuous(breaks=c(0,1.041393,2.004321,3.000434),labels=c("0","10","100","1000"))+scale_y_continuous(breaks=c(0,1.041393,2.004321,3.000434),labels=c("0","10","100","1000"))


gOC<-ggarrange(gC,real,
               labels = c("A","B"),
               ncol = 2, nrow = 1)
plot(gOC)
```

![](supp_OneCode_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/OneCode_cuscoTAS.png",width=8,height=6)
ggsave("/Users/filipwierzbicki/Desktop/trap_model/analysis/abu/figures/OneCode_cuscoTAS.pdf",width=8,height=6)
```
