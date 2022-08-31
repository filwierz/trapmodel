assembly-quality
================
Filip Wierzbicki
8/30/2022

\#CQ/ScQ analysis

``` bash
#example command (for Canton-S assembly) how busco was run:
export BUSCO_CONFIG_FILE=/home/linuxbrew/.linuxbrew/Cellar/busco/5.0.0/libexec/config.ini
busco -m genome -o Canton-S -i Canton-S.fasta -l diptera_odb10 -c 60 --augustus
#from full_table output (.tsv files) bed files were generated:
for i in *tsv;do n=${i%.tsv};cat $i|awk 'NR>3'|awk '$2=="Complete"'|awk '{print $3"\t"$4-1"\t"$5-1"\t"$1}' > ../busco_bed/${n}_busco.bed;done

#for cluster annotations (see obtain-clusters.Rmd)
cd /Volumes/Temp3/filip/trap_model/clusterscore/cluster_bed
for i in *bed;do python /Users/filipwierzbicki/Desktop/trap_model/github/trapmodel/helper-scripts/bed-cleaner_gaps.py --bed $i > ../${i};done
cd ../../bam
all:
nohup sh -c 'for i in *bam;do n=${i%.sort.bam};samtools view $i|python /Volumes/Temp3/filip/programs/roberts_scripts/cuscoquality-code/cluster-coverage-median.py --sam - --cluster ../cluster_bed/${n}_cluster.bed --reference ../busco_bed/${n}_busco.bed --min-mq 60 --min-len 5000 --output-cluster ../coverage/${n}_ml5k_mq60.cluster --output-reference ../coverage/${n}_ml5k_mq60.busco > ../coverage/${n}_ml5k_mq60.summary;done' &
nohup sh -c 'for i in *bam;do n=${i%.sort.bam};samtools view $i|python /Volumes/Temp3/filip/programs/roberts_scripts/cuscoquality-code/cluster-softclipcoverage-median.py --sam - --cluster ../cluster_bed/${n}_cluster.bed --reference ../busco_bed/${n}_busco.bed --min-mq 60 --min-len 5000 --output-cluster ../softclip/${n}_ml5k_mq60.cluster --output-reference ../softclip/${n}_ml5k_mq60.busco > ../softclip/${n}_ml5k_mq60.summary;done' &

gapless:
for i in *.bed;do cat $i|awk '$5=="1000"' > ../cluster_bed_gapless/${i};done
nohup sh -c 'for i in *bam;do n=${i%.sort.bam};samtools view $i|python /Volumes/Temp3/filip/programs/roberts_scripts/cuscoquality-code/cluster-coverage-median.py --sam - --cluster ../cluster_bed_gapless/${n}_cluster.bed --reference ../busco_bed/${n}_busco.bed --min-mq 60 --min-len 5000 --output-cluster ../coverage_gapless/${n}_ml5k_mq60.cluster --output-reference ../coverage_gapless/${n}_ml5k_mq60.busco > ../coverage_gapless/${n}_ml5k_mq60.summary;done' &
nohup sh -c 'for i in *bam;do n=${i%.sort.bam};samtools view $i|python /Volumes/Temp3/filip/programs/roberts_scripts/cuscoquality-code/cluster-softclipcoverage-median.py --sam - --cluster ../cluster_bed_gapless/${n}_cluster.bed --reference ../busco_bed/${n}_busco.bed --min-mq 60 --min-len 5000 --output-cluster ../softclip_gapless/${n}_ml5k_mq60.cluster --output-reference ../softclip_gapless/${n}_ml5k_mq60.busco > ../softclip_gapless/${n}_ml5k_mq60.summary;done' &

CQ-ScQ all/gapless:
cd /Volumes/Temp3/filip/trap_model/clusterscore/coverage/
for i in *_ml5k_mq60.summary;do n=${i%_ml5k_mq60.summary};head -1 $i|awk -v a="$n" '{print $3,a}';done > /Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/CQ-all.txt
cd ../softclip
for i in *_ml5k_mq60.summary;do n=${i%_ml5k_mq60.summary};head -1 $i|awk -v a="$n" '{print $4,a}';done > /Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/ScQ-all.txt
cd ../coverage_gapless
for i in *_ml5k_mq60.summary;do n=${i%_ml5k_mq60.summary};head -1 $i|awk -v a="$n" '{print $3,a}';done > /Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/CQ-gapless.txt
cd ../softclip_gapless
for i in *_ml5k_mq60.summary;do n=${i%_ml5k_mq60.summary};head -1 $i|awk -v a="$n" '{print $4,a}';done > /Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/ScQ-gapless.txt
```

\#Genome size

``` bash
#from samtools fasta index (.fai file)
cd /Users/filipwierzbicki/Desktop/trap_model/data/core/assemblies
for i in *fai;do n=${i%.fasta.fai};cat $i|awk -v a="$n" '{SUM+=$2}END{print SUM,a}';done > /Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/genomesize.txt
```

\#n50

``` bash
#from samtools fasta index (.fai file)
cd /Users/filipwierzbicki/Desktop/trap_model/data/core/assemblies
for i in *fai;do n=${i%.fasta.fai};cat $i|awk '{SUM+=$2}END{print SUM/2}'|while read h;do cat $i|awk '{print $2}'|sort -n -r|awk -v v="$h" '0>=(v-=$1)'|head -1|awk -v a="$n" '{print $0,a}';done;done > /Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/n50.txt
```

\#all cusco

``` bash
cd /Users/filipwierzbicki/Desktop/trap_model/analysis/cluster/cusco/bwasw/cluster_bed
for i in *bed;do n=${i%_cluster.bed};cat $i|wc -l|awk -v a="$n" '{print $0/85*100,a}';done > /Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/cusco.txt
```

\#ungapped cusco

``` bash
cd /Users/filipwierzbicki/Desktop/trap_model/analysis/cluster/cusco/bwasw/cluster_bed
for i in *bed;do n=${i%_cluster.bed};cat $i|awk '$5=="1000"'|wc -l|awk -v a="$n" '{print $0/85*100,a}';done > /Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/ucusco.txt
```

\#BUSCO scores

``` bash
#from BUSCO short summary output:
cd /Volumes/Temp3/filip/trap_model/clusterscore/busco-summary
for i in *.txt;do n=${i#short_summary.specific.diptera_odb10.};n=${n%.txt};cat $i|head -9|tail -1|awk -F "%" '{print $1}'|awk -F ":" '{print $2}'|awk -v a="$n" '{print $0,a}';done > /Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/busco.txt
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
#CQ/ScQ:
ca<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/CQ-all.txt")
names(ca)<-c("CQ-all","assembly")
sa<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/ScQ-all.txt")
names(sa)<-c("ScQ-all","assembly")
cu<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/CQ-gapless.txt")
names(cu)<-c("CQ-ungapped","assembly")
su<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/ScQ-gapless.txt")
names(su)<-c("ScQ-ungapped","assembly")

#Genome size:
gs<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/genomesize.txt")
gs$V1<-gs$V1/1000000
names(gs)<-c("assembly-size (MB)","assembly")

#n50:
n50<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/n50.txt")
n50$V1<-n50$V1/1000000
names(n50)<-c("N50 (MB)","assembly")

#(gapped)cusco:
c<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/cusco.txt")
names(c)<-c("g.CUSCO (%)","assembly")
#ungapped cusco:
uc<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/ucusco.txt")
names(uc)<-c("u.CUSCO (%)","assembly")

#BUSCO score:
b<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/output/busco.txt")
names(b)<-c("BUSCO (%)","assembly")

t<-inner_join(gs,n50,by="assembly")
t<-subset(t,select=c("assembly","assembly-size (MB)","N50 (MB)"))
t<-inner_join(t,b,by="assembly")
t<-inner_join(t,c,by="assembly")
t<-inner_join(t,uc,by="assembly")
t[,-1] <-round(t[,-1],1)
t<-inner_join(t,ca,by="assembly")
t<-inner_join(t,cu,by="assembly")
t<-inner_join(t,sa,by="assembly")
t<-inner_join(t,su,by="assembly")
t[,-1] <-round(t[,-1],3)
t<-t(t)
knitr::kable(t)
```

|                    |          |          |       |          |       |
| :----------------- | :------- | :------- | :---- | :------- | :---- |
| assembly           | Canton-S | DGRP-732 | Iso1  | Oregon-R | Pi2   |
| assembly-size (MB) | 149.1    | 141.6    | 143.7 | 136.3    | 167.8 |
| N50 (MB)           | 28.2     | 25.7     | 25.3  | 25.1     | 30.9  |
| BUSCO (%)          | 99.4     | 99.0     | 99.5  | 99.2     | 99.3  |
| g.CUSCO (%)        | 95.3     | 92.9     | 97.6  | 91.8     | 97.6  |
| u.CUSCO (%)        | 78.8     | 75.3     | 85.9  | 76.5     | 81.2  |
| CQ-all             | 0.069    | 0.107    | 0.141 | 0.139    | 0.066 |
| CQ-ungapped        | 0.141    | 0.120    | 0.192 | 0.194    | 0.099 |
| ScQ-all            | 0.512    | 0.145    | 0.099 | 0.329    | 0.351 |
| ScQ-ungapped       | 0.395    | 0.162    | 0.193 | 0.317    | 0.282 |

``` r
write.table(t,"/Users/filipwierzbicki/Desktop/trap_model/analysis/assembly-quality/final-table.txt",quote=FALSE,sep="\t")
```
