obtain-clusters
================
Filip Wierzbicki
4/6/2022

This scripts contains the pipeline for obtaining piRNA cluster and
reference region (control region) annotations.

Note: make-directory are usually not included and also change directory
only for rough orientation included but not completely

``` bash

cd /Users/filipwierzbicki/Desktop/trap_model/data/core/assemblies

for i in *fasta;do n=${i%.fasta};bwa bwasw $i /Users/filipwierzbicki/Desktop/AssembTech/cusco/resources/neighboring_flanks-mixed_240719.fasta > /Users/filipwierzbicki/Desktop/trap_model/analysis/cluster/cusco/bwasw/map/${n}.sam;done

for i in *fasta;do n=${i%.fasta};python /Users/filipwierzbicki/Documents/PhD_Project/Bioinformatics/cuscoquality-code/find-polyN.py --fasta $i > /Users/filipwierzbicki/Desktop/trap_model/analysis/cluster/cusco/bwasw/polyn/${n}.polyn;done

cd /Users/filipwierzbicki/Desktop/trap_model/analysis/cluster/cusco/bwasw/map

for i in *sam;do n=${i%.sam};python /Users/filipwierzbicki/Documents/PhD_Project/Bioinformatics/cuscoquality-code/cusco.py --sam $i --pic /Users/filipwierzbicki/Desktop/AssembTech/cusco/resources/mixed_flanks_mainChr_230719 --polyn ../polyn/${n}.polyn --output-cb ../cluster_bed/${n}_cluster.bed > ../cusco/${n}.cusco;done

for i in *cluster.bed;do n=${i};cat $i|awk '$5=="1000"'|sort -k1,1 -k2,2n > gapless/${n};done

for i in *bed;do n=${i};bedtools merge -i $i -o distinct -c 4,5,6 > distinct/${n};done
```

``` bash

#note here gapped cluster annotations are considered 

cd /Users/filipwierzbicki/Desktop/trap_model/analysis/cluster/cusco/bwasw/cluster_bed

for i in *bed;do n=${i};cat $i|sort -k1,1 -k2,2n|bedtools merge -i - -o distinct -c 4,5,6 > ../all_clusters/distinct-cluster_bed/${n};done

cp -r ../all_clusters/distinct-cluster_bed /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/ref/ref_recover/combine/gapped_cusco-distinct

cd /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/ref/ref_recover/combine/gapped_cusco-distinct

for i in *_cluster.bed;do n=${i%_cluster.bed}; cat $i ../TAS/${n}_cluster.bed|sort -k1,1 -k2,2n|bedtools merge -i - -o distinct -c 4,5,6 > /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/ref_recover/combined-distinct/$i;done

cd /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/ref_recover/combined-distinct

for i in *_cluster.bed;do n=${i%_cluster.bed};python /Users/filipwierzbicki/Desktop/trap_model/scripts/reference_recover.py --clu $i --hel /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/ref/ref_recover/helper_files/${n}_contigs.txt --polyn /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/ref/ref_recover/polyn/${n}.polyn > ../ref_bed/${n}_ref.bed;done
```

``` bash
cd /Users/filipwierzbicki/Desktop/trap_model/data/core/assemblies

for i in *fasta;do n=${i%.fasta};bwa bwasw $i /Users/filipwierzbicki/Desktop/trap_model/analysis/cluster/TAS/playground/first_genes/first_genes.fasta > /Users/filipwierzbicki/Desktop/trap_model/analysis/cluster/TAS/playground/first_genes/${n}.sam;done

cd /Users/filipwierzbicki/Desktop/trap_model/analysis/cluster/TAS/playground/first_genes/bwasw

for i in *sam;do n=${i%.sam};python /Users/filipwierzbicki/Desktop/trap_model/scripts/TAS_recover.py --sam $i --polyn ../polyn/${n}.polyn --fai ../fai/${n}.fasta.fai > ../cluster_bed/${n}_cluster.bed;done

cp -r ../cluster_bed /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/TAS_cluster_bed

cd /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/TAS_cluster_bed

for i in *_cluster.bed;do n=${i%_cluster.bed};cat $i /Users/filipwierzbicki/Desktop/trap_model/analysis/cluster/cusco/bwasw/cluster_bed/gapless/distinct/$i|sort -k1,1 -k2,2n|bedtools merge -i - -o distinct -c 4,5,6 > ../combined-distinct/$i;done
```

``` bash
#trimming for Canton-S, Pi2, Iso1:
/Users/fschwarz/.local/bin/cutadapt -j 20 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATTTTATCTCGTATGC -o trimmed/Pi2.fastq Pi2.fastq
 
/Users/fschwarz/.local/bin/cutadapt -j 20 -a AGATCGGAAGAGCACACGTCT -o trimmed/DGRP-732.fastq DGRP-732.fastq
  
#Oregon-R is trimmed

for i in *.fastq;do n=${i%.fastq};cat $i| paste - - - -|awk 'length($2) >17 && length($2) <36'|tr "\t" "\n" > ../trim-sf/${n}.sf.fastq;done

for i in *.fastq;do n=${i%.sf.fastq}; novoalign -d /Volumes/Temp3/filip/trap_model/proTRAC/ref/dmel_all_RNA.nvi -f $i -F STDFQ -o SAM -o FullNW -r RANDOM > map/${n}.sam;done

filter for piRNA reads (te und unmapped between 23-29nt):
for i in *.sam; do n=${i%.sam}; cat $i|grep -v '^@'|awk '$3 !~ /_mRNA|_miRNA|_rRNA|_snRNA|_snoRNA|_tRNA/'|awk '{if ((22<length($10)) && (length($10)<30)) print $0}'|awk '{print "@" $1; print $10; print "+" $1; print $11}' > piRNAs/${n}_piRNA.fastq; done

for i in *fastq;do n=${i};perl /Volumes/Temp3/filip/trap_model/proTRAC/TBr2/TBr2_collapse.pl -i $i -o collapsed/${n};done

for i in *fastq;do n=${i};perl /Volumes/Temp3/filip/trap_model/proTRAC/TBr2/TBr2_duster.pl -i $i;done

in for loop with specific assembly: 
perl sRNAmapper.pl -input piRNAs.fasta -genome genome.fasta -alignments best

for i in *map;do n=${i%_piRNA.fastq.no-dust.map}; cat $i|sort -k5,5|awk '{print $1,$2,$3,$4,$6,$7,$5}'|uniq -u -f 6|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$5"\t"$6}'|sort -k 1,1 -k 2n,2 > unique_V3_sorted/${n}.map;done

for i in *map;do n=${i%.map};perl /Volumes/Temp3/filip/trap_model/proTRAC/proTRAC_2.4.4.pl -map $i -genome ../genomes/${n}.fasta -pdens 0.05 -pimin 23 -pimax 29 -1Tor10A 0.3 -clsize 5000 -clstrand 0.5;done

for i in p0.*/proTRAC_*;do n=${i%.map*};n=${n#p0.*/proTRAC_};st=${i%/proTRAC*};cat $i/clusters.gtf|awk -v a=$n '{print $1"\t"$4-1"\t"$5-1"\t"$12}'|sed 's/;//g' > minimal_cluster_bed/${n}_${st}.bed;done

cd minimal_cluster_bed

for i in *.bed;do n=${i%.bed};g=${n%_*};python /Users/filipwierzbicki/Desktop/trap_model/scripts/merge_proTRAC_clusters.py --bed $i --polyn ../../merge_annotations/polyn/${g}.polyn > ../cluster_bed/${n}_cluster.bed;done
```

``` bash
cd /Users/filipwierzbicki/Desktop/trap_model/analysis/cluster/cusco/bwasw/cluster_bed

mkdir distinct
for i in *bed;do n=${i};cat $i|sort -k1,1 -k2,2n|bedtools merge -i - -o distinct -c 4,5,6 > distinct/${n};done

cd /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/TAS_cluster_bed

for i in *_cluster.bed;do n=${i%_cluster.bed};cat $i /Users/filipwierzbicki/Desktop/trap_model/analysis/cluster/cusco/bwasw/cluster_bed/distinct/$i|sort -k1,1 -k2,2n|bedtools merge -i - -o distinct -c 4,5,6 > ../gapped_combined-distinct/$i;done
```

``` bash
cd /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas/gapped_combined-distinct

for i in *_cluster.bed;do n=${i%_cluster.bed};cat $i /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/protrac/protrac_gapped_cluster_bed/${n}_p0.05_cluster.bed|sort -k1,1 -k2,2n|bedtools merge -i - -o distinct -c 4,5,6 > /Users/filipwierzbicki/Desktop/trap_model/analysis/abu/cusco_tas_protrac/gapped_combined-distinct/$i;done
```
