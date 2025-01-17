example\_butterfly-graph
================
Filip Wierzbicki
9/30/2022

``` bash
cd /Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/sort/
##remove "_RaGOO" from chromosome names to make it compatible with the downstream script:
samtools view Canton-S.sort.bam|sed 's/_RaGOO//g' > forGraph/Canton-S.sam

cd forGraph
##sum up piRNA reads for the 3R chromosome:
cat Canton-S.sam|awk '$3=="3R"'|python2.7 /Volumes/Temp3/filip/programs/roberts_scripts/te-tools-code/piRNA/graph-piRNA-distri-window-chromo.py --sam - --min-mq 5 --window-size 1 --min-length 23 --max-length 29 --sample-id Canton-S > forR/Canton-S_mq5-ws1bp_3Ronly.txt
```

``` r
library(ggplot2)
theme_set(theme_classic())
t<-read.table("/Volumes/Temp3/filip/trap_model/butterfly/filtered-reads/map/sort/forGraph/forR/Canton-S_mq5-ws1bp_3Ronly.txt")
names(t)<-c("sid","chr","start","end","id","signal")
t<-subset(t, id!="bias")
t[t$id=="as",]$signal<-t[t$id=="as",]$signal*-1

###example of a butterfly at a F-element
##3R_RaGOO:8,156,742-8,164,735
t<-subset(t,chr=="3R")
xmax=8156742
xmin=8164735


p<- ggplot(t,aes(x=start,y=signal,color=id))+geom_segment(aes(xend=start),yend=0)+xlim(xmax,xmin)+ylim(-250,250)
plot(p)
```

    ## Warning: Removed 356258 rows containing missing values (geom_segment).

![](example_butterfly-graph_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggsave("/Users/filipwierzbicki/Desktop/trap_model/figures/example_butterfly/graph_fromR.pdf",width=8,height=5)
```

    ## Warning: Removed 356258 rows containing missing values (geom_segment).

``` r
##this subfigure was manually added to the figure using Inkscape
```
