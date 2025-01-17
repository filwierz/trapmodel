main-simulation
================
Filip Wierzbicki
6/17/2022

``` bash
#cmds to run invade-simulation on vetlinux05
nohup sh -c 'python ../../../scripts/simstorm-constant-u-seed.py --number 300 --threads 60 --output output-constant_u- --invade ../../../invade-v0808.jar --silent' &
nohup sh -c 'python ../../scripts/simstorm-variable-u_unif.py --number 300 --threads 30 --output output-varu_unif- --invade ../../invade-v0808.jar --silent' &
nohup sh -c 'python ../../scripts/simstorm-site04xMULT.py --number 300 --threads 12 --output site04xMULT- --invade ../../invade-v0808.jar --silent' &
#followed by concatenation of output-* and tally-* files to load into R
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
gentx=2000

tally<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/constant-u/run3-seed/combined/tally-constant_u")

names(tally)<-c("replicate","generation","type","abundance","number")

t<-subset(tally,generation==gentx)
t<-subset(t,type=="cluster")

for (sid in unique(t$abundance)) { 
  i <- t$abundance == sid
  a = sum(t$number[i])
  t$sum[i] = a
} 

ts<-unique(subset(t,select=c("abundance","sum")))
ts$rel<-ts$sum/sum(ts$sum)

histo<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("trap insertions")+xlim(0,30)

t<-subset(tally,generation==gentx)
t<-subset(t,type=="ref")

for (sid in unique(t$abundance)) { 
  i <- t$abundance == sid
  a = sum(t$number[i])
  t$sum[i] = a
} 

ts<-unique(subset(t,select=c("abundance","sum")))
ts$rel<-ts$sum/sum(ts$sum)

historef<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("reference insertions")+xlim(0,30)


output<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/constant-u/run3-seed/combined/output-constant_u")


output<-output[,-28]
names(output)<-c("replicate","generation","delim1",
            "fwt",  "w",    "tes",  "popfreq",  "fixed","delim2",
            "fwcli",    "cluins",   "cluins_popfreq","cluins_fixed",    "phase","delim3",
            "fwrefi",   "refins",   "refins_popfreq", "refins_fixed","delim4",
            "novel",    "sites",    "clusites", "tes_stdev" ,"cluins_stdev" ,"fw0", "w_min","popsize")



ts<-subset(output,generation==gentx)

ts$tesc<-ts$tes-ts$cluins
ts$tesr<-ts$tes-ts$refins#-ts$cluins
cluster<-subset(ts,select = c("tesc","cluins"))
names(cluster)<-c("global","local")
cluster$type<-("cluster")
reference<-subset(ts,select = c("tesr","refins"))
names(reference)<-c("global","local")
reference$type<-c("reference")
cr<-rbind(cluster,reference)


correlationC<-ggplot(cluster, aes(x=global, y=local)) + geom_point()+stat_cor(method = "kendall", cor.coef.name="tau", label.x = 50, label.y = 12.5,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-trap insertions")+ylab("trap insertions")
correlationR<-ggplot(reference, aes(x=global, y=local)) + geom_point()+stat_cor(method = "kendall", cor.coef.name="tau", label.x = 50, label.y = 20.0,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-reference insertions")+ylab("reference insertions")

####trajectories:

trajectories<-ggplot(output, aes(x = generation, y = tes, group = replicate,col= factor(phase))) +
  geom_line() +geom_vline(xintercept =  gentx,col="black")+scale_color_manual(values = c("green", "yellow", "red"))+theme(legend.position = "None")+ylab("TE abundance")+xlim(0,7500) 



gconsu<-ggarrange(trajectories,histo,historef, correlationC, correlationR,
             #labels = c("A","B","C","D","E"),
             ncol = 5, nrow = 1)
```

    ## Warning: Removed 7500 row(s) containing missing values (geom_path).

    ## Warning: Removed 1 rows containing missing values (geom_bar).

    ## Warning: Removed 2 rows containing missing values (position_stack).

    ## Warning: Removed 2 rows containing missing values (geom_bar).

``` r
############################
#variable time points

output<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/constant-u/run3-seed/combined/output-constant_u")



output<-output[,-28]
names(output)<-c("replicate","generation","delim1",
                 "fwt", "w",    "tes",  "popfreq",  "fixed","delim2",
                 "fwcli",   "cluins",   "cluins_popfreq","cluins_fixed",    "phase","delim3",
                 "fwrefi",  "refins",   "refins_popfreq", "refins_fixed","delim4",
                 "novel",   "sites",    "clusites", "tes_stdev" ,"cluins_stdev" ,"fw0", "w_min","popsize")




tout<-subset(output,select=c("replicate","generation","tes","cluins","refins","phase"))
t<-tout
t<-subset(t,generation>2500&generation<7500) ##standard range
#t<-subset(t,generation>1000&generation<5000)  #special range
t$id<-paste(t$replicate,t$generation,sep="_")

reps<-sample(unique(t$replicate),300)#130
reps<-as.data.frame(reps)
names(reps)<-c("replicate")

id<-subset(t,select = c("id"))
id <- id[0,]

for (sid in reps$replicate) { 
  tx<-subset(t,t$replicate==sid)
  x<-sample(unique(tx$generation),1)
  x<-paste(sid,x,sep="_")
  idx<-as.data.frame(x)
  id<-rbind(id,idx)
} 
names(id)<-c("id")

t<-inner_join(t,id,by="id")
ts<-t

ts$tesc<-ts$tes-ts$cluins
ts$tesr<-ts$tes-ts$refins#-ts$cluins
cluster<-subset(ts,select = c("tesc","cluins"))
names(cluster)<-c("global","local")
cluster$type<-("cluster")
reference<-subset(ts,select = c("tesr","refins"))
names(reference)<-c("global","local")
reference$type<-c("reference")
cr<-rbind(cluster,reference)


correlationC<-ggplot(cluster, aes(x=global, y=local)) + geom_point()+stat_cor(method = "kendall", cor.coef.name="tau", label.x = 50, label.y = 15,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-trap insertions")+ylab("trap insertions")
correlationR<-ggplot(reference, aes(x=global, y=local)) + geom_point()+stat_cor(method = "kendall", cor.coef.name="tau", label.x = 50, label.y = 23.5,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-reference insertions")+ylab("reference insertions")

####trajectories:
t<-tout
idtr<-id
idtr$replicate<-gsub("_.*","",idtr$id)
idtr$generation<-gsub(".*_","",idtr$id)
t$replicate<-as.character(t$replicate)
t<-inner_join(t,idtr,by="replicate")

##for loop to get generation.x==generation.y!!
#tp<-unique(subset(t,select=c("replicate","tes","genera")))
#+ geom_point(data=h, aes(x=sample,y=ppm), color='red',size=3)  

h<-subset(t,generation.x==generation.y)

trajectories<-ggplot(t, aes(x = generation.x, y = tes, group = replicate,col= factor(phase))) +
  geom_line() +scale_color_manual(values = c("green", "yellow", "red"))+theme(legend.position = "None")+ylab("TE abundance") +
  geom_point(data=h, aes(x = generation.x, y = tes), color='black',size=1)+xlab("generation")+xlim(0,7500) 

#+geom_point(aes(x=2500,y=100))

#####


tally<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/constant-u/run3-seed/combined/tally-constant_u")

names(tally)<-c("replicate","generation","type","abundance","number")

t<-tally
t<-subset(t,type=="cluster")
t$id<-paste(t$replicate,t$generation,sep="_")
t<-inner_join(t,id,by="id")
for (sid in unique(t$abundance)) { 
  i <- t$abundance == sid
  a = sum(t$number[i])
  t$sum[i] = a
} 

ts<-unique(subset(t,select=c("abundance","sum")))
ts$rel<-ts$sum/sum(ts$sum)

histo<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("trap insertions")+xlim(0,30)



###historef:

t<-tally
t<-subset(t,type=="ref")
t$id<-paste(t$replicate,t$generation,sep="_")
t<-inner_join(t,id,by="id")
for (sid in unique(t$abundance)) { 
  i <- t$abundance == sid
  a = sum(t$number[i])
  t$sum[i] = a
} 

ts<-unique(subset(t,select=c("abundance","sum")))
ts$rel<-ts$sum/sum(ts$sum)

historef<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("reference insertions")+xlim(0,30)



gvart<-ggarrange(trajectories,histo,historef, correlationC, correlationR,
             #labels = c("A","B","C","D","E"),
             ncol = 5, nrow = 1)
```

    ## Warning: Removed 7500 row(s) containing missing values (geom_path).

    ## Warning: Removed 1 rows containing missing values (geom_bar).
    
    ## Warning: Removed 1 rows containing missing values (geom_bar).

``` r
###############
#variable U:

gentx=2000

tally<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/varu-unif/combined/tally-varu_unif.txt")



names(tally)<-c("replicate","generation","type","abundance","number")

t<-subset(tally,generation==gentx)
t<-subset(t,type=="cluster")

for (sid in unique(t$abundance)) { 
  i <- t$abundance == sid
  a = sum(t$number[i])
  t$sum[i] = a
} 

ts<-unique(subset(t,select=c("abundance","sum")))
ts$rel<-ts$sum/sum(ts$sum)

histo<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("trap insertions")+xlim(0,30)

t<-subset(tally,generation==gentx)
t<-subset(t,type=="ref")

for (sid in unique(t$abundance)) { 
  i <- t$abundance == sid
  a = sum(t$number[i])
  t$sum[i] = a
} 

ts<-unique(subset(t,select=c("abundance","sum")))
ts$rel<-ts$sum/sum(ts$sum)

historef<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("reference insertions")+xlim(0,30)




output<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/varu-unif/combined/output-varu_unif.txt")




output<-output[,-28]
names(output)<-c("replicate","generation","delim1",
            "fwt",  "w",    "tes",  "popfreq",  "fixed","delim2",
            "fwcli",    "cluins",   "cluins_popfreq","cluins_fixed",    "phase","delim3",
            "fwrefi",   "refins",   "refins_popfreq", "refins_fixed","delim4",
            "novel",    "sites",    "clusites", "tes_stdev" ,"cluins_stdev" ,"fw0", "w_min","popsize")



ts<-subset(output,generation==gentx)

ts$tesc<-ts$tes-ts$cluins
ts$tesr<-ts$tes-ts$refins#-ts$cluins
cluster<-subset(ts,select = c("tesc","cluins"))
names(cluster)<-c("global","local")
cluster$type<-("cluster")
reference<-subset(ts,select = c("tesr","refins"))
names(reference)<-c("global","local")
reference$type<-c("reference")
cr<-rbind(cluster,reference)


correlationC<-ggplot(cluster, aes(x=global, y=local)) + geom_point()+stat_cor(method = "kendall", cor.coef.name="tau", label.x = 50, label.y = 14.0,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-trap insertions")+ylab("trap insertions")
correlationR<-ggplot(reference, aes(x=global, y=local)) + geom_point()+stat_cor(method = "kendall", cor.coef.name="tau", label.x = 50, label.y = 18.7,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-reference insertions")+ylab("reference insertions")

####trajectories:

trajectories<-ggplot(output, aes(x = generation, y = tes, group = replicate,col= factor(phase))) +
  geom_line() +geom_vline(xintercept =  gentx,col="black")+scale_color_manual(values = c("green", "yellow", "red"))+theme(legend.position = "None")+ylab("TE abundance")+xlim(0,7500) 



gvaru<-ggarrange(trajectories,histo,historef, correlationC, correlationR,
             #labels = c("A","B","C","D","E"),
             ncol = 5, nrow = 1)
```

    ## Warning: Removed 7500 row(s) containing missing values (geom_path).
    
    ## Warning: Removed 1 rows containing missing values (geom_bar).
    
    ## Warning: Removed 1 rows containing missing values (geom_bar).

``` r
##############
#negative selection:

gentx=2000

tally<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/site04xMULT/tally-site04xMULT.txt")


names(tally)<-c("replicate","generation","type","abundance","number")

t<-subset(tally,generation==gentx)
t<-subset(t,type=="cluster")

for (sid in unique(t$abundance)) { 
  i <- t$abundance == sid
  a = sum(t$number[i])
  t$sum[i] = a
} 

ts<-unique(subset(t,select=c("abundance","sum")))
ts$rel<-ts$sum/sum(ts$sum)

histo<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("trap insertions")+xlim(0,30)

t<-subset(tally,generation==gentx)
t<-subset(t,type=="ref")

for (sid in unique(t$abundance)) { 
  i <- t$abundance == sid
  a = sum(t$number[i])
  t$sum[i] = a
} 

ts<-unique(subset(t,select=c("abundance","sum")))
ts$rel<-ts$sum/sum(ts$sum)

historef<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("reference insertions")+xlim(0,30)


output<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/site04xMULT/output_site04xMULT.txt")


output<-output[,-28]
names(output)<-c("replicate","generation","delim1",
            "fwt",  "w",    "tes",  "popfreq",  "fixed","delim2",
            "fwcli",    "cluins",   "cluins_popfreq","cluins_fixed",    "phase","delim3",
            "fwrefi",   "refins",   "refins_popfreq", "refins_fixed","delim4",
            "novel",    "sites",    "clusites", "tes_stdev" ,"cluins_stdev" ,"fw0", "w_min","popsize")



ts<-subset(output,generation==gentx)

ts$tesc<-ts$tes-ts$cluins
ts$tesr<-ts$tes-ts$refins#-ts$cluins
cluster<-subset(ts,select = c("tesc","cluins"))
names(cluster)<-c("global","local")
cluster$type<-("cluster")
reference<-subset(ts,select = c("tesr","refins"))
names(reference)<-c("global","local")
reference$type<-c("reference")
cr<-rbind(cluster,reference)


correlationC<-ggplot(cluster, aes(x=global, y=local)) + geom_point()+stat_cor(method = "kendall", cor.coef.name="tau", label.x = 25, label.y = 11.5,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-trap insertions")+ylab("trap insertions")
correlationR<-ggplot(reference, aes(x=global, y=local)) + geom_point()+stat_cor(method = "kendall", cor.coef.name="tau", label.x = 25, label.y = 11,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-reference insertions")+ylab("reference insertions")

####trajectories:

trajectories<-ggplot(output, aes(x = generation, y = tes, group = replicate,col= factor(phase))) +
  geom_line() +geom_vline(xintercept =  gentx,col="black")+scale_color_manual(values = c("green", "yellow", "red"))+theme(legend.position = "None")+ylab("TE abundance")+xlim(0,7500) 



gneg<-ggarrange(trajectories,histo,historef, correlationC, correlationR,
             #labels = c("A","B","C","D","E"),
             ncol = 5, nrow = 1)
```

    ## Warning: Removed 7500 row(s) containing missing values (geom_path).
    
    ## Warning: Removed 1 rows containing missing values (geom_bar).
    
    ## Warning: Removed 1 rows containing missing values (geom_bar).

``` r
#############
#combine plots:

g<-ggarrange(gconsu,gvart,gvaru,gneg,
             labels = c("A","B","C","D"),
             ncol = 1, nrow = 4)

plot(g)
```

![](main-simulation_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggsave("/Users/filipwierzbicki/Desktop/trap_model/figures/simulations/simulations_main.pdf",height = 12,width = 12)
ggsave("/Users/filipwierzbicki/Desktop/trap_model/figures/simulations/simulations_main.png",height = 12,width = 12)
```
