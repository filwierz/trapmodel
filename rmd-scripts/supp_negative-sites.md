simulations\_supp-negative-sites
================
Filip Wierzbicki
6/28/2022

``` bash
#cmds to run invade-simulation on vetlinux05
nohup sh -c 'python ../../scripts/simstorm-sitex01-IF_weaker.py --number 300 --threads 12 --output sitex01- --invade ../../invade-v0808.jar --silent' &
nohup sh -c 'python ../../scripts/simstorm-sitex03-IF_weaker.py --number 300 --threads 12 --output sitex03- --invade ../../invade-v0808.jar --silent' &
nohup sh -c 'python ../../scripts/simstorm-sitex05-IF_weaker.py --number 300 --threads 12 --output sitex05- --invade ../../invade-v0808.jar --silent' &
nohup sh -c 'python ../../scripts/simstorm-sitex07-IF_weaker.py --number 300 --threads 12 --output sitex07- --invade ../../invade-v0808.jar --silent' &
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

####10%


tally<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/sitex01-weaker/tally-sitex01-weaker.txt")


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

histo<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("trap insertions")

t<-subset(tally,generation==gentx)
t<-subset(t,type=="ref")

for (sid in unique(t$abundance)) { 
  i <- t$abundance == sid
  a = sum(t$number[i])
  t$sum[i] = a
} 

ts<-unique(subset(t,select=c("abundance","sum")))
ts$rel<-ts$sum/sum(ts$sum)

historef<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("reference insertions")


output<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/sitex01-weaker/output-sitex01-weaker.txt")

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
correlationR<-ggplot(reference, aes(x=global, y=local)) + geom_point()+stat_cor(method = "kendall", cor.coef.name="tau", label.x = 50, label.y = 17.0,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-reference insertions")+ylab("reference insertions")

####trajectories:

trajectories<-ggplot(output, aes(x = generation, y = tes, group = replicate,col= factor(phase))) +
  geom_line() +geom_vline(xintercept =  gentx,col="black")+scale_color_manual(values = c("green", "yellow", "red"))+theme(legend.position = "None")+ylab("TE abundance")+xlim(0,7500) 



g10<-ggarrange(trajectories,histo,historef, correlationC, correlationR,
             #labels = c("A","B","C","D","E"),
             ncol = 5, nrow = 1)
```

    ## Warning: Removed 7500 row(s) containing missing values (geom_path).

``` r
####30%

tally<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/sitex03-weaker/tally-sitex03-weaker.txt")


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

histo<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("trap insertions")

t<-subset(tally,generation==gentx)
t<-subset(t,type=="ref")

for (sid in unique(t$abundance)) { 
  i <- t$abundance == sid
  a = sum(t$number[i])
  t$sum[i] = a
} 

ts<-unique(subset(t,select=c("abundance","sum")))
ts$rel<-ts$sum/sum(ts$sum)

historef<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("reference insertions")


output<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/sitex03-weaker/output-sitex03-weaker.txt")

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


correlationC<-ggplot(cluster, aes(x=global, y=local)) + geom_point()+stat_cor(method = "kendall", cor.coef.name="tau", label.x = 25, label.y = 16,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-trap insertions")+ylab("trap insertions")
correlationR<-ggplot(reference, aes(x=global, y=local)) + geom_point()+stat_cor(method = "kendall", cor.coef.name="tau", label.x = 25, label.y = 12.0,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-reference insertions")+ylab("reference insertions")

####trajectories:

trajectories<-ggplot(output, aes(x = generation, y = tes, group = replicate,col= factor(phase))) +
  geom_line() +geom_vline(xintercept =  gentx,col="black")+scale_color_manual(values = c("green", "yellow", "red"))+theme(legend.position = "None")+ylab("TE abundance")+xlim(0,7500) 



g30<-ggarrange(trajectories,histo,historef, correlationC, correlationR,
             #labels = c("A","B","C","D","E"),
             ncol = 5, nrow = 1)
```

    ## Warning: Removed 7500 row(s) containing missing values (geom_path).

``` r
####50%

tally<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/sitex05-weaker/tally-sitex05-weaker.txt")


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

histo<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("trap insertions")

t<-subset(tally,generation==gentx)
t<-subset(t,type=="ref")

for (sid in unique(t$abundance)) { 
  i <- t$abundance == sid
  a = sum(t$number[i])
  t$sum[i] = a
} 

ts<-unique(subset(t,select=c("abundance","sum")))
ts$rel<-ts$sum/sum(ts$sum)

historef<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("reference insertions")


output<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/sitex05-weaker/output-sitex05-weaker.txt")

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


correlationC<-ggplot(cluster, aes(x=global, y=local)) + geom_point()+stat_cor(method = "kendall", cor.coef.name="tau", label.x = 12.5, label.y = 13,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-trap insertions")+ylab("trap insertions")
correlationR<-ggplot(reference, aes(x=global, y=local)) + geom_point()+stat_cor(method = "kendall", cor.coef.name="tau", label.x = 12.5, label.y = 10.0,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-reference insertions")+ylab("reference insertions")

####trajectories:

trajectories<-ggplot(output, aes(x = generation, y = tes, group = replicate,col= factor(phase))) +
  geom_line() +geom_vline(xintercept =  gentx,col="black")+scale_color_manual(values = c("green", "yellow", "red"))+theme(legend.position = "None")+ylab("TE abundance")+xlim(0,7500) 



g50<-ggarrange(trajectories,histo,historef, correlationC, correlationR,
             #labels = c("A","B","C","D","E"),
             ncol = 5, nrow = 1)
```

    ## Warning: Removed 7500 row(s) containing missing values (geom_path).

``` r
####70%

tally<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/sitex07-weaker/tally-sitex07-weaker.txt")


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

histo<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("trap insertions")

t<-subset(tally,generation==gentx)
t<-subset(t,type=="ref")

for (sid in unique(t$abundance)) { 
  i <- t$abundance == sid
  a = sum(t$number[i])
  t$sum[i] = a
} 

ts<-unique(subset(t,select=c("abundance","sum")))
ts$rel<-ts$sum/sum(ts$sum)

historef<-ggplot(ts, aes(x=abundance, y=rel)) + geom_bar(stat = "identity")+ylab("frequency")+xlab("reference insertions")


output<-read.table("/Users/filipwierzbicki/Desktop/trap_model/analysis/simulations/storm2/sitex07-weaker/output-sitex07-weaker.txt")

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


correlationC<-ggplot(cluster, aes(x=global, y=local)) + geom_point()+stat_cor(method = "kendall", cor.coef.name="tau", label.x = 8, label.y = 12.5,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-trap insertions")+ylab("trap insertions")
correlationR<-ggplot(reference, aes(x=global, y=local)) + geom_point()+stat_cor(method = "kendall", cor.coef.name="tau", label.x = 8, label.y = 5.5,size=3)+ 
  geom_smooth(method='lm', formula= y~x)+xlab("non-reference insertions")+ylab("reference insertions")

####trajectories:

trajectories<-ggplot(output, aes(x = generation, y = tes, group = replicate,col= factor(phase))) +
  geom_line() +geom_vline(xintercept =  gentx,col="black")+scale_color_manual(values = c("green", "yellow", "red"))+theme(legend.position = "None")+ylab("TE abundance")+xlim(0,7500) 



g70<-ggarrange(trajectories,histo,historef, correlationC, correlationR,
             #labels = c("A","B","C","D","E"),
             ncol = 5, nrow = 1)
```

    ## Warning: Removed 7500 row(s) containing missing values (geom_path).

``` r
#############
#combine plots:

g<-ggarrange(g10,g30,g50,g70,
             labels = c("A","B","C","D"),
             ncol = 1, nrow = 4)

plot(g)
```

![](supp_negative-sites_files/figure-gfm/cars-1.png)<!-- -->

``` r
ggsave("/Users/filipwierzbicki/Desktop/trap_model/figures/simulations/simulations_supp-negative-site.pdf",height = 12,width = 12)
ggsave("/Users/filipwierzbicki/Desktop/trap_model/figures/simulations/simulations_supp-negative-site.png",height = 12,width = 12)
```
