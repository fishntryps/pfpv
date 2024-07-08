
#########################################################################
### just a quick visualization of read-depths using ggplot2 and dplyr ###
#########################################################################

### run time is a couple of seconds

### load packages in R 4.2.2

library(dplyr)
library(ggplot2)

### this table has mean read-depth for filter-passing and non-passing samples

### the table also shows percent bases with >5 reads

x <- read.table("sample_seq_info.txt", header=T)

x %>% filter(passed %in% "YES") %>% ggplot(aes(x=spec, fill=spec, y=avdepth)) + geom_boxplot(outlier.color = NA, alpha=0.5) + 
  geom_jitter(shape=21, aes(bg=spec), stroke=0.001) + 
  facet_grid(.~passed) + 
  ylim(0,180) + 
  scale_fill_manual(values=c("#D41159","#1A85FF")) + 
  scale_color_manual(values=c("#D41159","#1A85FF")) + labs(y="Mean read-depth\n") +
  theme_minimal() + theme(axis.text = element_text(size=12)) + 
  theme(axis.title = element_text(size=14)) 


median(x$avdepth[x$spec %in% "Pf" & x$passed %in% "YES"])
median(x$avdepth[x$spec %in% "Pv" & x$passed %in% "YES"])


### this table has mean read-depth plus COI calls for filter-passing samples

x <- read.table("coil_calls.txt", header=F)
colnames(x) <- c("sample", "coi", "avdepth")
x$spec <- gsub("_.*", "", x$sample)

x %>% ggplot(aes(x=coi, group=coi, y=avdepth)) + geom_boxplot(outlier.color = NA, alpha=0.5) + 
  geom_jitter(shape=21, aes(bg=spec), stroke=0.001) + 
  facet_grid(.~spec) + 
  ylim(0,180) + 
  scale_fill_manual(values=c("#D41159","#1A85FF")) + 
  scale_color_manual(values=c("#D41159","#1A85FF")) + labs(y="Mean read-depth\n") +
  theme_minimal() + theme(axis.text = element_text(size=12)) + 
  theme(axis.title = element_text(size=14)) 
