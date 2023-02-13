library(Seurat)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)

rds = readRDS("rename.rds")
meta = rds@meta.data
meta = meta %>%
  mutate(time = if_else(sample %in% c("DYC0","WYP0","GC0"),  
                       true = "time0",
                       false = if_else(
                         sample %in% c("DYC1","WYP1","GC1"),
                         true = "time1",
                         false= "time2"
                       )))
rds@meta.data = meta

UMAPPlot(rds, group.by="celltype", label = TRUE,label.size = 5)
UMAPPlot(rds, group.by="time", label = TRUE,label.size = 5)

Idents(rds) <- "sample" 
rds1 <- subset(rds,idents = c("GC0","GC1","GC2"))

#BARPLOT

PP <- ggplot(data=rds1@meta.data, aes(x= celltype, fill= sample))+
  geom_bar(width=0.6, alpha=0.8 , colour="slategray",size=1 )+ 
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 16,color="black"),axis.text.y = element_text(size = 14,color="black"),
        legend.text = element_text(size=15)) +scale_fill_manual(values=c("lightskyblue","lightgreen", "#dc143c"))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  geom_text(aes(label=..count..),stat = 'count', size=6, position = position_stack(vjust = 0.8))
print(PP)


