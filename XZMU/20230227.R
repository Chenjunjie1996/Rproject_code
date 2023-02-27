library(Seurat)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library("ggsci")

rds = readRDS("demo.diff_PRO.rds")
meta = rds@meta.data
UMAPPlot(rds)

Idents(rds) <- "raw_cluster"
new.cluster.ids <- c('CD8_Teff','CD4_Tem','CD4_CD8_Tcm','CD4_Tn','CD4_Tn','CD8_Tcm','CD4_Tcm','CD4_Tcm','CD8_Tcm','CD4_CD8_Tcm','CD8_Teff','CD8_Teff','CD8_Teff','CD8_Tem')
names(new.cluster.ids) <- levels(rds)
rds <- RenameIdents(rds, new.cluster.ids)
rds <- StashIdent(object = rds, save.name = "celltype")
saveRDS(rds, "rename.rds")

rds = readRDS("rename.rds")
UMAPPlot(rds,group.by=("celltype"),label = TRUE,label.size = 5, pt.size=1)+ scale_color_npg()

# GC DYC WYP
Idents(rds) <- "sample"
rds_tmp = subset(rds,idents = c("WYP0","WYP1","WYP2"))
UMAPPlot(rds_tmp,group.by=("celltype"),label = TRUE,label.size = 5, pt.size=1)+ scale_color_npg()

# 显示数值 geom_text(aes(label=..count..),stat = 'count', size=6, position = position_stack(vjust = 0.8))
ggplot(data=rds_tmp@meta.data, aes(x= sample, fill= celltype))+
  geom_bar(width=0.6, alpha=0.8 , colour="slategray",size=0.2 )+ 
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 16,color="black"),axis.text.y = element_text(size = 14,color="black"),
        legend.text = element_text(size=15)) +   scale_fill_brewer(palette = "Set3")+ scale_y_continuous(expand = c(0,0)) +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

ggplot(data=rds_tmp@meta.data, aes(x= sample, fill= celltype))+
  geom_bar(width=0.6, alpha=0.8 , colour="slategray",size=0.2 , position = 'fill') + labs( y = 'Percent (%)') +
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 16,color="black"),axis.text.y = element_text(size = 14,color="black"),
        legend.text = element_text(size=15)) +   scale_fill_brewer(palette = "Set3")+ scale_y_continuous(expand = c(0,0)) +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))


# count 
Idents(rds) <- "sample"
for (i in c("GC0","GC1","GC2","DYC0","DYC1","DYC2","WYP0","WYP1","WYP2")){
  rds_tmp = subset(rds,idents = c(i))
  print(i)
  print(table(rds_tmp@meta.data$celltype))
  print(dim(rds_tmp@meta.data))
}

