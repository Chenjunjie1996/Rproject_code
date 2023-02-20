library(Seurat)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)

rds = readRDS("demo.diff_PRO.rds")
meta = rds@meta.data

outP = stringr::str_glue("cluster.png")
png(outP, height=1200, width=1200)
UMAPPlot(rds,group.by=("raw_cluster"),label = TRUE,label.size = 6)
dev.off()

outP = stringr::str_glue("sample.png")
png(outP, height=1200, width=1200)
UMAPPlot(rds,group.by=("sample"),label = TRUE,label.size = 6)
dev.off()

Idents(rds) <- "raw_cluster"
allmarkers <- FindAllMarkers(rds)
allmarkers <- subset(allmarkers,avg_log2FC>0) 
write.table(allmarkers,file="diffgenes.xls",sep='\t',quote=F,row.names=T)

target_genes= c('ANXA1', 'ANXA2', 'BATF', 'BCL6', 'BHLHE40', 'BTLA', 'CAPG', 'CCL5', 'CCR4', 'CCR6', 'CCR7', 'CCR8', 'CD160', 'CD200', 'CD27', 'CD28', 'CD4', 'CD44', 'CD6', 'CD69', 'CD8A', 'CD8B', 'CDC25B', 'CTLA4', 'CTSH', 'CTSW', 'CX3CR1', 'CXCR3', 'CXCR4', 'CXCR5', 'CXCR6', 'ENTPD1', 'EOMES', 'FCGR3A', 'FGFBP2', 'FOXP3', 'FURIN', 'GNLY', 'GPR183', 'GZMA', 'GZMB', 'GZMH', 'GZMK', 'HAVCR2', 'ICA1', 'ICAM2', 'ICOS', 'IFNG', 'IGFLR1', 'IKZF2', 'IL10', 'IL10RA', 'IL17A', 'IL21R', 'IL23R', 'IL2RA', 'IL6ST', 'IL7R', 'ITGAE', 'KIR2DL4', 'KLRB1', 'KLRC1', 'KLRC2', 'KLRC3', 'KLRG1', 'LAYN', 'LEF1', 'MAGEH1', 'MIR155HG', 'MYADM', 'NCR3', 'NKG7', 'NR4A1', 'NR4A2', 'NR4A3', 'PDCD1', 'PRF1', 'PTGER2', 'PTGER4', 'RGS1', 'RORA', 'RORC', 'RTKN2', 'RUNX3', 'S1PR1', 'S1PR4', 'S1PR5', 'SELL', 'SLC4A10', 'TBX21', 'TCF7', 'TMIGD2', 'TNFRSF18', 'TNFRSF9', 'TOX', 'TOX2', 'XCL1', 'XCL2', 'ZBTB16')

Idents(rds) <- "raw_cluster"
df = data.frame(genes = target_genes)
for (i in as.numeric(levels(rds))){
  tmp = c()
  rds1 <- subset(rds,idents = c(i))
  for (j in target_genes){
    print(j)
    percent = sum(GetAssayData(object = rds1, slot = "data")[j,]>0)/nrow(rds1@meta.data)
    percent = round(percent, 4)
    percent = scales::percent(percent, 0.01)
    print(percent)
    tmp <- c(tmp, percent)
  }
  df = cbind(df, tmp)
  colnames(df)[i+1] = paste('cluster',i,sep='_')
}
write.table(df,file="gene_exp_percent.tsv",sep='\t',quote=F,row.names=T)

# 文献中的marker
# CD8+TN
outP = stringr::str_glue("CD8_TN.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("CCR7","LEF1","SELL","TCF7","CD27", "CD28","S1PR1"))
dev.off()
# CD8+TCM
outP = stringr::str_glue("CD8_TCM.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("CCR7","SELL","IL7R","CD27","CD28","PRF1","GZMA","CCL5","GPR183","S1PR1"))
dev.off()
# CD8+TEMRA/TEFF
outP = stringr::str_glue("CD8_TEMRA_TEFF.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("KLRG1","CX3CR1","FCGR3A","FGFBP2","PRF1","GZMH","TBX21","EOMES","S1PR1","S1PR5"))
dev.off()
# CD8+TEM
outP = stringr::str_glue("CD8_TEM.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("GZMK","CXCR4","CXCR3","CD44"))
dev.off()
# CD8+TRM
outP = stringr::str_glue("CD8_TRM.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("CD6","XCL1","XCL2","MYADM","CAPG","RORA","NR4A1","NR4A2","NR4A3","CD69","ITGAE"))
dev.off()
# CD8+IEL
outP = stringr::str_glue("CD8_IEL.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("CD160","KIR2DL4","TMIGD2","KLRC1","KLRC2","KLRC3","NR4A1","NR4A2","NR4A3","IKZF2","ENTPD1","CD69","ITGAE"))
dev.off()
# CD8+TEx
outP = stringr::str_glue("CD8_TEx.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("HAVCR2","CXCL13","PDCD1","LAYN","TOX","IFNG","GZMB","MIR155HG","TNFRSF9","ITGAE"))
dev.off()
# MAIT
outP = stringr::str_glue("MAIT.png")
png(outP, height=1080, width=1080)
FeaturePlot(rds,features =c("SLC4A10","KLRB1","ZBTB16","NCR3","RORC","RORA"))
dev.off()


# CD4+TN
outP = stringr::str_glue("CD4_TN.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("CCR7","LEF1","SELL","TCF7","CD27","CD28","S1PR1"))
dev.off()
# CD4+Blood-TCM
outP = stringr::str_glue("CD4Blood_TCM.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("CCR7","SELL","PTGER2","ICAM2","ANXA1","ANXA2","S1PR1"))
dev.off()
# CD4+TEMRA/TEFF
outP = stringr::str_glue("CD4_TEMRA_TEFF.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("KLRG1", "CX3CR1","NKG7","PRF1","GNLY","GZMH","TBX21","CTSW","S1PR1","S1PR5"))
dev.off()
# CD4+Normal-Tcm
outP = stringr::str_glue("CD4_NormalTcm.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("CCR7","TCF7","RGS1","CD69"))
dev.off()
# CD4+TRM
outP = stringr::str_glue("CD4_TRM.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("CD69","KLRB1","PTGER4","IL7R","CXCR6","NR4A1","NR4A2","NR4A3","MYADM"))
dev.off()
# TFh
outP = stringr::str_glue("TFh.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("CXCR5","BCL6","ICA1","TOX","TOX2","IL6ST","MAGEH1","BTLA","ICOS","PDCD1", "CD200"))
dev.off()
# CD4+TEM/TH1-like
outP = stringr::str_glue("CD4_TEM_TH1like.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("GZMK","GZMA","CCL5","IFNG","RUNX3","EOMES","CXCR3","CXCR4","CD44"))
dev.off()
# TH17
outP = stringr::str_glue("TH17.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("IL23R","RORC","IL17A","FURIN","CTSH","CCR6","KLRB1","CAPG","ITGAE"))
dev.off()
# TH1-like T cells
outP = stringr::str_glue("TH1_like.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("CXCL13","IFNG","CXCR3","BHLHE40","GZMB","PDCD1","HAVCR2","ICOS","IGFLR1","ITGAE"))
dev.off()
# Blood-Treg
outP = stringr::str_glue("Blood_Treg.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("FOXP3","IL2RA","IL10RA","IKZF2","RTKN2","CDC25B","S1PR4"))
dev.off()
# TFR
outP = stringr::str_glue("TFR.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("FOXP3","IL2RA","CXCR5","PDCD1","IL10","CCR4","CD69"))
dev.off()
# Tumour-Treg
outP = stringr::str_glue("Tumour_Treg.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("FOXP3","CCR8","TNFRSF18","LAYN","TNFRSF9","IKZF2","RTKN2","CTLA4","BATF","IL21R"))
dev.off()

# CD4/CD8
outP = stringr::str_glue("CD4_CD8.png")
png(outP, height=1080, width=1920)
FeaturePlot(rds,features =c("CD4","CD8A","CD8B"))
dev.off()
