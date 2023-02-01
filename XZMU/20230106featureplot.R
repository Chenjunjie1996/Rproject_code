library(Seurat)

raw = readRDS("P22030204_230112_homo_PBMC.diff_PRO.rds")
UMAPPlot(raw)
Idents(raw) <- "cluster"
rds <- subset(raw,idents = c("TCells"))
UMAPPlot(rds)

rds <-  NormalizeData(object = rds)
rds <- FindVariableFeatures(object = rds)
rds <-  ScaleData(object = rds)
# genes.use<- head(HVFInfo(object = rds),2000)
rds <- RunPCA(object=rds,features = VariableFeatures(object = rds))
rds <- FindNeighbors(rds, reduction = "pca", dims = 1:20)
# 0.8 1.2 1.6
rds <- FindClusters(rds,resolution = 0.8, algorithm = 1)
rds <- RunTSNE(object=rds,dims.use=1:20,do.fast=TRUE,check_duplicates = FALSE)
rds <- RunUMAP(rds, reduction = "pca", dims = 1:20)

UMAPPlot(rds,group.by=("cluster"))
UMAPPlot(rds,group.by=("seurat_clusters"),label = TRUE,label.size = 6)

#MARKER 
Idents(rds) <- "seurat_clusters"
allmarkers <- FindAllMarkers(rds)
allmarkers <- subset(allmarkers,avg_log2FC>0) 
write.table(allmarkers,file="label/diffgenes.xls",sep='\t',quote=F,row.names=T)

# AVERAGE EXP
rt <- AverageExpression(rds,idents="seurat_clusters")
write.table(rt,"label/averageexpression.tsv",sep="\t",quote=F,row.names=T)

#  567 cd4naive 11,12 cd8naive 13 cd8Teff
# NKT
FeaturePlot(rds,features =c("FCGR3A","KLRD1","KLRC1","KLRF1","CD3D"))
FeaturePlot(rds1,features =c("FCGR3A","KLRD1","KLRC1","KLRF1","CD3D"))

FeaturePlot(rds1,features =c("CD4","CD8A","CD8B"))
# CD4Naive
FeaturePlot(rds1,features =c("CD4","CCR7","LEF1","SELL"))
# CD8Naive 
FeaturePlot(rds1,features =c("CD8A","CCR7","LEF1","SELL"))
# CD4+ tissue-resident memory T cells
FeaturePlot(rds1,features =c("CD4","ITGAE","ITGA1","CXCR6"))
# CD4+ regulatory T cells
FeaturePlot(rds1,features =c("CD4","FOXP3","CTLA4","IKZF2"))
# T-helper 1 cells
FeaturePlot(rds1,features =c("TBX21","IFNG","TNF,","CXCR3"))
# T-helper 2 cells
FeaturePlot(rds1,features =c("GATA3","IL4","IL5","CCR4","CCR6"))
# T-helper 17 cells
FeaturePlot(rds1,features =c("IL17A","IL17F","RORC","CCR6"))
# Follicular helper T cells
FeaturePlot(rds1,features =c("CD40LG","CXCR5","CD200"))
# CD8+ effector T cells
FeaturePlot(rds1,features =c("CD8A","NKG7","GZMA","GNLY"))
# CD8+ exhuasted T cells
FeaturePlot(rds1,features =c("CD8A","LAG3","CTLA4","PDCD1","HAVCR2"))
# memoryT
FeaturePlot(rds1,features =c("CD27","CD28","CD52"))
# Neutrophils
FeaturePlot(rds,features =c("CSF3R","CXCR2","FCGR3B","CAMP","LCN2"))


# 文献中的marker
# CD8+TN
FeaturePlot(rds,features =c("CCR7","LEF1","SELL","TCF7","CD27", "CD28","S1PR1"))
# CD8+TCM
FeaturePlot(rds,features =c("CCR7","SELL","IL7R","CD27","CD28","PRF1","GZMA","CCL5","GPR183","S1PR1"))
# CD8+TEMRA/TEFF
FeaturePlot(rds,features =c("KLRG1","CX3CR1","FCGR3A","FGFBP2","PRF1","GZMH","TBX21","EOMES","S1PR1","S1PR5"))
# CD8+TEM
FeaturePlot(rds,features =c("GZMK","CXCR4","CXCR3","CD44"))
# CD8+TRM
FeaturePlot(rds,features =c("CD6","XCL1","XCL2","MYADM","CAPG","RORA","NR4A1","NR4A2","NR4A3","CD69","ITGAE"))
# CD8+IEL
FeaturePlot(rds,features =c("CD160","KIR2DL4","TMIGD2","KLRC1","KLRC2","KLRC3","NR4A1","NR4A2","NR4A3","IKZF2","ENTPD1","CD69","ITGAE"))
# CD8+TEx
FeaturePlot(rds,features =c("HAVCR2","CXCL13","PDCD1","LAYN","TOX","IFNG","GZMB","MIR155HG","TNFRSF9","ITGAE"))
# MAIT
FeaturePlot(rds,features =c("SLC4A10","KLRB1","ZBTB16","NCR3","RORC","RORA"))

# CD4+TN
FeaturePlot(rds,features =c("CCR7","LEF1","SELL","TCF7","CD27","CD28","S1PR1"))
# CD4+Blood-TCM
FeaturePlot(rds,features =c("CCR7","SELL","PTGER2","ICAM2","ANXA1","ANXA2","S1PR1"))
# CD4+TEMRA/TEFF
FeaturePlot(rds,features =c("KLRG1", "CX3CR1","NKG7","PRF1","GNLY","GZMH","TBX21","CTSW","S1PR1","S1PR5"))
# CD4+Normal-Tcm
FeaturePlot(rds,features =c("CCR7","TCF7","RGS1","CD69"))
# CD4+TRM
FeaturePlot(rds,features =c("CD69","KLRB1","PTGER4","IL7R","CXCR6","NR4A1","NR4A2","NR4A3","MYADM"))
# TFh
FeaturePlot(rds,features =c("CXCR5","BCL6","ICA1","TOX","TOX2","IL6ST","MAGEH1","BTLA","ICOS","PDCD1", "CD200"))
# CD4+TEM/TH1-like
FeaturePlot(rds,features =c("GZMK","GZMA","CCL5","IFNG","RUNX3","EOMES","CXCR3","CXCR4","CD44"))
# TH17
FeaturePlot(rds,features =c("IL23R","RORC","IL17A","FURIN","CTSH","CCR6","KLRB1","CAPG","ITGAE"))
# TH1-like T cells
FeaturePlot(rds,features =c("CXCL13","IFNG","CXCR3","BHLHE40","GZMB","PDCD1","HAVCR2","ICOS","IGFLR1","ITGAE"))
# Blood-Treg
FeaturePlot(rds,features =c("FOXP3","IL2RA","IL10RA","IKZF2","RTKN2","CDC25B","S1PR4"))
# TFR
FeaturePlot(rds,features =c("FOXP3","IL2RA","CXCR5","PDCD1","IL10","CCR4","CD69"))
# Tumour-Treg
FeaturePlot(rds,features =c("FOXP3","CCR8","TNFRSF18","LAYN","TNFRSF9","IKZF2","RTKN2","CTLA4","BATF","IL21R"))

# CD4/CD8
FeaturePlot(rds,features =c("CD4","CD8A","CD8B"))


Idents(rds) <- "seurat_clusters"
rds1 <- subset(rds,idents = c(0))

cluster0 = vector()
cluster1 = vector()
cluster2 = vector()
cluster3 = vector()
cluster4 = vector()
cluster5 = vector()
cluster6 = vector()
cluster7 = vector()
cluster8 = vector()
cluster9 = vector()
cluster10 = vector()
cluster11 = vector()
cluster12 = vector()
cluster13 = vector()
cluster14 = vector()
cluster15 = vector()
cluster16 = vector()
cluster17 = vector()
cluster18 = vector()
cluster19 = vector()
cluster20 = vector()
cluster21 = vector()
cluster22 = vector()
cluster23 = vector()

target_genes= c('ANXA1', 'ANXA2', 'BATF', 'BCL6', 'BHLHE40', 'BTLA', 'CAPG', 'CCL5', 'CCR4', 'CCR6', 'CCR7', 'CCR8', 'CD160', 'CD200', 'CD27', 'CD28', 'CD4', 'CD44', 'CD6', 'CD69', 'CD8A', 'CD8B', 'CDC25B', 'CTLA4', 'CTSH', 'CTSW', 'CX3CR1', 'CXCR3', 'CXCR4', 'CXCR5', 'CXCR6', 'ENTPD1', 'EOMES', 'FCGR3A', 'FGFBP2', 'FOXP3', 'FURIN', 'GNLY', 'GPR183', 'GZMA', 'GZMB', 'GZMH', 'GZMK', 'HAVCR2', 'ICA1', 'ICAM2', 'ICOS', 'IFNG', 'IGFLR1', 'IKZF2', 'IL10', 'IL10RA', 'IL17A', 'IL21R', 'IL23R', 'IL2RA', 'IL6ST', 'IL7R', 'ITGAE', 'KIR2DL4', 'KLRB1', 'KLRC1', 'KLRC2', 'KLRC3', 'KLRG1', 'LAYN', 'LEF1', 'MAGEH1', 'MIR155HG', 'MYADM', 'NCR3', 'NKG7', 'NR4A1', 'NR4A2', 'NR4A3', 'PDCD1', 'PRF1', 'PTGER2', 'PTGER4', 'RGS1', 'RORA', 'RORC', 'RTKN2', 'RUNX3', 'S1PR1', 'S1PR4', 'S1PR5', 'SELL', 'SLC4A10', 'TBX21', 'TCF7', 'TMIGD2', 'TNFRSF18', 'TNFRSF9', 'TOX', 'TOX2', 'XCL1', 'XCL2', 'ZBTB16')
for (i in target_genes){
  print(i)
  percent = sum(GetAssayData(object = rds1, slot = "data")[i,]>0)/nrow(rds1@meta.data)
  print(percent)
  cluster0 <- c(cluster0,percent)
}

Idents(rds) <- "seurat_clusters"
df = data.frame(genes = target_genes)
for (i in 0:23){
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
  colnames(df)[i+2] = paste('cluster',i,sep='_')
}
write.table(df,file="label/gene_exp_percent.tsv",sep='\t',quote=F,row.names=T)
# check 
for (i in 0:23){
  rds1 <- subset(rds,idents = c(i))
  temp = sum(GetAssayData(object = rds1, slot = "data")["ANXA1",]>0)/nrow(rds1@meta.data)
  print(temp)
}

rds1 <- subset(rds,idents = c(23))
temp = GetAssayData(object = rds1, slot = "data")
temp = as.data.frame(temp)
aa = temp[which(rownames(temp) %in% c("ANXA1")),]
aa = t(aa)
