suppressMessages({library(tidyverse)
  library(Seurat)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library("ggsci")
  library("gridExtra")})
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tibble)
library(GSVA)
library(limma)
library(GSEABase)
library(presto)
library(DESeq2)
library(clusterProfiler)
library(pheatmap)
library(Seurat)
library(org.Hs.eg.db)
library(enrichplot)

# GSEA
rds = readRDS('C:/Users/admin/Desktop/daily work/GITHUB/徐州医科大施明/20230220subsetT_harmony/T_pos_1.2/rename.rds')
Idents(rds) <- "sample"
rds <- subset(rds,idents = c("GC0","GC1","GC2"))
exp <- rds@assays$RNA@data
dim(exp)
exp[1:4,1:4]
table(rds$sample)
Idents(rds) <- 'sample'
# , min.pct = 0,logfc.threshold = 0 
# DYC-0-X1 DYC-1-X1 DYC-2-X1 WYP-0-X2 WYP-1-X1 WYP-2-X2
deg = FindMarkers(object = rds,ident.1 = 'GC2', ident.2 = 'GC0', min.pct = 0.1,logfc.threshold = 0.1)
deg1 = FindMarkers(object = rds,ident.1 = 'GC1', ident.2 = 'GC0' ,min.pct = 0.1,logfc.threshold = 0.1)
deg2 = FindMarkers(object = rds,ident.1 = 'GC2', ident.2 = 'GC1' ,min.pct = 0.01,logfc.threshold = 0.01)

deg = deg1
deg = deg2
#deg%>%
#  filter(cluster %in% c('GC_0_ZL'))->GC0
#deg%>%
#  filter(cluster %in% c('GC_1_ZL'))->GC1
#deg%>%
#  filter(cluster %in% c('GC_2_ZL'))->GC2

geneList=deg$avg_log2FC
names(geneList)= toupper(rownames(deg))
geneList=sort(geneList,decreasing = T)
head(geneList)

gmtfile ='./h.all.v7.2.symbols1.gmt'
geneset <- read.gmt( gmtfile )
length(unique(geneset$term))
#egmt <- GSEA(geneList, TERM2GENE=geneset)
egmt <- GSEA(geneList, TERM2GENE=geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
#head(egmt)
#egmt@result 
gsea_results_df <- egmt@result 
#rownames(gsea_results_df)

m1=as.data.frame(egmt@result)

gseaplot2(egmt,geneSetID = c('HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_COMPLEMENT'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))

gseaplot2(egmt,geneSetID = c('HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_MYC_TARGETS_V1','HALLMARK_INTERFERON_GAMMA_RESPONSE'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))

gseaplot2(egmt,geneSetID = c('EXHAUSTION_RELATED_GENES'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('lightskyblue'))

write.table(m1,file="gc2vs1.tsv",sep ='\t',quote=F,row.names=T)


# 特殊
gseaplot2(egmt,geneSetID = c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_WNT_BETA_CATENIN_SIGNALING'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue'))





rds = readRDS('C:/Users/admin/Desktop/daily work/GITHUB/徐州医科大施明/20230220subsetT_harmony/T_pos_1.2/rename.rds')
Idents(rds) <- "celltype"
rds <- subset(rds,idents = c("CD8_Teff")) 
Idents(rds) <- "sample"
rds <- subset(rds,idents = c("WYP0","WYP1","WYP2"))
exp <- rds@assays$RNA@data
dim(exp)
exp[1:4,1:4]
table(rds$sample)
Idents(rds) <- 'sample'
# , min.pct = 0,logfc.threshold = 0 
# DYC-0-X1 DYC-1-X1 DYC-2-X1 WYP-0-X2 WYP-1-X1 WYP-2-X2
#deg = FindMarkers(object = rds,ident.1 = 'WYP2', ident.2 = 'WYP0')
deg = FindMarkers(object = rds,ident.1 = 'WYP2', ident.2 = 'WYP0' ,min.pct = 0,logfc.threshold = 0)
deg1 = FindMarkers(object = rds,ident.1 = 'WYP2', ident.2 = 'WYP1' ,min.pct = 0,logfc.threshold = 0)
deg2 = FindMarkers(object = rds,ident.1 = 'WYP1', ident.2 = 'WYP0' ,min.pct = 0,logfc.threshold = 0)
deg = deg1
deg = deg2
#deg = FindAllMarkers(rds)
head(deg)

#deg%>%
#  filter(cluster %in% c('GC_0_ZL'))->GC0
#deg%>%
#  filter(cluster %in% c('GC_1_ZL'))->GC1
#deg%>%
#  filter(cluster %in% c('GC_2_ZL'))->GC2

geneList=deg$avg_log2FC
names(geneList)= toupper(rownames(deg))
geneList=sort(geneList,decreasing = T)
head(geneList)

gmtfile ='./h.all.v7.2.symbols1.gmt'
geneset <- read.gmt( gmtfile )
length(unique(geneset$term))
#egmt <- GSEA(geneList, TERM2GENE=geneset)
egmt <- GSEA(geneList, TERM2GENE=geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
#head(egmt)
#egmt@result 
gsea_results_df <- egmt@result 
#rownames(gsea_results_df)

m1=as.data.frame(egmt@result)
write.table(m1,file="wyp1vs0.tsv",sep ='\t',quote=F,row.names=T)

gseaplot2(egmt,geneSetID = c('EXHAUSTION_RELATED_GENES'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('lightskyblue'))
