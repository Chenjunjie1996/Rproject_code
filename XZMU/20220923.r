memory.limit(1000000000)

# read rds and add car 
rds = readRDS('GC_1.2.rds')

# test
rds1 = readRDS('GC_1.2.rds')
Idents(rds1) <- "celltype"
rds1 <- subset(rds1,idents = c("NaiveT","NKT","CD8Tex","CD8Teff"))
# 1. 基因
rds = readRDS('GC_1.2_T.rds')

Idents(rds) <- "celltype"
NaiveT = subset(rds, idents=c('NaiveT'))
NKT = subset(rds, idents=c('NKT'))
CD8Tex = subset(rds, idents=c('CD8Tex'))
CD8Teff = subset(rds, idents=c('CD8Teff'))

Idents(NKT) <- "orig.ident"
rt <- AverageExpression(NKT)
write.table(rt,"./20220923/1diff_gene/NKT/averageexpression_nkt.txt",sep="\t",quote=F)

Idents(rds) <- "orig.ident"
# 一直上调
VlnPlot(NKT, pt.size=0, features = c('MALAT1','TMSB4X','CD247','RPS27','CCL5','CD52'))
# 先上后下
VlnPlot(NaiveT, pt.size=0, features = c('B2M','RPS18','RPS19','RPL27A','RPL13A','FTH1'))
# 先下后上
VlnPlot(NaiveT, pt.size=0, features = c('TMSB4X','PFN1','CD74','PPIA','HNRNPA2B1','MYL6'))
# 一直下
VlnPlot(NKT, pt.size=0, features = c('NKG7','GAPDH','ACTG1','TUBA1B','LGALS1','GZMB'))


Idents(NaiveT) <- "car"
markers = FindMarkers(NaiveT, ident.1 = 'CAR+')
Idents(NaiveT) <- "orig.ident"
VlnPlot(NaiveT, pt.size=0, features = c('TNFRSF9','HLA-DQA1','MIR155HG','BIRC3','NFKBIA','HLA-DRA'))

Idents(NKT) <- "car"
markers = FindMarkers(NKT, ident.1 = 'CAR+')
Idents(NKT) <- "orig.ident"
VlnPlot(NKT, pt.size=0, features = c('CD27','TRAC','GZMK','CD8B','CD247','CD52'))

Idents(CD8Tex) <- "car"
markers = FindMarkers(CD8Tex, ident.1 = 'CAR+')
Idents(CD8Tex) <- "orig.ident"
VlnPlot(CD8Tex, pt.size=0, features = c('TNFRSF9','FABP5','NME1','MRPL13','PRDX1','FTL'))

Idents(CD8Teff) <- "car"
markers = FindMarkers(CD8Teff, ident.1 = 'CAR+')
Idents(CD8Teff) <- "orig.ident"
VlnPlot(CD8Teff, pt.size=0, features = c('TNFRSF9','IFI6','TUBA1B','ISG15','LAG3','PDIA6'))

meta$celltype_car <-str_c(meta$celltype, "_", meta$car)
rds@meta.data = meta
Idents(rds) <- "celltype_car"
markers = FindAllMarkers(rds)
markers_top5 = markers%>%
  group_by(cluster)%>%
  top_n(n = 5,wt = avg_log2FC) 
DoHeatmap(object = rds,features = unique(markers_top5$gene),size=4)
?DoHeatmap()
# 4. 拟时序
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

#Extract data, phenotype data, and feature data from the SeuratObject
rds = readRDS('GC_1.2_T.rds')
data <- as.sparse(rds@assays$RNA@counts)
pd <- new('AnnotatedDataFrame', data = rds@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              expressionFamily = negbinomial.size())
HSMM<-monocle_cds
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
# HSMM <- detectGenes(HSMM, min_expr = 3 ) # 很多教程不用

##使用clusters差异表达基因
diff.wilcox = FindAllMarkers(rds)
all.markers = diff.wilcox %>% subset(p_val<0.05)
diff.genes <- subset(all.markers,p_val_adj<0.01)$gene
HSMM <- setOrderingFilter(HSMM, diff.genes)
p1 <- plot_ordering_genes(HSMM)
##使用seurat选择的高变基因
var.genes <- VariableFeatures(rds)
HSMM <- setOrderingFilter(HSMM, var.genes)
p2 <- plot_ordering_genes(HSMM)
##使用monocle选择的高变基因
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
p3 <- plot_ordering_genes(HSMM)
##结果对比
p3
p1|p2|p3

HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
HSMM <- orderCells(HSMM)

##分离CAR+和CAR-
meta = HSMM@phenoData@data
meta$car = 'CAR-'
meta[GC1_cells,'car'] = 'CAR+'
meta = meta %>% drop_na(orig.ident)

#write.table(meta,file="HSMMmeta.xls",sep = '\t',quote=F,row.names=T)
#write.table(meta,file="HSMMmetaT.xls",sep = '\t',quote=F,row.names=T)
meta=read.table("./HSMMmeta.xls")
meta=read.table("./HSMMmetaT.xls")
HSMM@phenoData@data=meta

plot_cell_trajectory(HSMM, color_by = "State")

plot_cell_trajectory(HSMM, color_by = "Pseudotime",cell_size = 1) +
  theme(legend.text = element_text(size=10),legend.background = element_rect(colour='black', size=1),
        legend.key.size = unit(50, "pt"))

plot_cell_trajectory(HSMM, color_by = "celltype",cell_size = 1.5) + 
  guides(colour = guide_legend(override.aes = list(size=4)))+
  theme(legend.text = element_text(size=15))+ scale_color_npg()


# 5. 富集分析
Idents(rds) <- 'car'
rds <- subset(rds,idents = c("CAR+"))
exp <- rds@assays$RNA@data
dim(exp)
exp[1:4,1:4]
table(rds$car)
Idents(rds) <- 'orig.ident'
# , min.pct = 0,logfc.threshold = 0
deg = FindMarkers(object = rds,ident.1 = 'GC_2_ZL', ident.2 = 'GC_0_ZL' ,min.pct = 0.1,logfc.threshold = 0.1)
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

gmtfile ='./h.all.v7.2.symbols.gmt'
geneset <- read.gmt( gmtfile )
length(unique(geneset$term))
#egmt <- GSEA(geneList, TERM2GENE=geneset)
egmt <- GSEA(geneList, TERM2GENE=geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
head(egmt)
egmt@result 
gsea_results_df <- egmt@result 
rownames(gsea_results_df)


m1=as.data.frame(egmt@result)

gseaplot2(egmt,geneSetID = c('HALLMARK_COMPLEMENT','HALLMARK_ALLOGRAFT_REJECTION','HALLMARK_TNFA_SIGNALING_VIA_NFKB'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))

gseaplot2(egmt,geneSetID = c('HALLMARK_MYC_TARGETS_V1','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))

# 特殊
gseaplot2(egmt,geneSetID = c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_ALLOGRAFT_REJECTION'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue'))


# GSVA
logFCcutoff=0
adjPvalueCutoff=1
library(GSVA)
library(limma)
library(GSEABase)

meta = rds@meta.data
write.table(meta,file="tmpmeta_car.xls",sep = '\t',quote=F,row.names=T)
meta=read.table("./tmpmeta_car.xls")
rds@meta.data = meta

Idents(rds) <- "orig.ident"
rt <- AverageExpression(rds)
write.table(rt,"averageexpression.xls",sep="\t",quote=F)

rt=read.table("./averageexpression.xls",sep="\t",header=T,check.names=F,row.names=1)
rt=as.matrix(rt)
dimnames=list(rownames(rt),colnames(rt))
mat=matrix(as.numeric(as.matrix(rt)),nrow=nrow(rt),dimnames=dimnames)
mat=avereps(rt)
mat=normalizeBetweenArrays(mat)
gmtFile = getGmt('./h.all.v7.2.symbols.gmt')

gsvaOut=gsva(mat,
             gmtFile,
             min.sz=10,
             max.sz=500,
             verbose=TRUE,
             parallel.sz=1)
gsvaOut=rbind(id=colnames(gsvaOut),gsvaOut)
write.table(gsvaOut,file="gsvaOut.txt",sep="\t",quote=F,col.names=F)
rt=read.table("gsvaOut.txt",sep="\t",header=T,check.names=F,row.names=1)

fit=lmFit(rt)
fit=eBayes(fit)
all=topTable(fit, number=Inf,adjust.method="holm")
all=rbind(id=colnames(all),all)
write.table(all,file="all.txt",sep="\t",quote=F,col.names=F)

diff <- topTable(fit, number=Inf,
                 p.value=adjPvalueCutoff, adjust="holm", lfc=logFCcutoff)
diffName=row.names(diff)
diff=rbind(id=colnames(diff),diff)
write.table(diff,file="diff.xls",sep="\t",quote=F,col.names=F)

hmExp=rt[diffName,]
hmExp=rbind(id=colnames(hmExp),hmExp)
write.table(hmExp,file="heatmap.txt",sep="\t",quote=F,col.names=F)


rt=read.table("heatmap.txt",sep="\t",header=T,row.names=1,check.names=F)
rt <- rt[c(1,2,9,8,3,4,6,5,7,13,12,15,11,10,14,16)]
library(pheatmap)

tiff(file="heatmap.tiff",
     width = 80,            #ͼƬ?Ŀ???
     height =140,            #ͼƬ?ĸ߶?
     units ="cm",
     compression="lzw",
     bg="white",
     res=500)
fontsize_new = 15
fontsize_new <- as.numeric(fontsize_new)

pheatmap(rt,
         color = colorRampPalette(c( "MistyRose", "white","FireBrick"))(10),
         cluster_cols =F,
         fontsize = fontsize_new,
         fontsize_row=20,
         fontsize_col=20,angle_col = 45,
         cellwidth = 50,
         cellheight = 50,
         legend_labels = 50
)
dev.off()
sessionInfo()
