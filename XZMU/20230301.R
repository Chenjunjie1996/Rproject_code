library(Seurat)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library("ggsci")

rds = readRDS("rename.rds")
meta = rds@meta.data
UMAPPlot(rds)
# 基因表达变化

Idents(rds) <- "sample"
rds = subset(rds,idents=c('WYP0','WYP1','WYP2'))
# DYC-0-X1 DYC-1-X1 DYC-2-X1 WYP-0-X2 WYP-1-X1 WYP-2-X2 
meta = rds@meta.data
#Idents(rds) <- "sample"
#markers = FindAllMarkers(rds)

GC10 = FindMarkers(object = rds,ident.1 = 'GC1', ident.2 = 'GC0')
GC20 = FindMarkers(object = rds,ident.1 = 'GC2', ident.2 = 'GC0')
GC21 = FindMarkers(object = rds,ident.1 = 'GC2', ident.2 = 'GC1')

GC10$gene = rownames(GC10)
GC20$gene = rownames(GC20)
GC21$gene = rownames(GC21)

GC10 = dplyr::select(GC10, -3,-4)
GC20 = dplyr::select(GC20, -3,-4)
GC21 = dplyr::select(GC21, -3,-4)

colnames(GC10) <- c("p_val_10", "avg_log2fc_10","p_val_adj_10","gene")
colnames(GC20) <- c("p_val_20", "avg_log2fc_20","p_val_adj_20","gene")
colnames(GC21) <- c("p_val_21", "avg_log2fc_21","p_val_adj_21","gene")

merge1 = merge(GC10, GC21, by='gene')
merge2 = merge(merge1, GC20, by='gene')

write.table(merge2,file="C:/Users/admin/Desktop/daily work/GITHUB/徐州医科大施明/20230220subsetT_harmony/20230301/GC_gene_exp_change.xls",sep="\t",quote=F, row.names = F)

# 一直上调
VlnPlot(rds, pt.size=0, features = c('GNLY','NKG7','S100A8','S100A9','CCL5','LYZ',"GZMH","CCL4","PLAAT4"))
# 先上后下
VlnPlot(rds, pt.size=0, features = c('HBB','HBA2','HBA1','ACTB','COTL1','ISG15','MCM5','FABP5','IFI6'))
# 先下后上
VlnPlot(rds, pt.size=0, features = c('MALAT1','XIST','RPS12','GBP5','RPS27','RPS29','CD247','STAT1','PLAC8'))
# 一直下
VlnPlot(rds, pt.size=0, features = c('TPI1','ENO1','GAPDH','PKM','H4C3','TUBA1B', "LDHA", "TUBB", "STMN1"))



# GSVA 
logFCcutoff=0
adjPvalueCutoff=1
library(GSVA)
library(limma)
library(GSEABase)

rds = readRDS('C:/Users/admin/Desktop/daily work/GITHUB/徐州医科大施明/20230220subsetT_harmony/T_pos_1.2/rename.rds')
Idents(rds) <- "sample"
#rds <- subset(rds,idents = c("DYC-0-X1","DYC-1-X1","DYC-2-X1"))
rds <- subset(rds,idents = c("WYP0","WYP1","WYP2"))

meta = rds@meta.data
write.table(meta,file="tmpmeta_car.xls",sep ='\t',quote=F,row.names=T)
meta=read.table("./tmpmeta_car.xls")

#文件中有-时 需要xls另存为csv文件
meta=read.csv("./tmpmeta_car.csv",header = T, row.names = 1)

rds@meta.data = meta

Idents(rds) <- "sample"
rt <- AverageExpression(rds)
write.table(rt,"averageexpression.xls",sep="\t",quote=F)

rt=read.table("./averageexpression.xls",sep="\t",header=T,check.names=F,row.names=1)
rt=as.matrix(rt)
dimnames=list(rownames(rt),colnames(rt))
mat=matrix(as.numeric(as.matrix(rt)),nrow=nrow(rt),dimnames=dimnames)
mat=avereps(rt)
mat=normalizeBetweenArrays(mat)
colnames(mat) <- c("WYP0", "WYP1","WYP2")
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
#rt <- rt[c(1,2,3,4,6,5)]
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
deg = FindMarkers(object = rds,ident.1 = 'GC2', ident.2 = 'GC1' ,min.pct = 0.01,logfc.threshold = 0.01)
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

gseaplot2(egmt,geneSetID = c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_KRAS_SIGNALING_UP','HALLMARK_INFLAMMATORY_RESPONSE'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))

gseaplot2(egmt,geneSetID = c('HALLMARK_E2F_TARGETS','HALLMARK_MYC_TARGETS_V1','HALLMARK_G2M_CHECKPOINT'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))

# 特殊
gseaplot2(egmt,geneSetID = c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_WNT_BETA_CATENIN_SIGNALING'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue'))


# 输出表达矩阵
rds = readRDS('C:/Users/admin/Desktop/daily work/GITHUB/徐州医科大施明/20230220subsetT_harmony/T_pos_1.2/rename.rds')
meta = rds@meta.data
UMAPPlot(rds)

Idents(rds) <- "sample"
rds1 = subset(rds,idents=c('WYP2'))

exprMatrix <- as.matrix(GetAssayData(rds1, slot='counts'))
write.table(exprMatrix,file = 'C:/Users/admin/Desktop/daily work/GITHUB/徐州医科大施明/20230220subsetT_harmony/20230313/WYP2_matrix.tsv',col.names=NA,sep="\t",quote =FALSE)
