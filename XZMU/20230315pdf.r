library(Seurat)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library("ggsci")

rds = readRDS("rename.rds")
UMAPPlot(rds,group.by=("celltype"),label = TRUE,label.size = 5, pt.size=1)+ scale_color_npg()

# GC DYC WYP
Idents(rds) <- "sample"
rds_tmp = subset(rds,idents = c("DYC0","DYC1","DYC2"))
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



rds = readRDS("rename.rds")
Idents(rds) <- "sample"
rds = subset(rds,idents=c('WYP0','WYP1','WYP2'))
# DYC-0-X1 DYC-1-X1 DYC-2-X1 WYP-0-X2 WYP-1-X1 WYP-2-X2 
# 一直上调
VlnPlot(rds, pt.size=0, features = c('GNLY','NKG7','S100A8','S100A9','CCL5','LYZ',"GZMH","CCL4","PLAAT4"))
# 先上后下
VlnPlot(rds, pt.size=0, features = c('HBB','HBA2','HBA1','ACTB','COTL1','ISG15','MCM5','FABP5','IFI6'))
# 先下后上
VlnPlot(rds, pt.size=0, features = c('MALAT1','XIST','RPS12','GBP5','RPS27','RPS29','CD247','STAT1','PLAC8'))
# 一直下
VlnPlot(rds, pt.size=0, features = c('TPI1','ENO1','GAPDH','PKM','H4C3','TUBA1B', "LDHA", "TUBB", "STMN1"))

rds = readRDS("rename.rds")
Idents(rds) <- "sample"
rds = subset(rds,idents=c('DYC0','DYC1','DYC2'))
# 一直上调
VlnPlot(rds, pt.size=0, features = c('S100A8','CCL5','NKG7','GZMK','S100A12','CST7',"TXNIP","PLEK","PRF1"))
# 先上后下
VlnPlot(rds, pt.size=0, features = c('HBB','HBA2','HBA1','HMGB2','PPIA','LYZ','NASP','RPSA','ID2'))
# 先下后上
VlnPlot(rds, pt.size=0, features = c('SH2D2A','JUNB','JUN','IER2','NFKBIA','SAMSN1','HSPA5','EZR','CD74'))
# 一直下
VlnPlot(rds, pt.size=0, features = c('TUBA1B','TUBB','PRDX1','H4C3','TPI1','UBE2C', "ENO1", "TXN", "TUBB4B"))


rds = readRDS("rename.rds")
Idents(rds) <- "sample"
rds = subset(rds,idents=c('GC0','GC1','GC2'))
# 一直上调
VlnPlot(rds, pt.size=0, features = c('CD247','HBB','CCL5','HBA2','GZMK','S100A8',"HBA1","S100A9","MALAT1"))
# 先上后下
VlnPlot(rds, pt.size=0, features = c('IFI6','GZMB','IFITM3','ISG15','FCER1G','PRF1','SPON2','LY6E','IFITM1'))
# 先下后上
VlnPlot(rds, pt.size=0, features = c('HSPA8','RPL7A','RPL3','RPL7','RPS25','RPS23','RPS12','RPS3A','RPL19'))
# 一直下
VlnPlot(rds, pt.size=0, features = c('TUBA1B','TUBB','LDHA','ENO1','H2AZ1','TOP2A', "STMN1", "MKI67", "HSP90AB1"))




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
rds <- subset(rds,idents = c("WYP0","WYP1","WYP2"))
exp <- rds@assays$RNA@data
dim(exp)
exp[1:4,1:4]
table(rds$sample)
Idents(rds) <- 'sample'
# , min.pct = 0,logfc.threshold = 0 
# DYC-0-X1 DYC-1-X1 DYC-2-X1 WYP-0-X2 WYP-1-X1 WYP-2-X2
deg = FindMarkers(object = rds,ident.1 = 'WYP2', ident.2 = 'WYP0')
deg = FindMarkers(object = rds,ident.1 = 'WYP2', ident.2 = 'WYP1' ,min.pct = 0.01,logfc.threshold = 0.01)
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

gseaplot2(egmt,geneSetID = c('HALLMARK_ALLOGRAFT_REJECTION','HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_KRAS_SIGNALING_UP'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))

gseaplot2(egmt,geneSetID = c('HALLMARK_E2F_TARGETS','HALLMARK_MYC_TARGETS_V1','HALLMARK_G2M_CHECKPOINT'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))

# 特殊
gseaplot2(egmt,geneSetID = c('HALLMARK_COMPLEMENT','HALLMARK_ALLOGRAFT_REJECTION'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue'))


