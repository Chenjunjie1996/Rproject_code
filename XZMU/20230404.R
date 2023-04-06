library(Seurat)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library("ggsci")

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
#deg = FindMarkers(object = rds,ident.1 = 'WYP2', ident.2 = 'WYP0')
deg = FindMarkers(object = rds,ident.1 = 'GC2', ident.2 = 'GC0' ,min.pct = 0.1,logfc.threshold = 0.1)
deg1 = FindMarkers(object = rds,ident.1 = 'GC2', ident.2 = 'GC1' ,min.pct = 0.1,logfc.threshold = 0.1)
deg2 = FindMarkers(object = rds,ident.1 = 'GC1', ident.2 = 'GC0' ,min.pct = 0.1,logfc.threshold = 0.1)
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
head(egmt)
egmt@result 
gsea_results_df <- egmt@result 
rownames(gsea_results_df)

m1=as.data.frame(egmt@result)

gseaplot2(egmt,geneSetID = c('EXHAUSTION_RELATED_GENES'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('lightskyblue'))





rds = readRDS('C:/Users/admin/Desktop/daily work/GITHUB/徐州医科大施明/20230220subsetT_harmony/T_pos_1.2/rename.rds')
Idents(rds) <- "celltype"
rds <- subset(rds,idents = c("CD8_Teff")) 
UMAPPlot(rds,group.by=("celltype"),label = TRUE,label.size = 5, pt.size=1)+ scale_color_npg()

Idents(rds) <- "sample"
#rds <- subset(rds,idents = c("DYC-0-X1","DYC-1-X1","DYC-2-X1"))
rds <- subset(rds,idents = c("DYC0","DYC1","DYC2"))
marker <- c("HAVCR2",
            "RGS16",
            "LAYN",
            "SRGAP3",
            "DUSP4",
            "CSF1",
            "TNFRSF9",
            "LYST",
            "TNFRSF18",
            "NDFIP2",
            "SQLE",
            "ID3",
            "SOX4",
            "CD9",
            "PHLDA1",
            "CCL3",
            "CCL4",
            "KLRC1",
            "KLRD1",
            "KLRB1",
            "KLRC2",
            "CDK6",
            "PLS3",
            "AFAP1L2",
            "CTSW",
            "IL2RA",
            "AHI1",
            "RBPJ",
            "GZMB",
            "GNLY",
            "IL7R",
            "TCF7",
            "SELL",
            "KLF2",
            "TUBA1B",
            "TOP2A",
            "PCNA",
            "CD8A",
            "CD3E"
)
marker = rev(marker)
DotPlot(rds, features = marker)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5, size=12, color="black",face="bold"),
        axis.text.y=element_text(size=12, color="black",face="bold"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))

# +scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))


# heatmap
library(pheatmap)

rds = readRDS('C:/Users/admin/Desktop/daily work/GITHUB/徐州医科大施明/20230220subsetT_harmony/T_pos_1.2/rename.rds')
Idents(rds) <- "celltype"
rds <- subset(rds,idents = c("CD8_Teff")) 
Idents(rds) <- "sample"
#rds <- subset(rds,idents = c("DYC-0-X1","DYC-1-X1","DYC-2-X1"))
rds <- subset(rds,idents = c("GC0","GC1","GC2"))


Idents(rds) <- "sample"
rt <- AverageExpression(rds)
target = c('ACP5', 'ADGRG1', 'AFAP1L2', 'AKAP5', 'ANXA5', 'APOBEC3C', 'APOBEC3G', 'ATP6V1C2', 'BATF', 'BST2', 'CCL3', 'CCL4', 'CCL4L1', 'CCND2', 'CCR1', 'CD27', 'CD27-AS1', 'CD2BP2', 'CD38', 'CD63', 'CD7', 'CD82', 'CDK2AP2', 'CHST12', 'CKS2', 'COTL1', 'COX5A', 'CREM', 'CSF1', 'CTLA4', 'CTSD', 'CTSW', 'CXCL13', 'CXCR6', 'DDIT4', 'DNPH1', 'DUSP4', 'DYNLL1', 'ENTPD1', 'ENTPD1-AS1', 'FABP5', 'FASLG', 'FKBP1A', 'FKBP1A-SDCBP2', 'FUT8', 'GALM', 'GPR25', 'GSTO1', 'GZMB', 'GZMH', 'HAVCR2', 'HLA-DMA', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DRB6', 'HMGN1', 'HMGN3', 'ICOS', 'ID3', 'IDH2', 'IFI27L2', 'IFI35', 'IFI6', 'IFNG', 'IGFLR1', 'IL2RB', 'ISG15', 'ITGAE', 'ITM2A', 'KRT81', 'KRT86', 'LAG3', 'LAYN', 'LINC00299', 'LYST', 'MIR155', 'MIR155HG', 'MIR3917', 'MIR4632', 'MIR497HG', 'MS4A6A', 'MTHFD1', 'MTHFD2', 'MYO1E', 'MYO7A', 'NAB1', 'NDFIP2', 'PARK7', 'PDCD1', 'PDIA6', 'PHLDA1', 'PKM', 'PRDM1', 'PRDX3', 'PRDX5', 'PRF1', 'PRKAR1A', 'PSMB3', 'PSMC3', 'PSMD4', 'PSMD8', 'PTTG1', 'RAB27A', 'RALGDS', 'RANBP1', 'RBPJ', 'RGS1', 'RGS2', 'SAMSN1', 'SARDH', 'SIRPG', 'SIT1', 'SNAP47', 'SNRPB', 'SNX9', 'STAT3', 'STMN1', 'STRA13', 'SYNGR2', 'TIGIT', 'TNFRSF18', 'TNFRSF1B', 'TNFRSF9', 'TNFSF4', 'TNIP3', 'TOX', 'TPI1', 'TRAFD1', 'UBE2F', 'UBE2F-SCLY', 'UBE2L6', 'VAPA', 'VCAM1', 'WARS', 'YARS')

write.table(rt,"averageexpression.xls",sep="\t",quote=F)
rt=read.table("./averageexpression.xls",sep="\t",header=T,check.names=F,row.names=1)
rt = rt[rownames(rt) %in% target, ]
colnames(rt) <- gsub("RNA.", "", colnames(rt))
rt = rt[apply(rt, 1, function(x) sd(x)!=0),]
pheatmap::pheatmap(as.matrix(rt),cluster_rows = T,cluster_cols = T,scale = "row",angle_col =0,  fontsize = 12,
                   color = colorRampPalette(c("navy","white","firebrick3"))(100),
                   border_color = "NA")
dev.off()


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
head(egmt)
egmt@result 
gsea_results_df <- egmt@result 
rownames(gsea_results_df)

m1=as.data.frame(egmt@result)

gseaplot2(egmt,geneSetID = c('EXHAUSTION_RELATED_GENES'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('lightskyblue'))
