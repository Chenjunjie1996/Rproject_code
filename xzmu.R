# 徐州医科大施明

suppressMessages({library(tidyverse)
  library(Seurat)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library("ggsci")
  library("gridExtra")})
memory.limit(1000000000)

# read rds and add car 
rds = readRDS('GC.diff_PRO.rds')
GC0 = read.table('220309010_404_FJ_filtered_UMI.csv',sep=',',header=T)
GC1 = read.table('220303012CD_FJ_filtered_UMI.csv',sep=',',header=T)
GC2 = read.table('220305023ABC_FJ_filtered_UMI.csv',sep=',',header=T)
GC0$shi<-paste('GC_0_ZL',"_",sep="")
GC0$barcode<-str_c(GC0$shi,GC0$barcode)
GC0_cells <- subset(GC0, sum_UMI!=0)
GC0_cells <- unique(GC0_cells$barcode)
GC1$shi<-paste('GC_1_ZL',"_",sep="")
GC1$barcode<-str_c(GC1$shi,GC1$barcode)
GC1_cells <- subset(GC1, sum_UMI!=0)
GC1_cells <- unique(GC1_cells$barcode)
GC2$shi<-paste('GC_2_ZL',"_",sep="")
GC2$barcode<-str_c(GC2$shi,GC2$barcode)
GC2_cells <- subset(GC2, sum_UMI!=0)
GC2_cells <- unique(GC2_cells$barcode)
cells = c(GC0_cells,GC1_cells,GC2_cells)
meta = rds@meta.data
meta$car = 'CAR-'
meta[cells,'car'] = 'CAR+'
rds@meta.data=meta

# tsneplot
Idents(rds) <- "orig.ident"
rds0 = subset(rds,idents='GC_0_ZL')
rds1 = subset(rds,idents='GC_1_ZL')
rds2 = subset(rds,idents='GC_2_ZL')

TSNEPlot(rds, group.by='celltype')+ scale_color_npg()
TSNEPlot(rds0, group.by='celltype')+ scale_color_npg()
TSNEPlot(rds1, group.by='celltype')+ scale_color_npg()
TSNEPlot(rds2, group.by='celltype')+ scale_color_npg()

TSNEPlot(rds, group.by='car',cols = c("grey", "red"))
TSNEPlot(rds0, group.by='car',cols = c("grey", "red"))
TSNEPlot(rds1, group.by='car',cols = c("grey", "red"))
TSNEPlot(rds2, group.by='car',cols = c("grey", "red"))

# barplot
PP <- ggplot(data=meta, aes(x= celltype, fill= car))+
  geom_bar(width=0.6, alpha=0.8 , colour="slategray",size=1 )+
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 16,color="black"),axis.text.y = element_text(size = 14,color="black"),
        legend.text = element_text(size=15)) +scale_fill_manual(values=c("lightskyblue","#dc143c"))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
print(PP)

# vlnplot
VlnPlot(rds, pt.size=0,features = c("TNFRSF9","PRDX1","HLA-DQA1","PKM","HSPD1",'BIRC3'),cols=c('grey','red'))
VlnPlot(rds, pt.size=0,features = c("HBA2","HBA1","HBB","S100A8","S100A9","FCER1G"),cols=c('grey','red'))

install.packages('dplyr')
BiocManager::install("tidyverse")
BiocManager::install("fgsea")
BiocManager::install("dplyr")
BiocManager::install("ggplot2")
BiocManager::install("tibble")
BiocManager::install("GSVA")
BiocManager::install("limma")
# GSEA
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

exp <- rds@assays$RNA@data
dim(exp)
exp[1:4,1:4]
table(rds$car)
Idents(rds) <- 'car'
deg = FindMarkers(object = rds,ident.1 = 'CAR+',min.pct = 0.1,logfc.threshold = 0.25)
head(deg)
geneList=deg$avg_log2FC
names(geneList)= toupper(rownames(deg))
geneList=sort(geneList,decreasing = T)
head(geneList)

gmtfile ='./h.all.v7.2.symbols1.gmt.txt'
# 31120 个基因集
#GSEA分析
library(GSEABase) # BiocManager::install('GSEABase')
geneset <- read.gmt( gmtfile )  
length(unique(geneset$term))
egmt <- GSEA(geneList, TERM2GENE=geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
head(egmt)
egmt@result 
gsea_results_df <- egmt@result 
rownames(gsea_results_df)


m1=as.data.frame(egmt@result)

gseaplot2(egmt,geneSetID = c('HALLMARK_MYC_TARGETS_V1','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))

#下调
gseaplot2(egmt,geneSetID = c('HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_TNFA_SIGNALING_VIA_NFKB'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))


# GSVA
rm(list = ls())
logFCcutoff=0
adjPvalueCutoff=1
library(GSVA)
library(limma)
library(GSEABase)

meta = rds@meta.data
#write.table(meta,file="meta_car.xls",sep = '\t',quote=F,row.names=T)
meta=read.table("./meta_car.xls")
rds@meta.data = meta

Idents(rds) <- "celltype"
rt <- AverageExpression(rds)
write.table(rt,"averageexpression.xls",sep="\t",quote=F)

rt=read.table("./averageexpression.xls",sep="\t",header=T,check.names=F,row.names=1)
rt=as.matrix(rt)
dimnames=list(rownames(rt),colnames(rt))
mat=matrix(as.numeric(as.matrix(rt)),nrow=nrow(rt),dimnames=dimnames)
mat=avereps(rt)
mat=normalizeBetweenArrays(mat)
gmtFile = getGmt('./h.all.v7.2.symbols1.gmt.txt')

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

library(pheatmap)

tiff(file="heatmap.tiff",
     width = 80,            #图片的宽度
     height =140,            #图片的高度
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

# 拟时序
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

##
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

rm(list = ls())


#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(GetAssayData(rds, slot='counts')), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = rds@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
HSMM<-monocle_cds
HSMM@phenoData@data
meta = HSMM@phenoData@data
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 3 )
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),num_cells_expressed >= 10))
print(head(pData(HSMM)))

#Clustering cells without marker genes 

disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
# Plots the percentage of variance explained by the each component based on PCA from the normalized expression data using the same procedure used in reduceDimension function.
# HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(HSMM, return_all = F) # norm_method='log'

HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 10,reduction_method = 'tSNE', verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 2)
HSMM <- clusterCells(HSMM, num_clusters = 5)
plot_cell_clusters(HSMM, 1, 2)

HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 2,
                        reduction_method = 'tSNE',
                        residualModelFormulaStr = "~Size_Factor + num_genes_expressed",
                        verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 5)
plot_cell_clusters(HSMM, 1, 2, color = "Cluster")  #

diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~percent.mt")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) ## 不要也写0.1 ，而是要写0.01。
#ordering_genes <- row.names(subset(diff_test_res,num_cells_expressed>300))

HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)

HSMM <- reduceDimension(HSMM, max_components = 2,method = 'DDRTree')
HSMM <- orderCells(HSMM)

plot_cell_trajectory(HSMM, color_by = "celltype",cell_size = 2 ) + 
  guides(colour = guide_legend(override.aes = list(size=4)))+
  theme(legend.text = element_text(size=15))+ scale_color_lancet()


plot_cell_trajectory(HSMM, color_by = "State",cell_size = 2)+ 
  guides(colour = guide_legend(override.aes = list(size=4)))+ scale_color_npg()

plot_cell_trajectory(HSMM, color_by = "Pseudotime",cell_size = 2) +
  theme(legend.text = element_text(size=10),legend.background = element_rect(colour='black', size=1),
        legend.key.size = unit(50, "pt"))


plot_cell_trajectory(HSMM, color_by = "celltype",cell_size = 4) + 
  geom_point(aes(size = Pseudotime , color = celltype)) + scale_color_npg() + 
  theme(legend.position = "right")


library(pheatmap)
p=plot_pseudotime_heatmap(HSMM[ordering_genes,],num_clusters = 6,cores = 1,return_heatmap=T,show_rownames = T)
print(p)


# Cellphone DB 
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(SeuratData)


write.table(as.matrix(rds@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(rds@meta.data), rds@meta.data[,'celltype', drop=F])
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

netf<- "count_network.txt"
mynet <- read.delim(paste0("./count_network.txt"), check.names = FALSE)
table(mynet$count)
mynet %>% filter(count>0) -> mynet  # 有零会报错
head(mynet)
net<- graph_from_data_frame(mynet,)
plot(net)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
            "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493",
            "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")


karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net, order =
                             order(membership(karate_groups)))  # 设置网络布局
E(net)$width  <- E(net)$count/10  # 边点权重（粗细）
plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7) 

net2 <- net  # 复制一份备用

for (i in 1: length(unique(mynet$SOURCE)) ){
  E(net)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
}  # 这波操作谁有更好的解决方案？ 

plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7) 


plot(net, edge.arrow.size=.1, 
     edge.curved=0.2, # 只是调了这个参数
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=1) 


dev.off()

length(unique(mynet$SOURCE)) # 查看需要绘制多少张图，以方便布局
par(mfrow=c(2,4), mar=c(.3,.3,.3,.3))

for (i in 1: length(unique(mynet$SOURCE)) ){
  net1<-net2
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
  
  plot(net1, edge.arrow.size=0.2, 
       edge.curved=0.4,
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1) 
  
}



dev.off()
length(unique(mynet$SOURCE))
par(mfrow=c(3,3), mar=c(.3,.3,.3,.3))

for (i in 1: length(unique(mynet$SOURCE)) ){
  net1<-net2
  E(net1)$count <- ""
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  <- E(net2)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  # 故技重施
  
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
  
  plot(net1, edge.arrow.size=.1, 
       edge.curved=0.2,
       edge.label = E(net1)$count, # 绘制边的权重
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1.2
  ) 
  
}
dev.off()



# 这些基因list很有意思啊，建议保存
mypvals <- read.delim(paste0(rds,"pvalues.txt"), check.names = FALSE)
mymeans <- read.delim(paste0(rds,"means.txt"), check.names = FALSE)



# 这些基因list很有意思啊，建议保存
chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR", mymeans$interacting_pair,value = T)
chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR", mymeans$interacting_pair,value = T)
th1 <- grep("IL2|IL12|IL18|IL27|IFNG|IL10|TNF$|TNF |LTA|LTB|STAT1|CCR5|CXCR3|IL12RB1|IFNGR1|TBX21|STAT4", 
            mymeans$interacting_pair,value = T)
th2 <- grep("IL4|IL5|IL25|IL10|IL13|AREG|STAT6|GATA3|IL4R", 
            mymeans$interacting_pair,value = T)
th17 <- grep("IL21|IL22|IL24|IL26|IL17A|IL17A|IL17F|IL17RA|IL10|RORC|RORA|STAT3|CCR4|CCR6|IL23RA|TGFB", 
             mymeans$interacting_pair,value = T)
treg <- grep("IL35|IL10|FOXP3|IL2RA|TGFB", mymeans$interacting_pair,value = T)
costimulatory <- grep("CD86|CD80|CD48|LILRB2|LILRB4|TNF|CD2|ICAM|SLAM|LT[AB]|NECTIN2|CD40|CD70|CD27|CD28|CD58|TSLP|PVR|CD44|CD55|CD[1-9]", 
                      mymeans$interacting_pair,value = T)
coinhibitory <- grep("SIRP|CD47|ICOS|TIGIT|CTLA4|PDCD1|CD274|LAG3|HAVCR|VSIR", 
                     mymeans$interacting_pair,value = T)
niche <- grep("CSF", mymeans$interacting_pair,value = T)


mymeans %>% dplyr::filter(interacting_pair %in% costimulatory)%>%
  dplyr::select("interacting_pair",starts_with("NK"),ends_with("NK"))  %>%  
  reshape2::melt() -> meansdf

colnames(meansdf)<- c("interacting_pair","CC","means")

mypvals %>% dplyr::filter(interacting_pair %in% costimulatory)%>%
  dplyr::select("interacting_pair",starts_with("NK"),ends_with("NK"))%>%  
  reshape2::melt()-> pvalsdf

colnames(pvalsdf)<- c("interacting_pair","CC","pvals")
pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")

summary((filter(pldf,means >1))$means)

pldf%>% filter(means >1) %>% 
  ggplot(aes(CC.x,interacting_pair.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
  scale_size_continuous(range = c(1,3))+
  scale_color_gradient2(high="red",mid = "yellow",low ="darkblue",midpoint = 25  )+ theme_bw()+ 
  theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 0.8)) 
