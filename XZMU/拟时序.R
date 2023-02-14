library('Seurat')
library('ggplot2')
library("ggsci")
library('tidyverse')
library('monocle')
options(warn=-1)

rds  =readRDS("/SGRNJ06/randd/USER/cjj/celedev/XZMU/20230213/rename.rds")

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

##使用monocle选择的高变基因
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
p3 <- plot_ordering_genes(HSMM)

HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
HSMM <- orderCells(HSMM)

HSMM@phenoData@data

# 时间轴
plot_cell_trajectory(HSMM, color_by = "Pseudotime",cell_size = 1) +
  theme(legend.text = element_text(size=10),legend.background = element_rect(colour='black', size=1),
        legend.key.size = unit(50, "pt"))

ggplot2::ggsave('pseudotime.png', height=10, width=10)

# 有时候（大多数时候），拟时序的方向或是根节点弄错了，还需要手动更改
HSMM=orderCells(HSMM,root_state = 2) 
plot_cell_trajectory(HSMM, color_by = "Pseudotime",cell_size = 1) +
  theme(legend.text = element_text(size=10),legend.background = element_rect(colour='black', size=1),
        legend.key.size = unit(50, "pt"))

plot_cell_trajectory(HSMM, color_by = "State",cell_size = 1.5) + 
  guides(colour = guide_legend(override.aes = list(size=4)))+
  theme(legend.text = element_text(size=15))+ scale_color_npg()

plot_cell_trajectory(HSMM, color_by = "celltype",cell_size = 1.5) + 
  guides(colour = guide_legend(override.aes = list(size=4)))+
  theme(legend.text = element_text(size=15))+ scale_color_npg()

plot_cell_trajectory(HSMM, color_by = "sample",cell_size = 1) + 
  guides(colour = guide_legend(override.aes = list(size=4)))+
  theme(legend.text = element_text(size=15))+ scale_color_npg()


##分离CAR+和CAR-
meta = HSMM@phenoData@data
write.table(meta,file="/SGRNJ06/randd/USER/cjj/celedev/XZMU/20230103/monocle/wyp//HSMMmeta.xls",sep = '\t',quote=F,row.names=T)

#文件中有-时 需要xls另存为csv文件
meta=read.csv("/SGRNJ06/randd/USER/cjj/celedev/XZMU/20230103/monocle/wyp/HSMMmeta.csv",header = T, row.names = 1)
HSMM@phenoData@data=meta

meta=read.table("/SGRNJ06/randd/USER/cjj/celedev/XZMU/20230103/monocle/wyp/HSMMmeta.xls")
HSMM@phenoData@data=meta