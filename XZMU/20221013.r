library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)

rds = readRDS('GC_1.2_T.rds')
meta = rds@meta.data
Idents(rds) <- "orig.ident"
markers = FindAllMarkers(rds)


GC10 = FindMarkers(object = rds,ident.1 = 'GC_1_ZL', ident.2 = 'GC_0_ZL')
GC20 = FindMarkers(object = rds,ident.1 = 'GC_2_ZL', ident.2 = 'GC_0_ZL')
GC21 = FindMarkers(object = rds,ident.1 = 'GC_2_ZL', ident.2 = 'GC_1_ZL')

GC10$gene = rownames(GC10)
GC20$gene = rownames(GC20)
GC21$gene = rownames(GC21)

GC10 = select(GC10, -3,-4)
GC20 = select(GC20, -3,-4)
GC21 = select(GC21, -3,-4)

colnames(GC10) <- c("p_val_10", "avg_log2fc_10","p_val_adj_10","gene")
colnames(GC20) <- c("p_val_20", "avg_log2fc_20","p_val_adj_20","gene")
colnames(GC21) <- c("p_val_21", "avg_log2fc_21","p_val_adj_21","gene")

merge1 = merge(GC10, GC21, by='gene')
merge2 = merge(merge1, GC20, by='gene')

write.table(merge2,file="20221013/allmarkers.txt",sep="\t",quote=F)

# barplot 
meta=read.table("./HSMMmeta.xls")
meta$State = as.character(meta$State)
meta$state = NA
meta = meta %>%
  mutate(state = if_else(State == 3,
                         true = "false",
                         false = "true")) 


PP <- ggplot(data=meta, aes(x= State, fill= car))+
  geom_bar(width=0.6, alpha=0.8 , colour="slategray",size=1 )+
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 16,color="black"),axis.text.y = element_text(size = 14,color="black"),
        legend.text = element_text(size=15))+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
print(PP)

PP <- ggplot(data=meta, aes(x= State, fill= celltype))+
  geom_bar(width=0.6, alpha=0.8 , colour="slategray",size=1 ) + scale_color_npg()
print(PP)

