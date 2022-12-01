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

{#~ read rds and add car  ####
  rds = readRDS('S22.rds')
  umi = read.table('C017T1-010-22FS-FJ_filtered_UMI.csv',sep=',',header=T)
  meta = rds@meta.data
  umi$shi<-paste('C017T1-010-22FS',"_",sep="")
  umi$barcode<-str_c(umi$shi,umi$barcode)
  umi_cells <- subset(umi, sum_UMI!=0)
  cells <- unique(umi_cells$barcode)
  
  meta$car = 'CAR-'
  meta$barcode = rownames(meta)
  meta = meta %>%
    mutate(car = if_else(barcode %in% cells,
                         true = "CAR+",
                         false = "CAR-"))
  meta = dplyr::select(meta, -c("barcode"))
  meta <- meta %>% drop_na(nCount_RNA)
  table(meta$car == 'CAR+')
  rds@meta.data=meta
  
  new.cluster.ids <- c('TCells','TCells','LowQuality','LowQuality','TCells','TCells','TCells','TCells','Monocytes','TCells','Erythrocytes')
  names(new.cluster.ids) <- levels(rds)
  rds <- RenameIdents(rds, new.cluster.ids)
  rds <- StashIdent(object = rds, save.name = "celltype")
  
  Idents(rds) <- "celltype"
  rds <- subset(rds,idents = c("TCells","Monocytes","Erythrocytes"))
  UMAPPlot(rds)
  saveRDS(rds, "S22.rds")  
}
  
{#~ umap and bar plot  ####
  rds = readRDS('S22.rds')
  meta= rds@meta.data
  UMAPPlot(rds, group.by='celltype')+ scale_color_npg()
  UMAPPlot(rds, group.by='car',cols = c("grey", "red"))
  
  
  PP <- ggplot(data=meta, aes(x= celltype, fill= car))+
    geom_bar(width=0.6, alpha=0.8 , colour="slategray",size=1 )+
    theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),
          axis.text.x = element_text(size = 16,color="black"),axis.text.y = element_text(size = 14,color="black"),
          legend.text = element_text(size=15)) +scale_fill_manual(values=c("lightskyblue","#dc143c"))+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
  print(PP)
  
  
  
} 

{#~ gsva ####
  logFCcutoff=0
  adjPvalueCutoff=1
  library(GSVA)
  library(limma)
  library(GSEABase)
  
  meta = rds@meta.data
  write.table(meta,file="tmpmeta_car.xls",sep = '\t',quote=F,row.names=T)
  meta=read.table("./tmpmeta_car.xls")
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
  rt <- rt[c(1,5,2,3,4)]
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
}

{#~ gsea ####
  
  # 单样本 CAR+ VS CAR-
  rds = readRDS("S22.rds")

  exp <- rds@assays$RNA@data
  dim(exp)
  exp[1:4,1:4]
  table(rds$car)
  Idents(rds) <- 'car'
  deg = FindMarkers(object = rds,ident.1 = 'CAR+',min.pct = 0,logfc.threshold = 0)
  head(deg)
  geneList=deg$avg_log2FC
  names(geneList)= toupper(rownames(deg))
  geneList=sort(geneList,decreasing = T)
  head(geneList)
  
  gmtfile ='./h.all.v7.2.symbols.gmt'
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
  
  gseaplot2(egmt,geneSetID = c('HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_MYC_TARGETS_V1','HALLMARK_PROTEIN_SECRETION'),
            pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))
  
  gseaplot2(egmt,geneSetID = c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_IL6_JAK_STAT3_SIGNALING','HALLMARK_HEME_METABOLISM'),
            pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))

}

{#~ 拟时序 ####
  library('monocle')
  options(warn=-1)
  rds = readRDS('S22.rds')
  Idents(rds) <- "celltype"
  rds <- subset(rds,idents = c("TCells"))
  
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
  p3
  
  HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
  HSMM <- orderCells(HSMM)
  # 有时候（大多数时候），拟时序的方向或是根节点弄错了，还需要手动更改
  # HSMM=orderCells(HSMM,root_state = 4) 
  
  plot_cell_trajectory(HSMM, color_by = "Pseudotime",cell_size = 1) +
    theme(legend.text = element_text(size=10),legend.background = element_rect(colour='black', size=1),
          legend.key.size = unit(50, "pt"))
  
  
  plot_cell_trajectory(HSMM, color_by = "State",cell_size = 1.5) + 
    guides(colour = guide_legend(override.aes = list(size=4)))+
    theme(legend.text = element_text(size=15))+ scale_color_npg()
  
  ##分离CAR+和CAR-
  meta = HSMM@phenoData@data
  write.table(meta,file="HSMMmeta.xls",sep = '\t',quote=F,row.names=T)
  meta=read.table("HSMMmeta.xls")
  HSMM@phenoData@data=meta
  
  plot_cell_trajectory(HSMM, color_by = "celltype",cell_size = 1) + 
    scale_color_manual(values=c("Green","purple","Gold","grey","red")) +
    guides(colour = guide_legend(override.aes = list(size=4)))+
    theme(legend.text = element_text(size=15))
  
  plot_cell_trajectory(HSMM, color_by = "car",cell_size = 1.5) + 
    guides(colour = guide_legend(override.aes = list(size=4)))+
    theme(legend.text = element_text(size=15))+ scale_color_npg()
  
  plot_cell_trajectory(HSMM, color_by = "car",cell_size = 1) + 
    scale_color_manual(values=c("grey","red")) +
    guides(colour = guide_legend(override.aes = list(size=4)))+
    theme(legend.text = element_text(size=15))
}

{#~ cellphone db  ####
  library(psych)
  library(qgraph)
  library(igraph)
  library(tidyverse)
  library(ggplot2)
  library(Seurat)
  library(SeuratData)
  
  rds = readRDS("S22.rds")
  meta = rds@meta.data
  
  write.table(meta,file="meta42.xls",sep = '\t',quote=F,row.names=T)
  meta=read.table("./meta42.xls")
  rds@meta.data = meta
  
  write.table(as.matrix(rds@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
  meta_data <- cbind(rownames(rds@meta.data), rds@meta.data[,'celltype', drop=F])
  meta_data <- as.matrix(meta_data)
  meta_data[is.na(meta_data)] = "Unkown" #  ϸ???????в?????NA
  write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
  
  netf<- "count_network.txt"
  mynet <- read.delim(paste0("./count_network.txt"), check.names = FALSE)
  table(mynet$count)
  mynet %>% filter(count>0) -> mynet  # ?????ᱨ??
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
                               order(membership(karate_groups)))  # ???????粼??
  E(net)$width  <- E(net)$count/10  # ?ߵ?Ȩ?أ???ϸ??
  plot(net, edge.arrow.size=.1, 
       edge.curved=0,
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=.7) 
  
  net2 <- net  # ????һ?ݱ???
  
  for (i in 1: length(unique(mynet$SOURCE)) ){
    E(net)[map(unique(mynet$SOURCE),function(x) {
      get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
    })%>% unlist()]$color <- allcolour[i]
  }  # ?Ⲩ????˭?и??õĽ????????? 
  
  plot(net, edge.arrow.size=.1, 
       edge.curved=0,
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=.7) 
  
  
  plot(net, edge.arrow.size=.1, 
       edge.curved=0.2, # ֻ?ǵ???????????
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1) 
  
  
  dev.off()
  
  length(unique(mynet$SOURCE)) # ?鿴??Ҫ???ƶ?????ͼ???Է??㲼??
  par(mfrow=c(2,2), mar=c(.3,.3,.3,.3))
  
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
         vertex.label.cex=1.2) 
    
  }
  
  
  
  dev.off()
  length(unique(mynet$SOURCE))
  par(mfrow=c(2,2), mar=c(.3,.3,.3,.3))
  
  for (i in 1: length(unique(mynet$SOURCE)) ){
    net1<-net2
    E(net1)$count <- ""
    E(net1)[map(unique(mynet$SOURCE),function(x) {
      get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
    })%>% unlist()]$count  <- E(net2)[map(unique(mynet$SOURCE),function(x) {
      get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
    })%>% unlist()]$count  # ?ʼ???ʩ
    
    E(net1)[map(unique(mynet$SOURCE),function(x) {
      get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
    })%>% unlist()]$color <- allcolour[i]
    
    plot(net1, edge.arrow.size=.1, 
         edge.curved=0.2,
         edge.label = E(net1)$count, # ???Ʊߵ?Ȩ??
         vertex.color=allcolour,
         vertex.frame.color="#555555",
         vertex.label.color="black",
         layout = coords,
         vertex.label.cex=1.2
    ) 
    
  }
  dev.off()
}


{#~ gene ####
  rds = readRDS("S21.rds")

  Idents(rds) <- 'car'
  deg = FindAllMarkers(object = rds)
  
  my_levels <- c("CAR-","CAR+")
  rds@active.ident <- factor(x = rds@active.ident, levels = my_levels)
  VlnPlot(rds, pt.size=0,features = c("GNLY","NFKBIA","CD27","GZMB","TRBC2","LGALS1"),cols=c('grey','red'))
  VlnPlot(rds, pt.size=0,features = c("GNLY","TRBC2","GZMK","CMC1","KLRC1","CCL5"),cols=c('grey','red'))
  #write.table(deg,file="C:/Users/admin/Desktop/daily work/GITHUB/徐州医科大施明/20221109/gene/单样本差异表达基因/gc2_diff_genes.xls",sep='\t',quote=F,row.names=T)
  
  # heatmap
  rds = readRDS("Tcells.rds")
  meta = rds@meta.data
  meta$celltype_car <-str_c(meta$cluster, "_", meta$car)
  rds@meta.data = meta
  Idents(rds) <- "celltype_car"
  markers = FindAllMarkers(rds)
  markers_top5 = markers%>%
    group_by(cluster)%>%
    top_n(n = 5,wt = avg_log2FC) 
  rds_scale = NormalizeData(rds)
  rds_scale = ScaleData(rds_scale)
  DoHeatmap(object = rds_scale,features = unique(markers_top5$gene),size=4)
}