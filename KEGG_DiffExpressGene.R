library(ggplot2)

mytheme_violin <-  theme_bw()+#theme(legend.position="none")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"))+
  theme(plot.title=element_text(size=8))+
  # theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))+
  theme(axis.text.x = element_text(colour="black", size=8))+
  theme(axis.text.y = element_text(colour="black", size=8))+
  theme(axis.title=element_text(size=8))+
  #theme(plot.title=element_text(size=rel(1.5), lineheight=.9))+
  theme(axis.ticks=element_line(colour="black",size=0.5))+
  theme(legend.title=element_text(size=8))+
  theme(legend.text=element_text(size=8))

library(tidyverse)
library(DESeq2)
library(dplyr)

##### Identify differentially expressed genes

path = '/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/4-DiffExpressGene/ColData_and_CountDataframe/'
sample = 'PAB1'

mycounts <- read.table(paste0(path,sample,'_CountData.txt'), header = TRUE)
mycounts <- na.omit(mycounts)
mycounts <- data.frame(mycounts, row.names = NULL)
head(mycounts)

metadata <- read.table(paste0(path,sample,'_ColData.txt'), header = TRUE)
metadata <- data.frame(metadata, row.names = NULL)
metadata

mycounts <- as.data.frame(mycounts) # Convert the data type to a dataframe

metadata <- as.data.frame(metadata)
head(mycounts)
metadata
class(mycounts)
class(metadata)

names(mycounts)[-1] # Check if the column names in the count data match the 'id' in the colData.

metadata$id
names(mycounts)[-1]==metadata$id
all(names(mycounts)[-1]==metadata$id)

dds <- DESeqDataSetFromMatrix(countData=mycounts, 
                              colData=metadata, 
                              design=~treatment, 
                              tidy=TRUE)
dds <- DESeq(dds)
res <- results(dds, tidy=TRUE)
res <- as_tibble(res)
write.table(res,paste0(path,sample,'_DEseq2.txt'),quote = FALSE)


Fold_pval <- read.table(paste0(path,sample,'_DEseq2.txt'))
Fold_pval <- na.omit(Fold_pval)

for (n in 1:length(Fold_pval$row)){
  if (Fold_pval$padj[n] <= 0.05 & Fold_pval$log2FoldChange[n] < 0){
    Fold_pval$group[n] = 'Down'
  } else if (Fold_pval$padj[n] <= 0.05 & Fold_pval$log2FoldChange[n] > 0){
    Fold_pval$group[n] = 'Up'
  }else {
    Fold_pval$group[n] = 'Not_significant'
  }
}

for (n in 1:length(Fold_pval$row)){
  if (Fold_pval$group[n] == 'Down' | Fold_pval$group[n] == 'Up'){
    Fold_pval$type[n] = 'Significant'
  } else {
    Fold_pval$type[n] = 'Not_significant'
  }
}

Fold_pval_TypeList <- Fold_pval[,c(1,8,9)]
head(Fold_pval_TypeList)

write.table(Fold_pval_TypeList,paste0(path,sample,'_DEseq2_TypeList.txt'),sep = '\t', quote=FALSE, row.names=FALSE)


##### KEGG enrichment analysis of differentially expressed genes

library(AnnotationDbi)
library(org.Sc.sgd.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)

sample <- 'PAB1'
path = 'E:/Q-lab/WJ-Ribo-seq/ColData_and_CountDataframe/'
diff <- read.table(paste0(path,sample,'_DEseq2_TypeList_unique.txt'),header = TRUE)
diff_Up <- diff[(diff$group == 'Up')&(diff$Unique == 'Yes'),]
diff_Down <- diff[(diff$group == 'Down')&(diff$Unique == 'Yes'),]

write.table(diff_Up$row,file = paste0(path,sample,'_UpGene_unique.txt'), row.names = F, col.names = F, quote = F,sep = '\t')
write.table(diff_Down$row,file = paste0(path,sample,'_DownGene_unique.txt'), row.names = F, col.names = F, quote = F,sep = '\t')



result_Up <- enrichKEGG(gene = diff_Up$row,
                        organism     = "sce",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        minGSSize = 1)

result_Up <- as.data.frame(result_Up)
write.table(result_Up,file = paste0(path,sample,'_KEGG_Up_unique.txt'), row.names = F, quote = F,sep = '\t')

result_Down <- enrichKEGG(gene = diff_Down$row,
                          organism     = "sce",
                          pvalueCutoff = 1,
                          qvalueCutoff = 1,
                          minGSSize = 1)

result_Down <- as.data.frame(result_Down)
write.table(result_Down,file = paste0(path,sample,'_KEGG_Down_unique.txt'), row.names = F, quote = F,sep = '\t')


##### KEGG heatmap

path = '/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/4-DiffExpressGene/KEGG_pathway/'

KEGG_FDR_Down <- read.table(paste0(path,'AllSample_KEGG_FDR_Down_0.05.txt'),sep = '\t',header = TRUE)
rownames(KEGG_FDR_Down) <- KEGG_FDR_Down$Term
KEGG_FDR_Down <- KEGG_FDR_Down[, -1]
KEGG_FDR_Down <- -log10(KEGG_FDR_Down)
KEGG_FDR_Down<- KEGG_FDR_Down[,c('eIF1','eIF1A','eIF4A','eIF4B','eIF4E','eIF4G1','PAB1','CAF20','EAP1')]
head(KEGG_FDR_Down)

color_palette <- colorRampPalette(colors = c('white','blue'))(n = 500)
breaks <- seq(0, 5, length.out = 500 + 1)

pheatmap_Down <- pheatmap(KEGG_FDR_Down,
                         cluster_cols = FALSE,
                         cluster_rows = FALSE,
                         main = 'Down-regulated Gene: KEGG pathway',
                         breaks = breaks,
                         color = color_palette,
                         border_color = 'grey')


KEGG_FDR_Up <- read.table(paste0(path,'AllSample_KEGG_FDR_Up_0.05.txt'),sep = '\t',header = TRUE)
rownames(KEGG_FDR_Up) <- KEGG_FDR_Up$Term
KEGG_FDR_Up <- KEGG_FDR_Up[, -1]
KEGG_FDR_Up <- -log10(KEGG_FDR_Up)
KEGG_FDR_Up<- KEGG_FDR_Up[,c('eIF1','eIF1A','eIF4A','eIF4B','eIF4E','eIF4G1','PAB1','CAF20','EAP1')]
head(KEGG_FDR_Up)

color_palette <- colorRampPalette(colors = c('white','red'))(n = 500)
breaks <- seq(0, 6, length.out = 500 + 1)

pheatmap_Up <- pheatmap(KEGG_FDR_Up,
                          cluster_cols = FALSE,
                          cluster_rows = FALSE,
                          main = 'Up-regulated Gene: KEGG pathway',
                          breaks = breaks,
                          color = color_palette,
                          border_color = 'grey')

pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/3-UMAP_KEGG_RNAandRPFfoldchangeCorrelation/HeatMap-KEGG-FDR-Down-0.05.pdf",width=6,height=5,useDingbats = F, colormodel = "rgb")
print(pheatmap_Down,vp=viewport(1,1,x=.5,y=.5))
dev.off()

pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/3-UMAP_KEGG_RNAandRPFfoldchangeCorrelation/HeatMap-KEGG-FDR-Up-0.05.pdf",width=6,height=5,useDingbats = F, colormodel = "rgb")
print(pheatmap_Up,vp=viewport(1,1,x=.5,y=.5))
dev.off()
#dev.new()