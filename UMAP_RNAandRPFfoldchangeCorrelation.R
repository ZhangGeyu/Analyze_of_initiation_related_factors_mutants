##### UMAP #####

library(Rmisc)
library(tidyverse)
library(ggsci)
library(scales)
library(Cairo)
library(grid)
options(tibble.width = Inf)
options(stringsAsFactors = F)
options(scipen=200)

# Set plot attributes
library(ggprism)

(p03 <-  ggplot()+theme_bw()
  + theme(text = element_text(colour="black", size=8),plot.title=element_text(colour="black",size=8)
          ,axis.text = element_text(colour="black", size=8),axis.title=element_text(colour="black",size=8)
          ,legend.title=element_text(colour="black",size=8),legend.text=element_text(colour="black",size=8)
          ,strip.text = element_text(colour="black", size=8))  
  + theme(axis.ticks=element_line(colour="black",size=0.47),axis.ticks.length =  unit(.08, "cm"))
  + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()
          ,plot.background = element_blank(),panel.background = element_blank()
          ,strip.background = element_blank(),legend.background = element_blank()
          ,panel.border = element_blank(),axis.line = element_line(colour = "black",size = 0.47)
  )
  # theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))+
  #theme(plot.title=element_text(size=rel(1.5), lineheight=.9))+
)

library(ggrepel)
library(umap)

### All nine mutants

d.0 <- read_delim("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/Ribosome_Density/TotalMutant-RibosomeDensity-FoldChange.txt",delim = "\t")
InitiationFactor_4EBP_df <- d.0[,-1]
InitiationFactor_4EBP_df <- d.0[,c('gene','CAF20','EAP1','eIF4E','eIF4G1','PAB1','eIF1','eIF1A','eIF4A','eIF4B')]
InitiationFactor_4EBP_df <- na.omit(InitiationFactor_4EBP_df)
head(InitiationFactor_4EBP_df)
plot(density(log2(InitiationFactor_4EBP_df$CAF20)))

# Take the logarithm
# UMAP was employed in this context, and given that linear data wasn't necessary, omitting the logarithm should have minimal impact
InitiationFactor_4EBP_df_log2 <- log2(InitiationFactor_4EBP_df[,2:ncol(InitiationFactor_4EBP_df)])

t_InitiationFactor_4EBP_df_log2 <- t(InitiationFactor_4EBP_df_log2) # Transpose the rows and columns of the dataframe
t_InitiationFactor_4EBP_df_log2 <- data.frame(t_InitiationFactor_4EBP_df_log2)
head(t_InitiationFactor_4EBP_df_log2[,1:20])

seeddd <- 60
methh <- "cosine"
#methh <- "euclidean"
mindi <- 0.001

n_nei <- 3 # Determine the number of neighboring points, usually set between 2 and 100.
spread_i <- 0.05

umap_1 <- umap(t_InitiationFactor_4EBP_df_log2,n_neighbors=n_nei, metric=methh, min_dist=mindi, random_state= seeddd, spread = spread_i) #$layout %>% as.data.frame
umap_p1 <- data.frame(umap_1$layout)
umap_p1$id <- rownames(umap_1$layout)
head(umap_p1)
#colors <- c("#E69F00","#E69F00","#E69F00","#56B4E9","#56B4E9","#56B4E9","#009E73","#009E73",
#            "#F0E442","#F0E442","#F0E442","#D55E00","#D55E00","#D55E00","#0072B2","#0072B2",
#            "#CC79A7","#CC79A7","#CC79A7","#999999","#999999","#999999","#000000","#000000","#000000")
colors <- c("#e60000","#4fd73d","#8256e9","#999999","#025e44","#F0E442","#D55E00","#0072B2","#CC79A7")
# umap_p1
InitiationFactor_9 <- p03%+%umap_p1 + aes(x = X1,y = X2,color = factor(id))+ 
  geom_point(size = 1.5) + 
  scale_color_manual(values = colors) + 
  labs(x = "UMAP 1",y="UMAP 2",title = paste(methh,mindi,n_nei,spread_i,seeddd,"Ribosome density foldchange [log2]",sep = "  ")) + 
  geom_text_repel(aes(label =factor(id))) + 
  theme(legend.position = "none")
InitiationFactor_9


### Five mutants
d.0 <- read_delim("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/Ribosome_Density/TotalMutant-RibosomeDensity-FoldChange.txt",delim = "\t")
InitiationFactor_4EBP_df <- d.0[,-1]
InitiationFactor_4EBP_df <- d.0[,c('gene','CAF20','EAP1','eIF4E','eIF4G1','PAB1')]
InitiationFactor_4EBP_df <- na.omit(InitiationFactor_4EBP_df)
#plot(density(log2(InitiationFactor_4EBP_df$CAF20)))

InitiationFactor_4EBP_df_log2 <- log2(InitiationFactor_4EBP_df[,2:ncol(InitiationFactor_4EBP_df)])

t_InitiationFactor_4EBP_df_log2 <- t(InitiationFactor_4EBP_df_log2)
t_InitiationFactor_4EBP_df_log2 <- data.frame(t_InitiationFactor_4EBP_df_log2)
head(t_InitiationFactor_4EBP_df_log2[,1:20])

seeddd <- 30
methh <- "cosine"
#methh <- "euclidean"
mindi <- 0.001

n_nei <- 2
spread_i <- 0.05

umap_1 <- umap(t_InitiationFactor_4EBP_df_log2,n_neighbors=n_nei, metric=methh, min_dist=mindi, random_state= seeddd, spread = spread_i) #$layout %>% as.data.frame
umap_p1 <- data.frame(umap_1$layout)
umap_p1$id <- rownames(umap_1$layout)
head(umap_p1)
colors <- c("#e60000","#4fd73d","#D55E00","#0072B2","#CC79A7")
# umap_p1
InitiationFactor_5 <- p03%+%umap_p1 + aes(x = X1,y = X2,color = factor(id))+ 
  geom_point(size = 1.5) + 
  scale_color_manual(values = colors) + 
  labs(x = "UMAP 1",y="UMAP 2",title = paste(methh,mindi,n_nei,spread_i,seeddd,"Ribosome density foldchange [log2]",sep = "  ")) + 
  geom_text_repel(aes(label =factor(id))) + 
  theme(legend.position = "none")
InitiationFactor_5

### Four mutants
d.0 <- read_delim("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/Ribosome_Density/TotalMutant-RibosomeDensity-FoldChange.txt",delim = "\t")
InitiationFactor_4EBP_df <- d.0[,-1]
InitiationFactor_4EBP_df <- d.0[,c('gene','eIF1','eIF1A','eIF4A','eIF4B')]
InitiationFactor_4EBP_df <- na.omit(InitiationFactor_4EBP_df)
#plot(density(log2(InitiationFactor_4EBP_df$eIF1)))

InitiationFactor_4EBP_df_log2 <- log2(InitiationFactor_4EBP_df[,2:ncol(InitiationFactor_4EBP_df)])

t_InitiationFactor_4EBP_df_log2 <- t(InitiationFactor_4EBP_df_log2)
t_InitiationFactor_4EBP_df_log2 <- data.frame(t_InitiationFactor_4EBP_df_log2)
head(t_InitiationFactor_4EBP_df_log2[,1:20])

seeddd <- 20
methh <- "cosine"
#methh <- "euclidean"
mindi <- 0.001

n_nei <- 2
spread_i <- 0.05

umap_1 <- umap(t_InitiationFactor_4EBP_df_log2,n_neighbors=n_nei, metric=methh, min_dist=mindi, random_state= seeddd, spread = spread_i) #$layout %>% as.data.frame
umap_p1 <- data.frame(umap_1$layout)
umap_p1$id <- rownames(umap_1$layout)
head(umap_p1)
colors <- c("#8256e9","#999999","#025e44","#F0E442")
# umap_p1
InitiationFactor_4 <- p03%+%umap_p1 + aes(x = X1,y = X2,color = factor(id))+ 
  geom_point(size = 1.5) + 
  scale_color_manual(values = colors) + 
  labs(x = "UMAP 1",y="UMAP 2",title = paste(methh,mindi,n_nei,spread_i,seeddd,"Ribosome density foldchange [log2]",sep = "  ")) + 
  geom_text_repel(aes(label =factor(id))) + 
  theme(legend.position = "none")
InitiationFactor_4


pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/3-UMAP_KEGG_RNAandRPFfoldchangeCorrelation/UMAP-RibosomeDensityFold.pdf",width=15,height=15,useDingbats = F)
print(InitiationFactor_9,vp=viewport(.22,.22,x=.25,y=.5))
print(InitiationFactor_5,vp=viewport(.22,.22,x=.5,y=.5))
print(InitiationFactor_4,vp=viewport(.22,.22,x=.75,y=.5))
dev.off()
#dev.new()


##### RNA fldchange ~ RPF foldchange #####

library('lmodel2')

path = '/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/Ribosome_Density/'

RPF_RNA_fold <- read.table(paste0(path,'TotalMutant-RNAandRPF-FoldChange-rbind.txt'),header = TRUE)
RPF_RNA_fold <- RPF_RNA_fold[RPF_RNA_fold$sample != 'HO',]
RPF_RNA_fold$RPF_foldchange <- log2(RPF_RNA_fold$RPF_foldchange)
RPF_RNA_fold$RNA_foldchange <- log2(RPF_RNA_fold$RNA_foldchange)
RPF_RNA_fold$sample <- factor(RPF_RNA_fold$sample, levels = c("CAF20","EAP1","eIF4E","eIF4G1","PAB1","eIF1","eIF1A","eIF4A","eIF4B"))
head(RPF_RNA_fold)

Ex5.res <- lmodel2(RPF_foldchange ~ RNA_foldchange,data=RPF_RNA_fold[RPF_RNA_fold$sample == 'eIF4B',], "interval", "interval", 99)
Ex5.res

Fold_slopes <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/Ribosome_Density/MajorAxis.txt',header = TRUE)
Fold_slopes <- as_tibble(Fold_slopes)
head(Fold_slopes)

RPF_RNA_Fold_corr <- ggplot(data = RPF_RNA_fold, 
                            aes(x = RNA_foldchange, y = RPF_foldchange)) +
  coord_cartesian(xlim = c(-3, 3),
                  ylim = c(-3, 3)) + 
  geom_point(size = 0.2, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = 'red',linetype="dashed") + 
  geom_abline(data = Fold_slopes, aes(intercept = Intercept, slope = slope), 
              color = 'blue') + 
  xlab('RNA abundance foldchange [log2]') + ylab('RPF abundance foldchange [log2]') + 
  mytheme_violin + 
  facet_wrap(~sample, scales = 'free', nrow = 2) + 
  stat_cor(label.x = -3, label.y = 2.5, method = "pearson", size = 3) + 
  geom_text(data = Fold_slopes, aes(label = paste0("Major Axis Slope: ", round(slope, 3))), 
            x = Inf, y = -2.5, hjust = 1, vjust = 0, size = 3) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
RPF_RNA_Fold_corr

RPF_RNA_Fold_corr_noWord <- ggplot(data = RPF_RNA_fold, 
                                   aes(x = RNA_foldchange, y = RPF_foldchange)) +
  coord_cartesian(xlim = c(-3, 3),
                  ylim = c(-3, 3)) + 
  geom_point(size = 0.2, alpha = 0.5) +
  xlab('RNA abundance foldchange [log2]') + ylab('RPF abundance foldchange [log2]') + 
  mytheme_violin + 
  facet_wrap(~sample, scales = 'free',nrow = 2) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
RPF_RNA_Fold_corr_noWord


ggsave(
  filename = "/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/RPFfold_RNAfold_corr_MajorAxisSlope-pearson.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 8, height = 3.5, units = "in", dpi = 4000
  )


library("grid")
pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/RPFfold_RNAfold_corr_MajorAxisSlope-pearson.pdf",width=20,height=10,useDingbats = F)
print(RPF_RNA_Fold_corr,vp=viewport(.6, .5, x=.5, y=.5))
dev.off()
#dev.new()

##### uORF: RNA fldchange ~ RPF foldchange #####


path = '/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/Ribosome_Density/'

RPF_RNA_fold <- read.table(paste0(path,'TotalMutant-RNAandRPF-FoldChange-rbind-uORF.txt'),header = TRUE)
RPF_RNA_fold <- RPF_RNA_fold[RPF_RNA_fold$sample != 'HO',]

RPF_RNA_fold$sample <- factor(RPF_RNA_fold$sample, levels = c("EAP1","CAF20","eIF4E","eIF4G1","PAB1","eIF1","eIF1A","eIF4A","eIF4B"))
RPF_RNA_fold$RPF_foldchange <- log2(RPF_RNA_fold$RPF_foldchange)
RPF_RNA_fold$RNA_foldchange <- log2(RPF_RNA_fold$RNA_foldchange)

RPF_RNA_fold_Yes <- RPF_RNA_fold[RPF_RNA_fold$uORF == 'Yes',]

Ex5.res <- lmodel2(RPF_foldchange ~ RNA_foldchange,data=RPF_RNA_fold_Yes[RPF_RNA_fold_Yes$sample == 'CAF20',], "interval", "interval", 99)
Ex5.res

Fold_slopes_Yes <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/Ribosome_Density/MajorAxis_uORF_Yes.txt',header = TRUE)
Fold_slopes_Yes <- as_tibble(Fold_slopes_Yes)
head(Fold_slopes_Yes)

RPF_RNA_Fold_corr_Yes <- ggplot(data = RPF_RNA_fold_Yes, 
                                aes(x = RNA_foldchange, y = RPF_foldchange)) +
  coord_cartesian(xlim = c(-3, 3),
                  ylim = c(-3, 3)) + 
  geom_point(size = 0.2, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = 'red',linetype="dashed") + 
  geom_abline(data = Fold_slopes_Yes, aes(intercept = Intercept, slope = slope), 
              color = 'blue') + 
  xlab('RNA abundance foldchange [log2]') + ylab('RPF abundance foldchange [log2]') + 
  ggtitle('uORF: Yes') + 
  mytheme_violin + 
  facet_wrap(~sample, scales = 'free',ncol = 5) + 
  stat_cor(label.x = -3, label.y = 2.5, method = "pearson",size = 2.5) + 
  geom_text(data = Fold_slopes_Yes, aes(label = paste0("Major Axis Slope: ", round(slope, 3))), 
            x = Inf, y = -2.5, hjust = 1, vjust = 0, size = 2.5) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
RPF_RNA_Fold_corr_Yes

RPF_RNA_Fold_corr_Yes_noWord <- ggplot(data = RPF_RNA_fold_Yes, 
                                aes(x = RNA_foldchange, y = RPF_foldchange)) +
  coord_cartesian(xlim = c(-3, 3),
                  ylim = c(-3, 3)) + 
  geom_point(size = 0.2, alpha = 0.5) +
  xlab('RNA abundance foldchange [log2]') + ylab('RPF abundance foldchange [log2]') + 
  ggtitle('uORF: Yes') + 
  mytheme_violin + 
  facet_wrap(~sample, scales = 'free',ncol = 5) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
RPF_RNA_Fold_corr_Yes_noWord

ggsave(
  filename = "/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/RPFfold_RNAfold_corr_MajorAxisSlope_uORFyes-pearson.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 8, height = 3.5, units = "in", dpi = 4000
)

RPF_RNA_fold_No <- RPF_RNA_fold[RPF_RNA_fold$uORF == 'No',]

Ex5.res <- lmodel2(RPF_foldchange ~ RNA_foldchange,data=RPF_RNA_fold_No[RPF_RNA_fold_No$sample == 'eIF4B',], "interval", "interval", 99)
Ex5.res

Fold_slopes_No <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/Ribosome_Density/MajorAxis_uORF_No.txt',header = TRUE)
Fold_slopes_No <- as_tibble(Fold_slopes_No)
head(Fold_slopes_No)

RPF_RNA_Fold_corr_No <- ggplot(data = RPF_RNA_fold_No, 
                               aes(x = RNA_foldchange, y = RPF_foldchange)) +
  coord_cartesian(xlim = c(-3, 3),
                  ylim = c(-3, 3)) + 
  geom_point(size = 0.2, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = 'red',linetype="dashed") + 
  geom_abline(data = Fold_slopes_No, aes(intercept = Intercept, slope = slope), 
              color = 'blue') + 
  xlab('RNA abundance foldchange [log2]') + ylab('RPF abundance foldchange [log2]') + 
  ggtitle('uORF: No') + 
  mytheme_violin + 
  facet_wrap(~sample, scales = 'free',ncol = 5) + 
  stat_cor(label.x = -3, label.y = 2.5, method = "pearson",size = 2.5) + 
  geom_text(data = Fold_slopes_No, aes(label = paste0("Major Axis Slope: ", round(slope, 3))), 
            x = Inf, y = -2.5, hjust = 1, vjust = 0, size = 2.5) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
RPF_RNA_Fold_corr_No

RPF_RNA_Fold_corr_No_noWord <- ggplot(data = RPF_RNA_fold_No, 
                               aes(x = RNA_foldchange, y = RPF_foldchange)) +
  coord_cartesian(xlim = c(-3, 3),
                  ylim = c(-3, 3)) + 
  geom_point(size = 0.2, alpha = 0.5) +
  xlab('RNA abundance foldchange [log2]') + ylab('RPF abundance foldchange [log2]') + 
  ggtitle('uORF: No') + 
  mytheme_violin + 
  facet_wrap(~sample, scales = 'free',ncol = 5) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
RPF_RNA_Fold_corr_No_noWord

ggsave(
  filename = "/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/RPFfold_RNAfold_corr_MajorAxisSlope_uORFno-pearson.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 8, height = 3.5, units = "in", dpi = 4000
)

library("grid")
pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/RPFfold_RNAfold_corr_MajorAxisSlope_uORF-pearson.pdf",width=15,height=23,useDingbats = F)
print(RPF_RNA_Fold_corr_Yes,vp=viewport(.7,.2,x=.5,y=.85))
print(RPF_RNA_Fold_corr_No,vp=viewport(.7,.2,x=.5,y=.65))
dev.off()
#dev.new()


##### slope and 95%CI

Fold_slopes <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/Ribosome_Density/MajorAxis.txt',header = TRUE)
Fold_slopes_Yes <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/Ribosome_Density/MajorAxis_uORF_Yes.txt',header = TRUE)
Fold_slopes_No <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/Ribosome_Density/MajorAxis_uORF_No.txt',header = TRUE)

Fold_slopes <- rbind(Fold_slopes,Fold_slopes_Yes,Fold_slopes_No)
Fold_slopes$group <- factor(Fold_slopes$group, levels = c("uORF_No","uORF_Yes","All"))
Fold_slopes$sample <- factor(Fold_slopes$sample, levels = c("CAF20","EAP1","eIF4E","eIF4G1","PAB1","eIF1","eIF1A","eIF4A","eIF4B"))
head(Fold_slopes)

Fold_slopes_bar <- ggplot(Fold_slopes, aes(x = sample, y = slope, fill = group)) +
  geom_point(aes(colour = group), size = 0.8, 
             stat = "identity", position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = X2.5..Slope, ymax = X97.5..Slope, colour = group), 
                 stat = "identity", position = position_dodge(width = 0.5), 
                 linewidth = 0.6) +
  scale_fill_grey(start = .6, end = .0) +
  scale_colour_grey(start = .6, end = .0) +
  labs(x = "", y = "Major axis slope\n(RPF fold change ~ RNA fold change)") +
  ggtitle("Major axis slope and 95% CI") +
  mytheme_violin + 
  theme(axis.text.x = element_text(angle=0, hjust=0.5))
Fold_slopes_bar

Fold_slopes_bar_2 <- ggplot(Fold_slopes, aes(x = sample, y = slope, fill = group)) +
  geom_point(aes(colour = group), size = 0.8, 
             stat = "identity", position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = X2.5..Slope, ymax = X97.5..Slope, colour = group), 
                 stat = "identity", position = position_dodge(width = 0.5), 
                 linewidth = 0.6) +
  scale_fill_grey(start = .6, end = .0) +
  scale_colour_grey(start = .6, end = .0) +
  labs(x = "", y = "Major axis slope\n(RPF fold change ~ RNA fold change)") +
  ggtitle("Major axis slope and 95% CI") +
  coord_flip() +
  mytheme_violin + 
  theme(axis.text.x = element_text(angle=0, hjust=0.5))
Fold_slopes_bar_2

library("grid")
pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/Fold_slopes_point.pdf",width=10,height=6,useDingbats = F)
print(Fold_slopes_bar,vp=viewport(.6,.4,x=.33,y=.5))
print(Fold_slopes_bar_2,vp=viewport(.35,.65,x=.8,y=.5))
dev.off()
#dev.new()