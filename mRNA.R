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


##### CDS Length and Ribosome Density foldchange

CDS_Length_RibosomeDensityfold <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/9-mRNA/CDSLength_and_RibosomeDensity-fold.txt',header = TRUE)
CDS_Length_RibosomeDensityfold <- CDS_Length_RibosomeDensityfold[CDS_Length_RibosomeDensityfold$sample != 'HO',]
CDS_Length_RibosomeDensityfold$sample <- factor(CDS_Length_RibosomeDensityfold$sample, levels = c("CAF20","EAP1","eIF4E","eIF4G1","PAB1","eIF1","eIF1A","eIF4A","eIF4B"))
colnames(CDS_Length_RibosomeDensityfold)[2] <- "CDS_length"
head(CDS_Length_RibosomeDensityfold)

CDSLength_RibosomeDensityfold <- ggplot(data = CDS_Length_RibosomeDensityfold,
                                          aes(x = CDS_length, 
                                              y = log2(RibosomeDensity_foldchange))) + 
  coord_cartesian(ylim = c(-2.5,2.5)) + 
  geom_boxplot(outlier.size = 0.5, fill = 'white',notch = TRUE,
               outlier.shape = NA, colour = 'black') + 
  ylab("TE foldchange [log2]") + xlab('bin') +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = 'red') +
  facet_wrap(~sample, ncol = 1, scales = 'free') + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))+ 
  theme(axis.text.x=element_blank())
CDSLength_RibosomeDensityfold

point_counts <- CDS_Length_RibosomeDensityfold %>%
  group_by(sample) %>%
  summarize(count = n())
point_counts

CDSLength_RibosomeDensityfold_corr <- ggplot(data = CDS_Length_RibosomeDensityfold,
                                               aes(x = log2(length), 
                                                   y = log2(RibosomeDensity_foldchange))) + 
  #coord_cartesian(ylim = c(-2,2)) + 
  geom_point(size = 0.2, alpha = 0.5) +
  xlab('CDS length [log2]') + ylab('TE foldchange [log2]') + 
  mytheme_violin + 
  facet_wrap(~sample, scales = 'free',ncol = 1) + 
  stat_cor(label.x = 4.85, label.y = 3, method = "spearman",size = 2) + 
  geom_text(data = point_counts, aes(label = paste("n =", count)),
            x = 6, y = 2,
            size = 2, color = "black") + 
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
CDSLength_RibosomeDensityfold_corr

CDS_Length_RibosomeDensityfold$log2_RibosomeDensity_foldchange <- log2(CDS_Length_RibosomeDensityfold$RibosomeDensity_foldchange)

CDS_Length_MedianAndCI <- file('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/CDS_Length_MedianAndCI.txt', "w")
writeLines(paste0('sample','\t','bin','\t','median','\t','CI_low','\t','CI_high'), con = CDS_Length_MedianAndCI)

for (i in unique(CDS_Length_RibosomeDensityfold$sample)){
  for (j in unique(CDS_Length_RibosomeDensityfold$CDS_length)){
    data_subset = CDS_Length_RibosomeDensityfold[(CDS_Length_RibosomeDensityfold$sample == i)&
                                                   (CDS_Length_RibosomeDensityfold$CDS_length == j),]$log2_RibosomeDensity_foldchange
    sample = i
    bin = j
    median = median(data_subset)
    CI_low = boxplot.stats(data_subset)$conf[1]
    CI_high = boxplot.stats(data_subset)$conf[2]
    cat(sample, bin, median,CI_low,CI_high,'\n', file = CDS_Length_MedianAndCI, sep = "\t")
  }
  }
close(CDS_Length_MedianAndCI)


CDS_Length_MedianAndCI <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/CDS_Length_MedianAndCI.txt',header = TRUE)
CDS_Length_MedianAndCI$sample <- factor(CDS_Length_MedianAndCI$sample, levels = c("CAF20","EAP1","eIF4E","eIF4G1","PAB1","eIF1","eIF1A","eIF4A","eIF4B"))
head(CDS_Length_MedianAndCI)

CDSLength_RibosomeDensityfold_point <- ggplot(CDS_Length_MedianAndCI, 
                                        aes(x = bin, y = median)) +
  coord_cartesian(ylim = c(-0.4,0.4)) + 
  ylab("TE foldchange [log2]") + xlab('bin') +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = 'red', lwd = 0.5) +
  geom_point(stat = "summary", fun = "median", shape = 16, size = 1, color = "black") + 
  geom_linerange(aes(ymin = CI_low, ymax = CI_high), 
                 stat = "identity", position = position_dodge(width = 0.5), 
                 linewidth = 0.5) +
  facet_wrap(~sample, ncol = 1, scales = 'free') + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))+ 
  theme(axis.text.x=element_blank())
CDSLength_RibosomeDensityfold_point


pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/CDS_length_8bins-TEfoldchange.pdf",width=7.5,height=13,useDingbats = F, colormodel = "cmyk")
print(CDSLength_RibosomeDensityfold,vp=viewport(.3,1,x=.2,y=.5))
print(CDSLength_RibosomeDensityfold_corr,vp=viewport(.3,1,x=.5,y=.5))
print(CDSLength_RibosomeDensityfold_point,vp=viewport(.3,1,x=.8,y=.5))
dev.off()

