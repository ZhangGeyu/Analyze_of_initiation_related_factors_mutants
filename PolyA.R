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

##### polyA Length and Ribosome Density fold change

PolyA_Length_RibosomeDensityfold <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/7-PolyA/PolyALength_and_RibosomeDensity-fold.txt',header = TRUE)
PolyA_Length_RibosomeDensityfold <- PolyA_Length_RibosomeDensityfold[PolyA_Length_RibosomeDensityfold$sample != 'HO',]
PolyA_Length_RibosomeDensityfold$sample <- factor(PolyA_Length_RibosomeDensityfold$sample, levels = c("CAF20","EAP1","eIF4E","eIF4G1","PAB1","eIF1","eIF1A","eIF4A","eIF4B"))
PolyA_Length_RibosomeDensityfold$group <- factor(PolyA_Length_RibosomeDensityfold$bin_label, levels = c('bin1','bin2','bin3','bin4','bin5','bin6','bin7','bin8'))
colnames(PolyA_Length_RibosomeDensityfold)[2] <- "PolyA_length"
head(PolyA_Length_RibosomeDensityfold)

PolyALength_RibosomeDensityfold <- ggplot(data = PolyA_Length_RibosomeDensityfold,
                                             aes(x = PolyA_length, 
                                                 y = log2(RibosomeDensity_foldchange))) + 
  coord_cartesian(ylim = c(-1.6,1.6)) + 
  geom_boxplot(outlier.size = 0.5, fill = 'white',notch = TRUE,
               outlier.shape = NA, colour = 'black') + 
  ylab("TE foldchange [log2]") + 
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = 'red') +
  facet_wrap(~sample, ncol = 1, scales = 'free') + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))+ 
  theme(axis.text.x=element_blank())
PolyALength_RibosomeDensityfold

point_counts <- PolyA_Length_RibosomeDensityfold %>%
  group_by(sample) %>%
  summarize(count = n())
point_counts

PolyALength_RibosomeDensityfold_corr <- ggplot(data = PolyA_Length_RibosomeDensityfold,
                                          aes(x = log2(poly_length_average), 
                                              y = log2(RibosomeDensity_foldchange))) + 
  #coord_cartesian(ylim = c(-2,2)) + 
  geom_point(size = 0.2, alpha = 0.5) +
  xlab('polyA length [log2]') + ylab('TE foldchange [log2]') + 
  mytheme_violin + 
  facet_wrap(~sample, scales = 'free',ncol = 1) + 
  stat_cor(label.x = 4.85, label.y = 1, method = "spearman",size = 2) + 
  geom_text(data = point_counts, aes(label = paste("n =", count)),
            x = 5, y = -2,
            size = 2, color = "black") + 
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
PolyALength_RibosomeDensityfold_corr


PolyA_Length_RibosomeDensityfold$log2_RibosomeDensity_foldchange <- log2(PolyA_Length_RibosomeDensityfold$RibosomeDensity_foldchange)
head(PolyA_Length_RibosomeDensityfold)

PolyA_Length_MedianAndCI <- file('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/PolyA_Length_MedianAndCI.txt', "w")
writeLines(paste0('sample','\t','bin','\t','median','\t','CI_low','\t','CI_high'), con = PolyA_Length_MedianAndCI)

for (i in unique(PolyA_Length_RibosomeDensityfold$sample)){
  for (j in unique(PolyA_Length_RibosomeDensityfold$PolyA_length)){
    data_subset = PolyA_Length_RibosomeDensityfold[(PolyA_Length_RibosomeDensityfold$sample == i)&
                                                     (PolyA_Length_RibosomeDensityfold$PolyA_length == j),]$log2_RibosomeDensity_foldchange
    sample = i
    bin = j
    median = median(data_subset)
    CI_low = boxplot.stats(data_subset)$conf[1]
    CI_high = boxplot.stats(data_subset)$conf[2]
    cat(sample, bin, median,CI_low,CI_high,'\n', file = PolyA_Length_MedianAndCI, sep = "\t")
  }
}
close(PolyA_Length_MedianAndCI)


PolyA_Length_MedianAndCI <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/PolyA_Length_MedianAndCI.txt',header = TRUE)
PolyA_Length_MedianAndCI$sample <- factor(PolyA_Length_MedianAndCI$sample, levels = c("CAF20","EAP1","eIF4E","eIF4G1","PAB1","eIF1","eIF1A","eIF4A","eIF4B"))
head(PolyA_Length_MedianAndCI)

PolyALength_RibosomeDensityfold_point <- ggplot(PolyA_Length_MedianAndCI, 
                                              aes(x = bin, y = median)) +
  coord_cartesian(ylim = c(-0.35,0.5)) + 
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
PolyALength_RibosomeDensityfold_point

pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/PolyA_length_8bins-TEfoldchange.pdf",width=7.5,height=13,useDingbats = F, colormodel = "cmyk")
print(PolyALength_RibosomeDensityfold,vp=viewport(.3,1,x=.2,y=.5))
print(PolyALength_RibosomeDensityfold_corr,vp=viewport(.3,1,x=.5,y=.5))
print(PolyALength_RibosomeDensityfold_point,vp=viewport(.3,1,x=.8,y=.5))
dev.off()


##### polyA Length distribution

polyA_length_nonA = read.table('/Dell/Dell15/zhanggy/project/3-mRNA_looping/tmp/Pacbio_data/PolyA_analyze/H_polyLength_nonAcontent.txt')
head(polyA_length_nonA)

PolyALength_distribution <- ggplot(polyA_length_nonA, 
                             aes(x = poly_length_average)) + 
  scale_color_manual(values = c('red','blue')) + 
  geom_density(lwd = 0.3) + 
  xlab('PolyA Length') + ggtitle('Wild Type') + 
  mytheme_violin + 
  scale_x_continuous(breaks = seq(25, 70, by = 5)) + 
  theme(axis.text.x = element_text(angle=0, hjust=0.5))
PolyALength_distribution

PolyA_length_histogram <- ggplot(data = polyA_length_nonA, 
                                 aes(poly_length_average)) +
  geom_histogram(aes(y = after_stat(density)), position = 'identity', 
                 alpha = 0.7, bins = 30, fill = 'black') +
  xlab('Poly(A) length') +
  scale_x_continuous(breaks = seq(25, 70, by = 5)) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
PolyA_length_histogram

pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/PolyA-Length_distribution.pdf",width=5,height=1.5,useDingbats = F, colormodel = "cmyk")
print(PolyALength_distribution,vp=viewport(.5,1,x=.25,y=.5))
print(PolyA_length_histogram,vp=viewport(.5,1,x=.75,y=.5))
dev.off()


##### polyA Length and nonA content

polyA_length_nonA = read.table('/Dell/Dell15/zhanggy/project/3-mRNA_looping/tmp/Pacbio_data/PolyA_analyze/H_polyLength_nonAcontent.txt')
head(polyA_length_nonA)

#polyA_length_nonA$purity <- (1 - polyA_length_nonA$nonA_freq_average)

PolyALength_nonA_purity <- ggplot(data = polyA_length_nonA,
                                  aes(x = poly_length_average, y = log2(nonA_freq_average))) + 
  #coord_cartesian(ylim = c(0.992,1)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  xlab("polyA length [log2]") + ylab("nonA frequency [log2]") + 
  #stat_cor(label.x = 26, label.y = 0.001, method = "spearman",size = 2.8) + 
  scale_x_continuous(breaks = seq(25, 70, by = 5)) + 
  scale_y_continuous(breaks = seq(-14, -6, by = 1)) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
PolyALength_nonA_purity

pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/PolyA-Length_and_nonA.pdf",width=2.2,height=2,useDingbats = F, colormodel = "rgb")
print(PolyALength_nonA_purity,vp=viewport(1,1,x=.5,y=.5))
dev.off()
#dev.new()
