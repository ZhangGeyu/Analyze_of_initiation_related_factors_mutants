mytheme_violin <-  theme_bw()+
  #theme(legend.position="none")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  theme(plot.title=element_text(size=8))+
  theme(axis.text.x = element_text(angle=0, hjust=0.5))+
  theme(axis.text.x = element_text(colour="black", size=8))+
  theme(axis.text.y = element_text(colour="black", size=8))+
  theme(axis.title=element_text(size=8))+
  #theme(plot.title=element_text(size=rel(1.5), lineheight=.9))+
  theme(axis.ticks=element_line(colour="black",size=0.5))+
  theme(axis.line=element_line(colour="black")) + 
  theme(legend.title=element_text(size=8))+
  theme(legend.text=element_text(size=8))


path = '/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/8-CooperativeGene/'


##### Five_vs_Four
## CDS Length

Cluster_length <- read.table(paste0(path,'Five_vs_Four','/Cluster_Length_and_GC.txt'),header = TRUE)

Cluster_CDS_length <- Cluster_length[Cluster_length$Region == 'CDS',]
Cluster_CDS_length <- Cluster_CDS_length[Cluster_CDS_length$group != 'Not_significant',]
head(Cluster_CDS_length)

Feature_CDS_length <- ggplot(Cluster_CDS_length, 
       aes(x = log2(length), colour = group)) + 
  scale_color_manual(values = c('red','blue')) + 
  geom_density(lwd = 0.3, adjust = 1.5) + 
  xlab('CDS length [log2]') + ggtitle('looping vs. scanning') + 
  scale_x_continuous(breaks = seq(6, 15, by = 1)) + 
  mytheme_violin + 
  theme(axis.text.x = element_text(angle=0, hjust=0.5))
Feature_CDS_length

Cluster_CDS_length$log2_length <- log2(Cluster_CDS_length$length)
UTR5_length_cumulative <- Cluster_CDS_length %>%
  group_by(group) %>%
  mutate(cumulative = ecdf(log2_length)(log2_length))
head(UTR5_length_cumulative)

Feature_CDS_length_CumulativeCurve <- ggplot(UTR5_length_cumulative, 
                                             aes(log2_length, cumulative, color = group)) +
  geom_step() + 
  scale_color_manual(values = c('red','blue')) + 
  xlab('CDS length [log2]') + 
  ylab('Cumulative Probability') + 
  scale_x_continuous(breaks = seq(6, 15, by = 1)) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
Feature_CDS_length_CumulativeCurve

Feature_CDS_length_histogram <- ggplot(data = Cluster_CDS_length, 
                                       aes(log2(length), fill = group)) +
  geom_histogram(aes(fill = group, y = after_stat(density)), position = 'identity', 
                 alpha = 0.85, bins = 30) +
  #geom_density(alpha = 0.5) +
  scale_fill_manual(values = c('red','blue')) + 
  xlab('CDS length [log2]') + 
  scale_x_continuous(breaks = seq(6, 15, by = 1)) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
Feature_CDS_length_histogram


## 5'UTR Length

Cluster_5UTR_length <- Cluster_length[Cluster_length$Region == "5'UTR",]
Cluster_5UTR_length <- Cluster_5UTR_length[Cluster_5UTR_length$group != 'Not_significant',]
head(Cluster_5UTR_length)

Feature_5UTR_length <- ggplot(Cluster_5UTR_length, 
                             aes(x = log2(length), colour = group)) + 
  scale_color_manual(values = c('red','blue')) + 
  geom_density(lwd = 0.3, adjust = 1.5) + 
  xlab('5UTR length [log2]') + ggtitle('looping vs. scanning') + 
  scale_x_continuous(breaks = seq(1, 11, by = 1)) + 
  mytheme_violin + 
  theme(axis.text.x = element_text(angle=0, hjust=0.5))
Feature_5UTR_length

Cluster_5UTR_length$log2_length <- log2(Cluster_5UTR_length$length)
UTR5_length_cumulative <- Cluster_5UTR_length %>%
  group_by(group) %>%
  mutate(cumulative = ecdf(log2_length)(log2_length))
head(UTR5_length_cumulative)

Feature_5UTR_length_CumulativeCurve <- ggplot(UTR5_length_cumulative, 
                              aes(log2_length, cumulative, color = group)) +
  geom_step() + 
  scale_color_manual(values = c('red','blue')) + 
  xlab('5UTR length [log2]') +
  ylab('Cumulative Probability') + 
  scale_x_continuous(breaks = seq(1, 11, by = 1)) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
Feature_5UTR_length_CumulativeCurve

Feature_5UTR_length_histogram <- ggplot(data = Cluster_5UTR_length, 
                                        aes(log2(length), fill = group)) +
  geom_histogram(aes(fill = group, y = after_stat(density)), position = 'identity', 
                 alpha = 0.85, bins = 30) +
  #geom_density(alpha = 0.5) +
  scale_fill_manual(values = c('red','blue')) + 
  xlab('5UTR length [log2]') +
  scale_x_continuous(breaks = seq(1, 11, by = 1)) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
Feature_5UTR_length_histogram


## Transcript Length

Cluster_transcript_length <- read.table(paste0(path,'Five_vs_Four','/Cluster_transcript_length.txt'),header = TRUE)
Cluster_transcript_length <- Cluster_transcript_length[Cluster_transcript_length$group != 'Not_significant',]
head(Cluster_transcript_length)

Feature_transcript_length <- ggplot(Cluster_transcript_length, 
                             aes(x = log2(transcript_length), colour = group)) + 
  scale_color_manual(values = c('red','blue')) + 
  geom_density(lwd = 0.3, adjust = 1.5) + 
  xlab('Transcript length [log2]') + ggtitle('looping vs. scanning') + 
  scale_x_continuous(breaks = seq(9, 16, by = 1)) + 
  mytheme_violin + 
  theme(axis.text.x = element_text(angle=0, hjust=0.5))
Feature_transcript_length

Cluster_transcript_length$log2_transcript_length <- log2(Cluster_transcript_length$transcript_length)
transcript_length_cumulative <- Cluster_transcript_length %>%
  group_by(group) %>%
  mutate(cumulative = ecdf(log2_transcript_length)(log2_transcript_length))
head(transcript_length_cumulative)

Feature_transcript_length_CumulativeCurve <- ggplot(transcript_length_cumulative, 
                                              aes(log2_transcript_length, cumulative, color = group)) +
  geom_step() + 
  scale_color_manual(values = c('red','blue')) + 
  xlab('Transcript length [log2]') +
  ylab('Cumulative Probability') + 
  scale_x_continuous(breaks = seq(9, 16, by = 1)) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
Feature_transcript_length_CumulativeCurve

Feature_transcript_length_histogram <- ggplot(data = Cluster_transcript_length, 
                                        aes(log2(transcript_length), fill = group)) +
  geom_histogram(aes(fill = group, y = after_stat(density)), position = 'identity', 
                 alpha = 0.85, bins = 30) +
  #geom_density(alpha = 0.5) +
  scale_fill_manual(values = c('red','blue')) + 
  xlab('Transcript length [log2]') +
  scale_x_continuous(breaks = seq(9, 16, by = 1)) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
Feature_transcript_length_histogram


## polyA length

Cluster_PolyA <- read.table(paste0(path,'Five_vs_Four','/Cluster_PolyA.txt'),header = TRUE)
Cluster_PolyA <- Cluster_PolyA[Cluster_PolyA$group != 'Not_significant',]
head(Cluster_PolyA)
length(Cluster_PolyA[Cluster_PolyA$group == 'Cluster2',]$gene)

Feature_PolyA_length <- ggplot(Cluster_PolyA, 
                               aes(x = poly_length_average, colour = group)) + 
  scale_color_manual(values = c('red','blue')) + 
  geom_density(lwd = 0.3, adjust = 1.5) + 
  xlab('Poly(A) length [log2]') + ggtitle('looping vs. scanning') + 
  scale_x_continuous(breaks = seq(25, 70, by = 5)) + 
  mytheme_violin + 
  theme(axis.text.x = element_text(angle=0, hjust=0.5))
Feature_PolyA_length

Cluster_PolyA$log2_poly_length <- log2(Cluster_PolyA$poly_length_average)
PolyA_length_cumulative <- Cluster_PolyA %>%
  group_by(group) %>%
  mutate(cumulative = ecdf(log2_poly_length)(log2_poly_length))
head(PolyA_length_cumulative)

PolyA_length_CumulativeCurve <- ggplot(PolyA_length_cumulative, 
                                       aes(log2_poly_length, cumulative, color = group)) +
  geom_step() + 
  scale_color_manual(values = c('red','blue')) + 
  xlab('Poly(A) length [log2]') +
  ylab('Cumulative Probability') + 
  scale_x_continuous(breaks = seq(5, 6, by = 0.2)) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
PolyA_length_CumulativeCurve

PolyA_length_histogram <- ggplot(data = Cluster_PolyA, 
                                 aes(poly_length_average, fill = group)) +
  geom_histogram(aes(fill = group, y = after_stat(density)), position = 'identity', 
                 alpha = 0.85, bins = 30) +
  #geom_density(alpha = 0.5) +
  scale_fill_manual(values = c('red','blue')) + 
  xlab('Poly(A) length') +
  scale_x_continuous(breaks = seq(25, 70, by = 5)) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
PolyA_length_histogram

## uORF enrichment

uORF_data <- data.frame(
  looping = c(190,15),
  scanning = c(291,73),
  group = c("uORF_No", "uORF_Yes")
)

uORF_data$looping <- uORF_data$looping / sum(uORF_data$looping)
uORF_data$scanning <- uORF_data$scanning / sum(uORF_data$scanning)

uORF_data_looping <- uORF_data[,c('looping','group')]
uORF_data_looping$cluster <- 'looping'
colnames(uORF_data_looping) <- c('fraction','group','cluster')

uORF_data_scanning <- uORF_data[,c('scanning','group')]
uORF_data_scanning$cluster <- 'scanning'
colnames(uORF_data_scanning) <- c('fraction','group','cluster')

uORF_data <- rbind(uORF_data_looping,uORF_data_scanning)
uORF_data$group <- factor(uORF_data$group, levels = c("uORF_Yes","uORF_No"))
head(uORF_data)

uORF_enrichment <- ggplot(uORF_data, aes(x = cluster, y = fraction, fill = group)) +
  geom_bar(stat='identity', alpha = 0.7, width = 0.6) +
  labs(x = "Cluster", y = "Fraction of genes have uORF") +
  scale_fill_grey(start = .3, end = .6) +
  scale_colour_grey(start = .3, end = .6) +
  mytheme_violin + 
  theme(axis.text.x = element_text(angle=0, hjust=0.5))
uORF_enrichment

library("grid")
pdf(paste0('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/Cluster_','Five_vs_Four','.pdf'),width=18,height=6,useDingbats = F, colormodel = "rgb")
print(Feature_CDS_length,vp=viewport(.2,.25,x=.1,y=.8))
print(Feature_transcript_length,vp=viewport(.2,.25,x=.3,y=.8))
print(Feature_5UTR_length,vp=viewport(.2,.25,x=.5,y=.8))
print(Feature_PolyA_length,vp=viewport(.2,.25,x=.7,y=.8))
print(uORF_enrichment,vp=viewport(.15,.25,x=.9,y=.8))

print(Feature_CDS_length_CumulativeCurve,vp=viewport(.2,.25,x=.1,y=.5))
print(Feature_transcript_length_CumulativeCurve,vp=viewport(.2,.25,x=.3,y=.5))
print(Feature_5UTR_length_CumulativeCurve,vp=viewport(.2,.25,x=.5,y=.5))
print(PolyA_length_CumulativeCurve,vp=viewport(.2,.25,x=.7,y=.5))

print(Feature_CDS_length_histogram,vp=viewport(.2,.25,x=.1,y=.2))
print(Feature_transcript_length_histogram,vp=viewport(.2,.25,x=.3,y=.2))
print(Feature_5UTR_length_histogram,vp=viewport(.2,.25,x=.5,y=.2))
print(PolyA_length_histogram,vp=viewport(.2,.25,x=.7,y=.2))
dev.off()


##### Two_vs_Two
## 5'UTR Length and GC%

Cluster_length_GC <- read.table(paste0(path,'Two_vs_Two','/Cluster_Length_and_GC.txt'),header = TRUE)
Cluster_length_GC <- Cluster_length_GC[Cluster_length_GC$Region == "5'UTR",]
Cluster_length_GC <- Cluster_length_GC[Cluster_length_GC$group != 'Not_significant',]
head(Cluster_length_GC)

length(Cluster_length_GC[Cluster_length_GC$group == 'Cluster1',]$gene)


Feature_5UTR_length <- ggplot(Cluster_length_GC, 
                             aes(x = log2(length), colour = group)) + 
  scale_color_manual(values = c('orange','purple')) + 
  geom_density(lwd = 0.3, adjust = 1.5) + 
  xlab("5'UTR length [log2]") + ggtitle('1 1A vs. 4A 4B') + 
  scale_x_continuous(breaks = seq(2, 10, by = 1)) + 
  mytheme_violin + 
  theme(axis.text.x = element_text(angle=0, hjust=0.5))
Feature_5UTR_length

Cluster_length_GC$log2_length <- log2(Cluster_length_GC$length)
Feature_5UTR_length_cumulative <- Cluster_length_GC %>%
  group_by(group) %>%
  mutate(cumulative = ecdf(log2_length)(log2_length))
head(Feature_5UTR_length_cumulative)

Feature_5UTR_length_CumulativeCurve <- ggplot(Feature_5UTR_length_cumulative, 
                                       aes(log2_length, cumulative, color = group)) +
  geom_step() + 
  scale_color_manual(values = c('orange','purple')) + 
  xlab("5'UTR length [log2]") +
  ylab('Cumulative Probability') + 
  scale_x_continuous(breaks = seq(2, 10, by = 1)) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
Feature_5UTR_length_CumulativeCurve

Feature_5UTR_length_histogram <- ggplot(data = Cluster_length_GC, 
                                 aes(log2(length), fill = group)) +
  geom_histogram(aes(fill = group, y = after_stat(density)), position = 'identity', 
                 alpha = 0.85, bins = 30) +
  #geom_density(alpha = 0.5) +
  scale_fill_manual(values = c('orange','purple')) + 
  xlab("5'UTR length [log2]") +
  scale_x_continuous(breaks = seq(2, 10, by = 1)) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
Feature_5UTR_length_histogram


Feature_5UTR_GC <- ggplot(Cluster_length_GC, 
                              aes(x = GC, colour = group)) + 
  scale_color_manual(values = c('orange','purple')) + 
  geom_density(lwd = 0.3, adjust = 1.5) + 
  xlab("5'UTR GC%") + ggtitle('1 1A vs. 4A 4B') + 
  mytheme_violin + 
  theme(axis.text.x = element_text(angle=0, hjust=0.5))
Feature_5UTR_GC

Feature_GC_cumulative <- Cluster_length_GC %>%
  group_by(group) %>%
  mutate(cumulative = ecdf(GC)(GC))
head(Feature_GC_cumulative)

Feature_GC_CumulativeCurve <- ggplot(Feature_GC_cumulative, 
                                     aes(GC, cumulative, color = group)) +
  geom_step() + 
  scale_color_manual(values = c('orange','purple')) + 
  xlab("5'UTR GC%") +
  ylab('Cumulative Probability') + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
Feature_GC_CumulativeCurve

Feature_GC_histogram <- ggplot(data = Cluster_length_GC, 
                               aes(GC, fill = group)) +
  geom_histogram(aes(fill = group, y = after_stat(density)), position = 'identity', 
                 alpha = 0.85, bins = 30) +
  #geom_density(alpha = 0.5) +
  scale_fill_manual(values = c('orange','purple')) + 
  xlab("5'UTR GC%") +
  #scale_x_continuous(breaks = seq(2, 10, by = 1)) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
Feature_GC_histogram


library("grid")
pdf(paste0('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/Cluster_','Two_vs_Two','.pdf'),width=12,height=6,useDingbats = F, colormodel = "rgb")
print(Feature_5UTR_length,vp=viewport(.3,.25,x=.2,y=.8))
print(Feature_5UTR_GC,vp=viewport(.3,.25,x=.5,y=.8))

print(Feature_5UTR_length_CumulativeCurve,vp=viewport(.3,.25,x=.2,y=.5))
print(Feature_GC_CumulativeCurve,vp=viewport(.3,.25,x=.5,y=.5))

print(Feature_5UTR_length_histogram,vp=viewport(.3,.25,x=.2,y=.2))
print(Feature_GC_histogram,vp=viewport(.3,.25,x=.5,y=.2))
dev.off()


##### Feature conclusion

### Utilize the wilcox.test to calculate whether there are significant differences in all features

path = '/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/8-CooperativeGene/'
sample = 'Five_vs_Four'

Feature_conclusion_data <- read.table(paste0(path,sample,'/Feature_rbind_3kmer.txt'))
head(Feature_conclusion_data)
unique(Feature_conclusion_data$Feature)

Kmer_mean_pval <- file(paste0(path,sample,'/Feature_mean_pval_3mer_wilcox.txt'), "w")
writeLines(paste0('Feature',"\t",'normalized_scale',"\t",'pval'), con = Kmer_mean_pval)

for (i in unique(Feature_conclusion_data$Feature)) {
  Feature_rbind_subset <- Feature_conclusion_data %>%
    filter(Feature == i)
  
  Feature_Cluster1 <- Feature_rbind_subset %>%
    filter(group == 'Cluster1') %>%
    pull(value)
  mean_Cluster1 <- mean(Feature_Cluster1)
  
  Feature_Cluster2 <- Feature_rbind_subset %>%
    filter(group == 'Cluster2') %>%
    pull(value)
  mean_Cluster2 <- mean(Feature_Cluster2)
  
  #mean_Feature <- mean(c(Feature_Cluster1, Feature_Cluster2))
  #difference_Feature <- mean_Cluster2 - mean_Cluster1
  #normalized_scale <- difference_Feature / mean_Feature
  normalized_scale <- log2(mean_Cluster2/mean_Cluster1)
  
  wilcox_test_result <- wilcox.test(Feature_Cluster1,Feature_Cluster2,exact = FALSE)
  pval <- wilcox_test_result$p.value
  
  cat(i, normalized_scale, pval,'\n', file = Kmer_mean_pval, sep = "\t")
}

close(Kmer_mean_pval)

Kmer_mean_pval <- read.table(paste0(path,sample,'/Feature_mean_pval_3mer_wilcox.txt'),header = TRUE)

for (n in 1:length(Kmer_mean_pval$Feature)){
  if (Kmer_mean_pval$pval[n] <= 0.05){
    Kmer_mean_pval$group[n] = 'significant'
  }else {
    Kmer_mean_pval$group[n] = 'Not_significant'
  }
}

Kmer_mean_pval[1:13,]$group = 'not_kmer'
head(Kmer_mean_pval)

write.table(Kmer_mean_pval,paste0(path,sample,'/Feature_mean_pval_3mer_wilcox.txt'),sep = '\t', quote = FALSE)


### volcano plot

Feature_3mer <- read.table(paste0(path,sample,'/Feature_mean_pval_3mer_wilcox.txt'))
head(Feature_3mer)

Feature_3mer_point <- ggplot(data = Feature_3mer, 
                        aes(x = normalized_scale, y = -log10(pval),
                            colour = group)) + 
  #coord_cartesian(xlim = c(-2, 2)) + 
  scale_color_manual(values = c('blue','grey','black')) + 
  geom_point(size = 1) + 
  geom_hline(aes(yintercept = -log10(0.05)), linetype = 'dashed') + 
  geom_label_repel(data = subset(Feature_3mer, group != "Not_significant"),
                   aes(label = Feature, colour = group), 
                   size = 2) + 
  xlab('log2(looping-dependent genes average / scanning-dependent genes average)') + ylab('p value [-log10]') + 
  ggtitle('3-mer') + 
  mytheme_violin + 
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"),
        legend.position="none")
Feature_3mer_point


library("grid")
pdf(paste0('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/Five_vs_Four_Feature_conclusion_wilcox.test.pdf'),width=15,height=7,useDingbats = F, colormodel = "rgb")
print(Feature_3mer_point,vp=viewport(.5,1,x=.25,y=.5))
dev.off()

##### correlation: 3-mer% and mRNA abundance

path = '/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/8-CooperativeGene/'
sample = 'Five_vs_Four'

kmer_data <- read.table(paste0(path,sample,'/Cluster_3mer.txt'))
head(kmer_data)

CCC_freq <- kmer_data[kmer_data$Kmer == 'CCC',]
CGG_freq <- kmer_data[kmer_data$Kmer == 'CGG',]
GCG_freq <- kmer_data[kmer_data$Kmer == 'GCG',]
CGC_freq <- kmer_data[kmer_data$Kmer == 'CGC',]
head(CCC_freq)

mRNA_abundance <- read.table(paste0(path,sample,'/Cluster_HO_mRNA_abundance.txt'))
mRNA_abundance <- mRNA_abundance[,c(1,2)]
head(mRNA_abundance)

CCC_freq <- merge(CCC_freq, mRNA_abundance, by = 'gene')
CGG_freq <- merge(CGG_freq, mRNA_abundance, by = 'gene')
GCG_freq <- merge(GCG_freq, mRNA_abundance, by = 'gene')
CGC_freq <- merge(CGC_freq, mRNA_abundance, by = 'gene')
head(CCC_freq)

CCC_freq_mRNA_corr <- ggplot(data = CCC_freq, 
                             aes(x = value, y = log2(RNA_abundance))) +
  #coord_cartesian(ylim = c(-2, 2.1)) + 
  geom_point(size = 0.2, alpha = 0.5) +
  geom_text(label = paste("n = ", length(CCC_freq$gene)),
            x = 0.025, y = -30,
            size = 3, fontface = "plain") + 
  xlab('CCC%') + ylab('mRNA abundance [log2]') + 
  mytheme_violin + 
  stat_cor(label.x = 0, label.y = -15, method = "spearman",size = 2.8) + 
  #scale_x_continuous(breaks = seq(25, 70, by = 5)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
CCC_freq_mRNA_corr

CGG_freq_mRNA_corr <- ggplot(data = CGG_freq, 
                             aes(x = value, y = log2(RNA_abundance))) +
  #coord_cartesian(ylim = c(-2, 2.1)) + 
  geom_point(size = 0.2, alpha = 0.5) +
  geom_text(label = paste("n = ", length(CGG_freq$gene)),
            x = 0.015, y = -30,
            size = 3, fontface = "plain") + 
  xlab('CGG%') + ylab('mRNA abundance [log2]') + 
  mytheme_violin + 
  stat_cor(label.x = 0, label.y = -15, method = "spearman",size = 2.8) + 
  #scale_x_continuous(breaks = seq(25, 70, by = 5)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
CGG_freq_mRNA_corr

GCG_freq_mRNA_corr <- ggplot(data = GCG_freq, 
                             aes(x = value, y = log2(RNA_abundance))) +
  #coord_cartesian(ylim = c(-2, 2.1)) + 
  geom_point(size = 0.2, alpha = 0.5) +
  geom_text(label = paste("n = ", length(GCG_freq$gene)),
            x = 0.015, y = -30,
            size = 3, fontface = "plain") + 
  xlab('GCG%') + ylab('mRNA abundance [log2]') + 
  mytheme_violin + 
  stat_cor(label.x = 0, label.y = -15, method = "spearman",size = 2.8) + 
  #scale_x_continuous(breaks = seq(25, 70, by = 5)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
GCG_freq_mRNA_corr

CGC_freq_mRNA_corr <- ggplot(data = CGC_freq, 
                             aes(x = value, y = log2(RNA_abundance))) +
  #coord_cartesian(ylim = c(-2, 2.1)) + 
  geom_point(size = 0.2, alpha = 0.5) +
  geom_text(label = paste("n = ", length(CGC_freq$gene)),
            x = 0.015, y = -30,
            size = 3, fontface = "plain") + 
  xlab('CGC%') + ylab('mRNA abundance [log2]') + 
  mytheme_violin + 
  stat_cor(label.x = 0, label.y = -15, method = "spearman",size = 2.8) + 
  #scale_x_continuous(breaks = seq(25, 70, by = 5)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
CGC_freq_mRNA_corr

CDS_GC <- read.table(paste0(path,sample,'/Cluster_Length_and_GC.txt'))
CDS_GC <- CDS_GC[CDS_GC$Region == 'CDS',]
CDS_GC <- CDS_GC[,c(1,4)]
CDS_GC <- merge(CDS_GC, mRNA_abundance, by = 'gene')
head(CDS_GC)

CDS_GC_mRNA_corr <- ggplot(data = CDS_GC, 
                             aes(x = GC, y = log2(RNA_abundance))) +
  #coord_cartesian(ylim = c(-2, 2.1)) + 
  geom_point(size = 0.2, alpha = 0.5) +
  geom_text(label = paste("n = ", length(CDS_GC$gene)),
            x = 50, y = -30,
            size = 3, fontface = "plain") + 
  xlab('CDS GC%') + ylab('mRNA abundance [log2]') + 
  mytheme_violin + 
  stat_cor(label.x = 24, label.y = -15, method = "spearman",size = 2.8) + 
  #scale_x_continuous(breaks = seq(25, 70, by = 5)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
CDS_GC_mRNA_corr


library("grid")
pdf(paste0('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/Five_vs_Four_Kmer_mRNA_corr.pdf'),width=12,height=2,useDingbats = F, colormodel = "rgb")
print(CCC_freq_mRNA_corr,vp=viewport(.2,1,x=.1,y=.5))
print(CGG_freq_mRNA_corr,vp=viewport(.2,1,x=.3,y=.5))
print(GCG_freq_mRNA_corr,vp=viewport(.2,1,x=.5,y=.5))
print(CGC_freq_mRNA_corr,vp=viewport(.2,1,x=.7,y=.5))
print(CDS_GC_mRNA_corr,vp=viewport(.2,1,x=.9,y=.5))
dev.off()

