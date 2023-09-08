
##### RNA abundance: example----between replication/samples

RNA_abundance <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/RNA_Abundance/TotalMutant-RNA_Abundance.txt')
head(RNA_abundance)

eIF4A_rep_corr <- ggplot(data = RNA_abundance, 
                      aes(x = log2(eIF4A.1), y = log2(eIF4A.2))) +
  coord_cartesian(xlim = c(-35, -12),
                  ylim = c(-35, -12)) + 
  geom_point(size = 0.4, alpha = 0.5, shape = 16) +
  xlab('eIF4A-1: RNA Abundance [log2]') + ylab('eIF4A-2: RNA Abundance [log2]') + 
  mytheme_violin + 
  #stat_cor(label.x = -34, label.y = -14, method = "spearman",size = 2.5) + 
  #geom_text(label = paste0("n = ", length(RNA_abundance$eIF4A.1) - 329),
  #          x = -32, y = -16,
  #          size = 2.5, fontface = "plain") + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
eIF4A_rep_corr

ggsave(
  filename = "/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/RNA_abundance_example_spearman-eIF4A_rep_corr.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 2, height = 2, units = "in", dpi = 4000
)

eIF4A1_eIF1A1_corr <- ggplot(data = RNA_abundance, 
                      aes(x = log2(eIF4A.1), y = log2(eIF1A.1))) +
  coord_cartesian(xlim = c(-35, -12),
                  ylim = c(-35, -12)) + 
  geom_point(size = 0.4, alpha = 0.5, shape = 16) +
  xlab('eIF4A-1: RNA Abundance [log2]') + ylab('eIF1A-1: RNA Abundance [log2]') + 
  mytheme_violin + 
  #stat_cor(label.x = -34, label.y = -14, method = "spearman",size = 2.5) + 
  #geom_text(label = paste0("n = ", length(RNA_abundance$eIF4A.1) - 313),
  #          x = -32, y = -16,
  #          size = 2.5, fontface = "plain") + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
eIF4A1_eIF1A1_corr

ggsave(
  filename = "/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/RNA_abundance_example_spearman-eIF4A1_eIF1A1_corr.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 2, height = 2, units = "in", dpi = 4000
)

pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/RNA_abundance_example_spearman.pdf",width=5.4,height=2.5,useDingbats = F, colormodel = "rgb")
print(eIF4A_rep_corr,vp=viewport(.5,1,x=.25,y=.5))
print(eIF4A1_eIF1A1_corr,vp=viewport(.5,1,x=.75,y=.5))
dev.off()


##### RPF abundance: example----between replication/samples


RPF_abundance <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/Ribo_Abundance/TotalMutant-Ribo_Abundance.txt')
head(RPF_abundance)

eIF4A_rep_corr <- ggplot(data = RPF_abundance, 
                      aes(x = log2(eIF4A.1), y = log2(eIF4A.2))) +
  coord_cartesian(xlim = c(-32, -12),
                  ylim = c(-32, -12)) + 
  geom_point(size = 0.4, alpha = 0.5, shape = 16) +
  xlab('eIF4A-1: RPF Abundance [log2]') + ylab('eIF4A-2: RPF Abundance [log2]') + 
  mytheme_violin + 
  stat_cor(label.x = -32, label.y = -14, method = "spearman", size = 2.5) + 
  geom_text(label = paste0("n = ", length(RPF_abundance$eIF4A.1) - 1074),
            x = -30, y = -16,
            size = 2.5, fontface = "plain") + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
eIF4A_rep_corr

eIF4A1_eIF4G11_corr <- ggplot(data = RPF_abundance, 
                          aes(x = log2(eIF4A.1), y = log2(eIF4G1.1))) +
  coord_cartesian(xlim = c(-32, -12),
                  ylim = c(-32, -12)) + 
  geom_point(size = 0.4, alpha = 0.5, shape = 16) +
  xlab('eIF4A-1: RPF Abundance [log2]') + ylab('eIF4G1-1: RPF Abundance [log2]') + 
  mytheme_violin + 
  stat_cor(label.x = -32, label.y = -14, method = "spearman",size = 2.5) + 
  geom_text(label = paste0("n = ", length(RPF_abundance$eIF4A.1) - 1001),
            x = -30, y = -16,
            size = 2.5, fontface = "plain") + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
eIF4A1_eIF4G11_corr

pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/RPF_abundance_example_spearman.pdf",width=5.4,height=2.5,useDingbats = F, colormodel = "rgb")
print(eIF4A_rep_corr,vp=viewport(.5,1,x=.25,y=.5))
print(eIF4A1_eIF4G11_corr,vp=viewport(.5,1,x=.75,y=.5))
dev.off()

##### TE fold change: eIF4E ~ CAF20/EAP1

TEfoldchange <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/Ribosome_Density/TotalMutant-RibosomeDensity-FoldChange.txt')
head(TEfoldchange)

eIF4E_CAF20_corr <- ggplot(data = TEfoldchange, 
                               aes(x = log2(eIF4E), y = log2(CAF20))) +
  coord_cartesian(xlim = c(-3, 3),
                  ylim = c(-3, 3)) + 
  geom_text(label = paste0("n = ", length(TEfoldchange$eIF4E) - 182),
            x = -2, y = 2,
            size = 3, fontface = "plain") + 
  geom_point(size = 0.2, alpha = 0.5) +
  xlab('eIF4E TE fold change [log2]') + ylab('CAF20 TE fold change [log2]') + 
  mytheme_violin + 
  stat_cor(label.x = -3, label.y = 3, method = "spearman",size = 2.3) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
eIF4E_CAF20_corr

eIF4E_EAP1_corr <- ggplot(data = TEfoldchange, 
                           aes(x = log2(eIF4E), y = log2(EAP1))) +
  coord_cartesian(xlim = c(-3, 3),
                  ylim = c(-3, 3)) + 
  geom_text(label = paste0("n = ", length(TEfoldchange$eIF4E) - 173),
            x = -2, y = 2,
            size = 3, fontface = "plain") + 
  geom_point(size = 0.2, alpha = 0.5) +
  xlab('eIF4E TE fold change [log2]') + ylab('EAP1 TE fold change [log2]') + 
  mytheme_violin + 
  stat_cor(label.x = -3, label.y = 3, method = "spearman", size = 2.3) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
eIF4E_EAP1_corr


pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/eIF4E_CAF20_EAP1_spearman.pdf",width=4.3,height=2,useDingbats = F, colormodel = "rgb")
print(eIF4E_CAF20_corr,vp=viewport(.5,1,x=.25,y=.5))
print(eIF4E_EAP1_corr,vp=viewport(.5,1,x=.75,y=.5))
dev.off()


##### Correlation histogram

RNA_cor_group <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/RNA_Abundance/TotalMutant-RNA-Density-corr.txt')
head(RNA_cor_group)

RNA_histogram <- ggplot(data = RNA_cor_group, aes(x = rho, fill = group)) +
  geom_histogram(aes(fill = group, y = after_stat(density)), position = 'identity', 
                 alpha = 0.85, bins = 20) +
  scale_fill_manual(values = c('gold2','cyan4')) + 
  xlab("Spearman's correlation coefficient between two rna-seq samples") + 
  ggtitle('RNA-seq') + 
  scale_x_continuous(breaks = seq(0.88, 1, by = 0.02)) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
RNA_histogram


RPF_cor_group <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/3-RibosomeDensity/Ribo_Abundance/TotalMutant-RPF-Density-corr.txt')
head(RPF_cor_group)

RPF_histogram <- ggplot(data = RPF_cor_group, aes(x = rho, fill = group)) +
  geom_histogram(aes(fill = group, y = after_stat(density)), position = 'identity', 
                 alpha = 0.85, bins = 20) +
  #geom_density(alpha = 0.5) +
  scale_fill_manual(values = c('gold2','cyan4')) + 
  xlab("Spearman's correlation coefficient between two ribo-seq samples") + 
  ggtitle('Ribo-seq') + 
  scale_x_continuous(breaks = seq(0.88, 1, by = 0.02)) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
RPF_histogram


pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/Correlation_Rho.pdf",width=10,height=5.25,useDingbats = F, colormodel = "rgb")
print(RNA_histogram,vp=viewport(.3,.3,x=.25,y=.5))
print(RPF_histogram,vp=viewport(.3,.3,x=.75,y=.5))
dev.off()