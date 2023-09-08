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


##### Binding Affinity Over/Under represented ~ RNA/RPF/Ribosome Density foldchange

BindAffi_fold <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/6-BindingAffinity/TotalMutant-RNA_RPF_TE-Density-Fold-BindAffi.txt')
BindAffi_fold$group <- factor(BindAffi_fold$group, levels = c('Over-Represented','Not_Significant','Under-Represented'))
head(BindAffi_fold)

BindAffi_TEfold_corr <- ggplot(data = BindAffi_fold, 
                               aes(x = log2(RibosomeDensity_foldchange), y = logFC)) +
  coord_cartesian(xlim = c(-3, 3),
                  ylim = c(-3, 3)) + 
  geom_point(size = 0.2, alpha = 0.5) +
  xlab('Ribosome Density foldchange [log2]') + ylab('Binding Affinity [log2]') + 
  mytheme_violin + 
  facet_wrap(~sample, scales = 'free',ncol = 5) + 
  stat_cor(label.x = -3, label.y = 3, method = "spearman",size = 2) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
BindAffi_TEfold_corr


pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/6-BindingAffinity/BindingAffinity_TEfoldchange.pdf",width=8,height=3,useDingbats = F, colormodel = "rgb")
print(BindAffi_TEfold_corr,vp=viewport(.95,.6,x=.5,y=.5))
dev.off()


##### Binding Affinity ~ PolyA length and nonA

BindAffi_PolyA <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/6-BindingAffinity/PolyA_length_nonA_BindAffi.txt')
BindAffi_PolyA_PAB1 <- BindAffi_PolyA[BindAffi_PolyA$sample == 'PAB1',]
head(BindAffi_PolyA_PAB1)

BindAffi_PolyA_length_corr <- ggplot(data = BindAffi_PolyA_PAB1, 
                                    aes(x = poly_length_average, y = logFC)) +
  coord_cartesian(ylim = c(-2, 2.1)) + 
  geom_point(size = 0.2, alpha = 0.5) +
  geom_text(label = paste0("n = ", length(BindAffi_PolyA_PAB1$poly_length_average)),
            x = 4.8, y = 1.5,
            size = 3, fontface = "plain") + 
  xlab('polyA length [log2]') + ylab('PAB1 Binding Affinity') + 
  mytheme_violin + 
  facet_wrap(~sample, scales = 'free',ncol = 5) + 
  stat_cor(label.x = 26, label.y = 2, method = "spearman",size = 2.8) + 
  scale_x_continuous(breaks = seq(25, 70, by = 5)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"),
        strip.background = element_rect(color = "white", fill = "white"))
BindAffi_PolyA_length_corr

pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/10-20230831/BindingAffinity_PolyA.pdf",width=8,height=5,useDingbats = F, colormodel = "rgb")
print(BindAffi_PolyA_length_corr,vp=viewport(.29,.5,x=.25,y=.75))
dev.off()
