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

##### Reads Fraction in different region

### gene have uORF

AllSample_RegionReadsCount <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/5-ReadsFraction-Region/Ribo/AllSample_Reads_Fraction_28nt_uORF.txt')
AllSample_RegionReadsCount$sample <- factor(AllSample_RegionReadsCount$sample, levels = c('eIF4B','eIF4A','eIF1A','eIF1','PAB1','eIF4G1','eIF4E','EAP1','CAF20','HO'))
AllSample_RegionReadsCount$region <- factor(AllSample_RegionReadsCount$region, levels = c('5UTR','CDS','3UTR'))
#AllSample_RegionReadsCount <- AllSample_RegionReadsCount[AllSample_RegionReadsCount$region == '5UTR',]
head(AllSample_RegionReadsCount)

RegionReadsCount_bar <- ggplot(data = AllSample_RegionReadsCount,
                               aes(x = sample, 
                                   y = fraction, fill = region)) + 
  stat_summary(fun = match.fun(mean), geom = "bar", fun.args = list(mult = 1), 
               position = position_dodge(0.6), width = 0.6, alpha = 0.5) +
  stat_summary(aes(colour = region),fun.data = mean_se, 
               position = position_dodge(0.5), geom = "errorbar", width = 0.1) + 
  geom_jitter(aes(colour = region),shape = 21, size = 0.6, 
              position = position_dodge(0.5)) +
  scale_color_manual(values = c('orange','grey','purple')) +
  scale_fill_manual(values = c('orange','grey','purple')) +
  xlab("sample") + ylab('Fraction') + 
  ggtitle('Ribo-seq: 28nt; genes have uORF') + 
  facet_wrap(~region,scales = 'free') + 
  mytheme_violin + coord_flip()+
  theme(legend.position="none")+
  theme(strip.background = element_rect(color = "white", fill = "white"))
RegionReadsCount_bar

RegionReadsCount_bar_5UTR <- ggplot(data = AllSample_RegionReadsCount[AllSample_RegionReadsCount$region == '5UTR',],
                                    aes(x = sample, y = fraction, fill = region)) + 
  stat_summary(fun = match.fun(mean), geom = "bar", fun.args = list(mult = 1), 
               position = position_dodge(0.6), width = 0.6, alpha = 0.5) +
  stat_summary(aes(colour = region),fun.data = mean_se, 
               position = position_dodge(0.5), geom = "errorbar", width = 0.1) + 
  geom_jitter(aes(colour = region),shape = 21, size = 0.6, 
              position = position_dodge(0.5)) +
  scale_color_manual(values = c('orange','grey','purple')) +
  scale_fill_manual(values = c('orange','grey','purple')) +
  xlab("sample") + ylab('Fraction') + 
  ggtitle('Ribo-seq: 28nt; genes have uORF') + 
  geom_signif(test = "wilcox.test",
              comparisons =list(c('HO','CAF20'),c('HO','EAP1'),c('HO','eIF4E'),
                                c('HO','eIF4G1'),c('HO','PAB1'),c('HO','eIF1'),
                                c('HO','eIF1A'),c('HO','eIF4A'),c('HO','eIF4B')),
              map_signif_level = F, color = "black",step_increase = 0.1,
              textsize = 2.5, tip_length = 0.005, y_position = 0.02) + 
  facet_wrap(~region,scales = 'free') + 
  mytheme_violin + coord_flip()+
  theme(legend.position="none")+
  theme(strip.background = element_rect(color = "white", fill = "white"))
RegionReadsCount_bar_5UTR

pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/5-Region_ReadsFraction/RegionReadsFraction_uORF.pdf",width=6,height=8,useDingbats = F)
print(RegionReadsCount_bar,vp=viewport(1,.4,x=.5,y=.75))
print(RegionReadsCount_bar_5UTR,vp=viewport(.5,.4,x=.3,y=.25))
dev.off()
#dev.new()

### gene no uORF


AllSample_RegionReadsCount <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/5-ReadsFraction-Region/Ribo/AllSample_Reads_Fraction_28nt_NOuORF.txt')
AllSample_RegionReadsCount$sample <- factor(AllSample_RegionReadsCount$sample, levels = c('eIF4B','eIF4A','eIF1A','eIF1','PAB1','eIF4G1','eIF4E','EAP1','CAF20','HO'))
AllSample_RegionReadsCount$region <- factor(AllSample_RegionReadsCount$region, levels = c('5UTR','CDS','3UTR'))
#AllSample_RegionReadsCount <- AllSample_RegionReadsCount[AllSample_RegionReadsCount$region == '5UTR',]
head(AllSample_RegionReadsCount)

RegionReadsCount_bar <- ggplot(data = AllSample_RegionReadsCount,
                               aes(x = sample, 
                                   y = fraction, fill = region)) + 
  stat_summary(fun = match.fun(mean), geom = "bar", fun.args = list(mult = 1), 
               position = position_dodge(0.6), width = 0.6, alpha = 0.5) +
  stat_summary(aes(colour = region),fun.data = mean_se, 
               position = position_dodge(0.5), geom = "errorbar", width = 0.1) + 
  geom_jitter(aes(colour = region),shape = 21, size = 0.6, 
              position = position_dodge(0.5)) +
  scale_color_manual(values = c('orange','grey','purple')) +
  scale_fill_manual(values = c('orange','grey','purple')) +
  xlab("sample") + ylab('Fraction') + 
  ggtitle('Ribo-seq: 28nt; genes without uORF') + 
  facet_wrap(~region,scales = 'free') + 
  mytheme_violin + coord_flip()+
  theme(legend.position="none")+
  theme(strip.background = element_rect(color = "white", fill = "white"))
RegionReadsCount_bar

RegionReadsCount_bar_5UTR <- ggplot(data = AllSample_RegionReadsCount[AllSample_RegionReadsCount$region == '5UTR',],
                                    aes(x = sample, y = fraction, fill = region)) + 
  stat_summary(fun = match.fun(mean), geom = "bar", fun.args = list(mult = 1), 
               position = position_dodge(0.6), width = 0.6, alpha = 0.5) +
  stat_summary(aes(colour = region),fun.data = mean_se, 
               position = position_dodge(0.5), geom = "errorbar", width = 0.1) + 
  geom_jitter(aes(colour = region),shape = 21, size = 0.6, 
              position = position_dodge(0.5)) +
  scale_color_manual(values = c('orange','grey','purple')) +
  scale_fill_manual(values = c('orange','grey','purple')) +
  xlab("sample") + ylab('Fraction') + 
  ggtitle('Ribo-seq: 28nt; genes without uORF') + 
  geom_signif(test = "wilcox.test",
              comparisons =list(c('HO','CAF20'),c('HO','EAP1'),c('HO','eIF4E'),
                                c('HO','eIF4G1'),c('HO','PAB1'),c('HO','eIF1'),
                                c('HO','eIF1A'),c('HO','eIF4A'),c('HO','eIF4B')),
              map_signif_level = F, color = "black",step_increase = 0.1,
              textsize = 2.5, tip_length = 0.005, y_position = 0.02) + 
  facet_wrap(~region,scales = 'free') + 
  mytheme_violin + coord_flip()+
  theme(legend.position="none")+
  theme(strip.background = element_rect(color = "white", fill = "white"))
RegionReadsCount_bar_5UTR

pdf("/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/5-Region_ReadsFraction/RegionReadsFraction_NOuORF.pdf",width=6,height=8,useDingbats = F)
print(RegionReadsCount_bar,vp=viewport(1,.4,x=.5,y=.75))
print(RegionReadsCount_bar_5UTR,vp=viewport(.5,.4,x=.3,y=.25))
dev.off()
#dev.new()