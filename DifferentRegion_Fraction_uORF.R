
detach('package:Rmisc')
detach('package:plyr')

##### RPF Fraction Cumulative curve——uAUG 'HO','eIF1','eIF1A','eIF4A','eIF4B'

Fraction_Four <- Fraction[(Fraction$sample == 'HO')|(Fraction$sample == 'eIF1')|
                            (Fraction$sample == 'eIF1A'),]
head(Fraction_Four)

Fraction_cumulative <- Fraction_Four %>%
  group_by(sample) %>%
  mutate(cumulative = ecdf(log2_Fraction)(log2_Fraction))
head(Fraction_cumulative)

Fraction_CumulativeCurve <- ggplot(Fraction_cumulative, 
                                   aes(log2_Fraction, cumulative, color = sample)) +
  geom_step() + 
  scale_color_manual(values = c("black", "red", "#0000FF")) + 
  labs(title = 'RPF Fraction',
       x = 'RPF Fraction [log2]',
       y = "Cumulative Probability") + 
  facet_wrap(~region, ncol = 1,scales = 'free') + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
Fraction_CumulativeCurve

Fraction_cumulative$Fraction <- Fraction_cumulative$Fraction + 1
Fraction_boxplot <- ggplot(Fraction_cumulative, 
                           aes(sample, log2(Fraction), color = sample)) +
  geom_boxplot(outlier.size = 0.5, fill = 'white', notch=TRUE) + 
  coord_cartesian(ylim = c(0, 1.4)) + 
  scale_color_manual(values = c("black", "red", "#0000FF")) + 
  labs(title = 'RPF Fraction',
       x = 'samples',
       y = '(RPF Fraction + 1) [log2]') + 
  geom_signif(comparisons = list(c('HO','eIF1'),c('HO','eIF1A'),c('HO','eIF4A'),c('HO','eIF4B')), 
              step_increase = 0.08,
              map_signif_level = F, color = "black", textsize = 2.5, tip_length = 0.005,y_position = 1) +
  facet_wrap(~region, ncol = 1, scales = 'free') + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"),
        axis.text.x = element_blank())
Fraction_boxplot

Fraction_cumulative_uORF <- Fraction_Four %>%
  group_by(sample,uORF) %>%
  mutate(cumulative = ecdf(log2_Fraction)(log2_Fraction))

Fraction_CumulativeCurve_uORF <- ggplot(Fraction_cumulative_uORF, 
                                        aes(log2_Fraction, cumulative, color = sample)) +
  geom_step() + 
  scale_color_manual(values = c("black", "red", "#0000FF")) + 
  labs(title = 'RPF Fraction',
       x = 'RPF Fraction [log2]',
       y = "Cumulative Probability") + 
  facet_wrap(uORF~region, ncol = 1,scales = 'free') + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
Fraction_CumulativeCurve_uORF

Fraction_cumulative_uORF$Fraction <- Fraction_cumulative_uORF$Fraction + 1
Fraction_boxplot_uORF <- ggplot(Fraction_cumulative_uORF, 
                                aes(sample, log2(Fraction), color = sample)) +
  geom_boxplot(outlier.size = 0.5, fill = 'white', notch=TRUE) + 
  coord_cartesian(ylim = c(0, 1.4)) + 
  scale_color_manual(values = c("black", "red", "#0000FF")) + 
  labs(title = 'RPF Fraction',
       x = 'samples',
       y = '(RPF Fraction + 1) [log2]') + 
  geom_signif(comparisons = list(c('HO','eIF1'),c('HO','eIF1A'),c('HO','eIF4A'),c('HO','eIF4B')), 
              step_increase = 0.08,
              map_signif_level = F, color = "black", textsize = 2.5, tip_length = 0.005,y_position = 1) +
  facet_wrap(uORF~region, ncol = 1,scales = 'free') + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"),
        axis.text.x = element_blank())
Fraction_boxplot_uORF

Fraction_cumulative_uORF <- Fraction_Four %>%
  group_by(sample,uORF) %>%
  mutate(cumulative = ecdf(log2_Fraction)(log2_Fraction))

Fraction_CumulativeCurve_sample <- ggplot(Fraction_cumulative_uORF, 
                                          aes(log2_Fraction, cumulative, color = uORF)) +
  geom_step() + 
  scale_color_manual(values = c("black", "red")) + 
  labs(title = 'RPF Fraction',
       x = 'RPF Fraction [log2]',
       y = "Cumulative Probability",
       color = "Group") + 
  facet_wrap(sample~region, scales = 'free',ncol = 1) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"))
Fraction_CumulativeCurve_sample

Fraction_cumulative_uORF$Fraction <- Fraction_cumulative_uORF$Fraction + 1
Fraction_boxplot_sample <- ggplot(Fraction_cumulative_uORF, 
                                  aes(uORF, log2(Fraction), color = uORF)) +
  geom_boxplot(outlier.size = 0.5, fill = 'white', notch=TRUE) + 
  coord_cartesian(ylim = c(0, 1.4)) + 
  scale_color_manual(values = c('black','red')) + 
  labs(title = 'RPF Fraction',
       x = 'have uORF ?',
       y = '(RPF Fraction + 1) [log2]') + 
  geom_signif(comparisons = list(c('No','Yes')), 
              step_increase = 0.08,
              map_signif_level = F, color = "black", textsize = 2.5, tip_length = 0.005,y_position = 1) +
  facet_wrap(sample~region, ncol = 1,scales = 'free') + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"),
        axis.text.x = element_blank())
Fraction_boxplot_sample

library("grid")
pdf('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/5-Region_ReadsFraction-Gene/uORF/CumulativeCurve-Fraction-Four.pdf',width=10,height=10,useDingbats = F)
print(Fraction_CumulativeCurve,vp=viewport(.35,.25,x=.3,y=.8))
print(Fraction_CumulativeCurve_uORF,vp=viewport(.35,.5,x=.3,y=.3))
print(Fraction_CumulativeCurve_sample,vp=viewport(.25,.95,x=.75,y=.5))
dev.off()

library("grid")
pdf('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/5-Region_ReadsFraction-Gene/uORF/Boxplot-Fraction-Four.pdf',width=10,height=10,useDingbats = F)
print(Fraction_boxplot,vp=viewport(.35,.27,x=.3,y=.78))
print(Fraction_boxplot_uORF,vp=viewport(.35,.55,x=.3,y=.3))
print(Fraction_boxplot_sample,vp=viewport(.25,.95,x=.72,y=.5))
dev.off()
