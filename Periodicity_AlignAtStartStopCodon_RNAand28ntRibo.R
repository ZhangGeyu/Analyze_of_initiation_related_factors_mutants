library(ggplot2)

mytheme_violin <-  theme_bw()+#theme(legend.position="none")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black", fill=NA, linewidth=.5))+
  theme(plot.title=element_text(size=5))+
  # theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))+
  theme(axis.text.x = element_text(colour="black", size=5))+
  theme(axis.text.y = element_text(colour="black", size=5))+
  theme(axis.title=element_text(size=5))+
  #theme(plot.title=element_text(size=rel(1.5), lineheight=.9))+
  theme(legend.title=element_text(size=5))+
  theme(legend.text=element_text(size=5))+
  theme(axis.ticks = element_line(size = 0.2))

##### Relative to Start Codon

path = '/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/2-TripletPeriodicity/Periodicity_AlignAtStartStopCodon/Ribo_RNA_Merge/'
sample = 'PAB1'

### rep merge

Start_count <- read.table(paste0(path,sample,'/',sample,'_StartFrame.txt'),header = TRUE)
Start_count$Ribo_counts <- Start_count$Ribo_counts / sum(Start_count$Ribo_counts)
Start_count$RNA_counts <- Start_count$RNA_counts / sum(Start_count$RNA_counts)
Start_count <- Start_count[((Start_count$StartPos_to_StartCodon >= -21)&(Start_count$StartPos_to_StartCodon <= 15)),]
ymax <- max(Start_count$Ribo_counts)
head(Start_count)

Start_count_bar <- ggplot(data = Start_count) + 
  #geom_bar(aes(x = StartPos_to_StartCodon,y = Ribo_counts,fill = factor(frame)),
  #         position = 'dodge', stat = 'identity') + 
  geom_line(aes(x = StartPos_to_StartCodon,y = Ribo_counts),colour = 'blue',linewidth = 0.4) + 
  geom_line(aes(x = StartPos_to_StartCodon,y = RNA_counts),colour = 'darkgreen',linewidth = 0.3) + 
  xlab('Position Relative to StartCodon') + ylab('Count normalized') + 
  coord_cartesian(ylim = c(0, ymax)) + 
  scale_color_manual(values = c('blue','grey','grey')) +
  scale_fill_manual(values = c('blue','grey','grey')) +
  ggtitle(paste0(sample," 5'end to start codon")) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.2, colour = "black")) +
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-22, 16),breaks=seq(-24,16,6)) + 
  scale_y_continuous(limits = c(0,ymax),breaks=seq(0,ymax,0.0015))
Start_count_bar

Stop_count <- read.table(paste0(path,sample,'/',sample,'_StopFrame.txt'),header = TRUE)
Stop_count$Ribo_counts <- Stop_count$Ribo_counts / sum(Stop_count$Ribo_counts)
Stop_count$RNA_counts <- Stop_count$RNA_counts / sum(Stop_count$RNA_counts)
Stop_count <- Stop_count[((Stop_count$StartPos_to_StopCodon >= -36)&(Stop_count$StartPos_to_StopCodon <=0)),]
head(Stop_count)

Stop_count_bar <- ggplot(data = Stop_count) + 
  #geom_bar(aes(x = StartPos_to_StopCodon,y = Ribo_counts,fill = factor(frame)),
  #         position = 'dodge', stat = 'identity') + 
  geom_line(aes(x = StartPos_to_StopCodon,y = Ribo_counts),colour = 'blue',linewidth = 0.4) + 
  geom_line(aes(x = StartPos_to_StopCodon,y = RNA_counts),colour = 'darkgreen',linewidth = 0.3) + 
  xlab('Position Relative to StopCodon') + ylab('Count normalized') + 
  coord_cartesian(ylim = c(0, ymax)) + 
  scale_color_manual(values = c('blue','grey','grey')) +
  scale_fill_manual(values = c('blue','grey','grey')) +
  ggtitle(paste0(sample," 5'end to stop codon")) + 
  mytheme_violin + 
  theme(strip.background = element_rect(color = "white", fill = "white"),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.2, colour = "black")) +
  theme(legend.position="none")+
  scale_x_continuous(limits = c(-37, 1),breaks=seq(-36,1,6)) + 
  scale_y_continuous(limits = c(0,ymax),breaks=seq(0,ymax,0.0015),position = "right")
Stop_count_bar


library("grid")
pdf(paste0('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/1-Periodicity/Periodicity_AlignAtStartStopCodon_RNAand28ntRibo/',sample,'-Periodicity.pdf'),width=4,height=1.5,useDingbats = F)
print(Start_count_bar,vp=viewport(.35,.5,x=.3,y=.5))
print(Stop_count_bar,vp=viewport(.35,.5,x=.65,y=.5))

dev.off()
#dev.new()

##### Frame Fraction

FrameFraction <- read.table('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/2-TripletPeriodicity/Frame_Fraction/AllSample_RegionFrame.txt')
FrameFraction$sample <- factor(FrameFraction$sample, levels = c('HO','eIF1','eIF1A','eIF4A','eIF4B','eIF4E','eIF4G1','PAB1','EAP1','CAF20'))
FrameFraction$frame <- factor(FrameFraction$frame, levels = c('0','1','2'))
head(FrameFraction)

FrameFraction_bar <- ggplot(data = FrameFraction,
                            aes(x = frame, 
                                y = fraction, fill = Type)) + 
  stat_summary(fun = match.fun(mean), geom = "bar", fun.args = list(mult = 1), 
               position = position_dodge(0.6), width = 0.6, alpha = 0.7) +
  stat_summary(aes(colour = Type),fun.data = mean_se, 
               position = position_dodge(0.5), geom = "errorbar", width = 0.1) + 
  geom_jitter(aes(colour = Type),shape = 21, size = 0.6, 
              position = position_dodge(0.5)) +
  coord_cartesian(ylim = c(0, 1)) + 
  scale_color_manual(values = c('blue','grey')) +
  scale_fill_manual(values = c('blue','grey')) +
  xlab("frame") + ylab('Fraction') + 
  ggtitle('Ribo-seq: 28nt and RNA: all length') + 
  facet_wrap(~sample,scales = 'free') + 
  scale_y_continuous(limits = c(0,1),breaks=seq(0,1,0.2)) + 
  mytheme_violin +
  theme(strip.background = element_rect(color = "white", fill = "white"),
        axis.text.x = element_text(angle=0, hjust=0.5))
FrameFraction_bar

library("grid")
pdf('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/1-Periodicity/Periodicity_AlignAtStartStopCodon_RNAand28ntRibo/FrameFraction.pdf',width=7,height=5,useDingbats = F)
print(FrameFraction_bar,vp=viewport(1,1,x=.5,y=.5))

dev.off()
#dev.new()