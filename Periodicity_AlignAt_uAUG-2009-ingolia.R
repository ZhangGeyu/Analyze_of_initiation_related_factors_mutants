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
  
##### First uAUG: Rep merge, All sample, 28nt

path = '/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/2-TripletPeriodicity/Periodicity_AlignAt_uAUG-2009-ingolia/First/'

HO_uAUG_Start_count <- read.table(paste0(path,'HO','/','HO','_StartCodonFrame_Length.txt'),header = TRUE)
HO_uAUG_Start_count <- HO_uAUG_Start_count[(HO_uAUG_Start_count$length == 28),]
HO_uAUG_Start_count$count <- HO_uAUG_Start_count$count / sum(HO_uAUG_Start_count$count)
HO_uAUG_Start_count <- HO_uAUG_Start_count[((HO_uAUG_Start_count$StartPos_to_uAUG_start >= -25)&(HO_uAUG_Start_count$StartPos_to_uAUG_start <= 91)),]
HO_uAUG_Start_count$frame <- HO_uAUG_Start_count$StartPos_to_uAUG_start %% 3
head(HO_uAUG_Start_count)

HO_28_length_Start <- ggplot(data = HO_uAUG_Start_count,
                                aes(x = StartPos_to_uAUG_start, 
                                    y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('HO'," 5'end to uAUG")) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.0021),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
HO_28_length_Start

CAF20_uAUG_Start_count <- read.table(paste0(path,'CAF20','/','CAF20','_StartCodonFrame_Length.txt'),header = TRUE)
CAF20_uAUG_Start_count <- CAF20_uAUG_Start_count[(CAF20_uAUG_Start_count$length == 28),]
CAF20_uAUG_Start_count$count <- CAF20_uAUG_Start_count$count / sum(CAF20_uAUG_Start_count$count)
CAF20_uAUG_Start_count <- CAF20_uAUG_Start_count[((CAF20_uAUG_Start_count$StartPos_to_uAUG_start >= -25)&(CAF20_uAUG_Start_count$StartPos_to_uAUG_start <= 91)),]
CAF20_uAUG_Start_count$frame <- CAF20_uAUG_Start_count$StartPos_to_uAUG_start %% 3
head(CAF20_uAUG_Start_count)

CAF20_28_length_Start <- ggplot(data = CAF20_uAUG_Start_count,
                                aes(x = StartPos_to_uAUG_start, 
                                    y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('CAF20')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.0021),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
CAF20_28_length_Start

EAP1_uAUG_Start_count <- read.table(paste0(path,'EAP1','/','EAP1','_StartCodonFrame_Length.txt'),header = TRUE)
EAP1_uAUG_Start_count <- EAP1_uAUG_Start_count[(EAP1_uAUG_Start_count$length == 28),]
EAP1_uAUG_Start_count$count <- EAP1_uAUG_Start_count$count / sum(EAP1_uAUG_Start_count$count)
EAP1_uAUG_Start_count <- EAP1_uAUG_Start_count[((EAP1_uAUG_Start_count$StartPos_to_uAUG_start >= -25)&(EAP1_uAUG_Start_count$StartPos_to_uAUG_start <= 91)),]
EAP1_uAUG_Start_count$frame <- EAP1_uAUG_Start_count$StartPos_to_uAUG_start %% 3
head(EAP1_uAUG_Start_count)

EAP1_28_length_Start <- ggplot(data = EAP1_uAUG_Start_count,
                               aes(x = StartPos_to_uAUG_start, 
                                   y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('EAP1')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.0021),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
EAP1_28_length_Start

eIF4E_uAUG_Start_count <- read.table(paste0(path,'eIF4E','/','eIF4E','_StartCodonFrame_Length.txt'),header = TRUE)
eIF4E_uAUG_Start_count <- eIF4E_uAUG_Start_count[(eIF4E_uAUG_Start_count$length == 28),]
eIF4E_uAUG_Start_count$count <- eIF4E_uAUG_Start_count$count / sum(eIF4E_uAUG_Start_count$count)
eIF4E_uAUG_Start_count <- eIF4E_uAUG_Start_count[((eIF4E_uAUG_Start_count$StartPos_to_uAUG_start >= -25)&(eIF4E_uAUG_Start_count$StartPos_to_uAUG_start <= 91)),]
eIF4E_uAUG_Start_count$frame <- eIF4E_uAUG_Start_count$StartPos_to_uAUG_start %% 3
head(eIF4E_uAUG_Start_count)

eIF4E_28_length_Start <- ggplot(data = eIF4E_uAUG_Start_count,
                                aes(x = StartPos_to_uAUG_start, 
                                    y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('eIF4E')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.0021),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
eIF4E_28_length_Start

eIF4G1_uAUG_Start_count <- read.table(paste0(path,'eIF4G1','/','eIF4G1','_StartCodonFrame_Length.txt'),header = TRUE)
eIF4G1_uAUG_Start_count <- eIF4G1_uAUG_Start_count[(eIF4G1_uAUG_Start_count$length == 28),]
eIF4G1_uAUG_Start_count$count <- eIF4G1_uAUG_Start_count$count / sum(eIF4G1_uAUG_Start_count$count)
eIF4G1_uAUG_Start_count <- eIF4G1_uAUG_Start_count[((eIF4G1_uAUG_Start_count$StartPos_to_uAUG_start >= -25)&(eIF4G1_uAUG_Start_count$StartPos_to_uAUG_start <= 91)),]
eIF4G1_uAUG_Start_count$frame <- eIF4G1_uAUG_Start_count$StartPos_to_uAUG_start %% 3
head(eIF4G1_uAUG_Start_count)

eIF4G1_28_length_Start <- ggplot(data = eIF4G1_uAUG_Start_count,
                                 aes(x = StartPos_to_uAUG_start, 
                                     y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('eIF4G1')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.0021),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
eIF4G1_28_length_Start


PAB1_uAUG_Start_count <- read.table(paste0(path,'PAB1','/','PAB1','_StartCodonFrame_Length.txt'),header = TRUE)
PAB1_uAUG_Start_count <- PAB1_uAUG_Start_count[(PAB1_uAUG_Start_count$length == 28),]
PAB1_uAUG_Start_count$count <- PAB1_uAUG_Start_count$count / sum(PAB1_uAUG_Start_count$count)
PAB1_uAUG_Start_count <- PAB1_uAUG_Start_count[((PAB1_uAUG_Start_count$StartPos_to_uAUG_start >= -25)&(PAB1_uAUG_Start_count$StartPos_to_uAUG_start <= 91)),]
PAB1_uAUG_Start_count$frame <- PAB1_uAUG_Start_count$StartPos_to_uAUG_start %% 3
head(PAB1_uAUG_Start_count)

PAB1_28_length_Start <- ggplot(data = PAB1_uAUG_Start_count,
                               aes(x = StartPos_to_uAUG_start, 
                                   y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('PAB1')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.0021),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
PAB1_28_length_Start

eIF1_uAUG_Start_count <- read.table(paste0(path,'eIF1','/','eIF1','_StartCodonFrame_Length.txt'),header = TRUE)
eIF1_uAUG_Start_count <- eIF1_uAUG_Start_count[(eIF1_uAUG_Start_count$length == 28),]
eIF1_uAUG_Start_count$count <- eIF1_uAUG_Start_count$count / sum(eIF1_uAUG_Start_count$count)
eIF1_uAUG_Start_count <- eIF1_uAUG_Start_count[((eIF1_uAUG_Start_count$StartPos_to_uAUG_start >= -25)&(eIF1_uAUG_Start_count$StartPos_to_uAUG_start <= 91)),]
eIF1_uAUG_Start_count$frame <- eIF1_uAUG_Start_count$StartPos_to_uAUG_start %% 3
head(eIF1_uAUG_Start_count)

eIF1_28_length_Start <- ggplot(data = eIF1_uAUG_Start_count,
                               aes(x = StartPos_to_uAUG_start, 
                                   y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('eIF1')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.0021),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
eIF1_28_length_Start


eIF1A_uAUG_Start_count <- read.table(paste0(path,'eIF1A','/','eIF1A','_StartCodonFrame_Length.txt'),header = TRUE)
eIF1A_uAUG_Start_count <- eIF1A_uAUG_Start_count[(eIF1A_uAUG_Start_count$length == 28),]
eIF1A_uAUG_Start_count$count <- eIF1A_uAUG_Start_count$count / sum(eIF1A_uAUG_Start_count$count)
eIF1A_uAUG_Start_count <- eIF1A_uAUG_Start_count[((eIF1A_uAUG_Start_count$StartPos_to_uAUG_start >= -25)&(eIF1A_uAUG_Start_count$StartPos_to_uAUG_start <= 91)),]
eIF1A_uAUG_Start_count$frame <- eIF1A_uAUG_Start_count$StartPos_to_uAUG_start %% 3
head(eIF1A_uAUG_Start_count)

eIF1A_28_length_Start <- ggplot(data = eIF1A_uAUG_Start_count,
                                aes(x = StartPos_to_uAUG_start, 
                                    y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('eIF1A')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.0021),breaks=seq(0,0.0025,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
eIF1A_28_length_Start


eIF4A_uAUG_Start_count <- read.table(paste0(path,'eIF4A','/','eIF4A','_StartCodonFrame_Length.txt'),header = TRUE)
eIF4A_uAUG_Start_count <- eIF4A_uAUG_Start_count[(eIF4A_uAUG_Start_count$length == 28),]
eIF4A_uAUG_Start_count$count <- eIF4A_uAUG_Start_count$count / sum(eIF4A_uAUG_Start_count$count)
eIF4A_uAUG_Start_count <- eIF4A_uAUG_Start_count[((eIF4A_uAUG_Start_count$StartPos_to_uAUG_start >= -25)&(eIF4A_uAUG_Start_count$StartPos_to_uAUG_start <= 91)),]
eIF4A_uAUG_Start_count$frame <- eIF4A_uAUG_Start_count$StartPos_to_uAUG_start %% 3
head(eIF4A_uAUG_Start_count)

eIF4A_28_length_Start <- ggplot(data = eIF4A_uAUG_Start_count,
                                aes(x = StartPos_to_uAUG_start, 
                                    y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('eIF4A')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.0021),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
eIF4A_28_length_Start


eIF4B_uAUG_Start_count <- read.table(paste0(path,'eIF4B','/','eIF4B','_StartCodonFrame_Length.txt'),header = TRUE)
eIF4B_uAUG_Start_count <- eIF4B_uAUG_Start_count[(eIF4B_uAUG_Start_count$length == 28),]
eIF4B_uAUG_Start_count$count <- eIF4B_uAUG_Start_count$count / sum(eIF4B_uAUG_Start_count$count)
eIF4B_uAUG_Start_count <- eIF4B_uAUG_Start_count[((eIF4B_uAUG_Start_count$StartPos_to_uAUG_start >= -25)&(eIF4B_uAUG_Start_count$StartPos_to_uAUG_start <= 91)),]
eIF4B_uAUG_Start_count$frame <- eIF4B_uAUG_Start_count$StartPos_to_uAUG_start %% 3
head(eIF4B_uAUG_Start_count)

eIF4B_28_length_Start <- ggplot(data = eIF4B_uAUG_Start_count,
                                aes(x = StartPos_to_uAUG_start, 
                                    y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('eIF4B')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.0021),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
eIF4B_28_length_Start



library("grid")
pdf(paste0('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/1-Periodicity/Periodicity_AlignAt_uAUG-2009-ingolia/First_uAUG_Periodicity-RepMerge-28.pdf'),width=4,height=10,useDingbats = F)
print(HO_28_length_Start,vp=viewport(.9,.1,x=.5,y=.94))
print(CAF20_28_length_Start,vp=viewport(.9,.1,x=.5,y=.86))
print(EAP1_28_length_Start,vp=viewport(.9,.1,x=.5,y=.78))
print(eIF4E_28_length_Start,vp=viewport(.9,.1,x=.5,y=.7))
print(eIF4G1_28_length_Start,vp=viewport(.9,.1,x=.5,y=.62))
print(PAB1_28_length_Start,vp=viewport(.9,.1,x=.5,y=.54))
print(eIF1_28_length_Start,vp=viewport(.9,.1,x=.5,y=.46))
print(eIF1A_28_length_Start,vp=viewport(.9,.1,x=.5,y=.38))
print(eIF4A_28_length_Start,vp=viewport(.9,.1,x=.5,y=.30))
print(eIF4B_28_length_Start,vp=viewport(.9,.1,x=.5,y=.22))

dev.off()
#dev.new()

##### Second uAUG: Rep merge, All sample, 28nt

path = '/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/2-TripletPeriodicity/Periodicity_AlignAt_uAUG-2009-ingolia/Second/'

HO_uAUG_Start_count <- read.table(paste0(path,'HO','/','HO','_StartCodonFrame_Length.txt'),header = TRUE)
HO_uAUG_Start_count <- HO_uAUG_Start_count[(HO_uAUG_Start_count$length == 28),]
HO_uAUG_Start_count$count <- HO_uAUG_Start_count$count / sum(HO_uAUG_Start_count$count)
HO_uAUG_Start_count <- HO_uAUG_Start_count[((HO_uAUG_Start_count$StartPos_to_Second_uAUG_start >= -25)&(HO_uAUG_Start_count$StartPos_to_Second_uAUG_start <= 91)),]
HO_uAUG_Start_count$frame <- HO_uAUG_Start_count$StartPos_to_Second_uAUG_start %% 3
head(HO_uAUG_Start_count)

HO_28_length_Start <- ggplot(data = HO_uAUG_Start_count,
                             aes(x = StartPos_to_Second_uAUG_start, 
                                 y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('HO'," 5'end to Second uAUG")) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.002),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
HO_28_length_Start

CAF20_uAUG_Start_count <- read.table(paste0(path,'CAF20','/','CAF20','_StartCodonFrame_Length.txt'),header = TRUE)
CAF20_uAUG_Start_count <- CAF20_uAUG_Start_count[(CAF20_uAUG_Start_count$length == 28),]
CAF20_uAUG_Start_count$count <- CAF20_uAUG_Start_count$count / sum(CAF20_uAUG_Start_count$count)
CAF20_uAUG_Start_count <- CAF20_uAUG_Start_count[((CAF20_uAUG_Start_count$StartPos_to_Second_uAUG_start >= -25)&(CAF20_uAUG_Start_count$StartPos_to_Second_uAUG_start <= 91)),]
CAF20_uAUG_Start_count$frame <- CAF20_uAUG_Start_count$StartPos_to_Second_uAUG_start %% 3
head(CAF20_uAUG_Start_count)

CAF20_28_length_Start <- ggplot(data = CAF20_uAUG_Start_count,
                                aes(x = StartPos_to_Second_uAUG_start, 
                                    y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('CAF20')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.002),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
CAF20_28_length_Start

EAP1_uAUG_Start_count <- read.table(paste0(path,'EAP1','/','EAP1','_StartCodonFrame_Length.txt'),header = TRUE)
EAP1_uAUG_Start_count <- EAP1_uAUG_Start_count[(EAP1_uAUG_Start_count$length == 28),]
EAP1_uAUG_Start_count$count <- EAP1_uAUG_Start_count$count / sum(EAP1_uAUG_Start_count$count)
EAP1_uAUG_Start_count <- EAP1_uAUG_Start_count[((EAP1_uAUG_Start_count$StartPos_to_Second_uAUG_start >= -25)&(EAP1_uAUG_Start_count$StartPos_to_Second_uAUG_start <= 91)),]
EAP1_uAUG_Start_count$frame <- EAP1_uAUG_Start_count$StartPos_to_Second_uAUG_start %% 3
head(EAP1_uAUG_Start_count)

EAP1_28_length_Start <- ggplot(data = EAP1_uAUG_Start_count,
                               aes(x = StartPos_to_Second_uAUG_start, 
                                   y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('EAP1')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.002),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
EAP1_28_length_Start

eIF4E_uAUG_Start_count <- read.table(paste0(path,'eIF4E','/','eIF4E','_StartCodonFrame_Length.txt'),header = TRUE)
eIF4E_uAUG_Start_count <- eIF4E_uAUG_Start_count[(eIF4E_uAUG_Start_count$length == 28),]
eIF4E_uAUG_Start_count$count <- eIF4E_uAUG_Start_count$count / sum(eIF4E_uAUG_Start_count$count)
eIF4E_uAUG_Start_count <- eIF4E_uAUG_Start_count[((eIF4E_uAUG_Start_count$StartPos_to_Second_uAUG_start >= -25)&(eIF4E_uAUG_Start_count$StartPos_to_Second_uAUG_start <= 91)),]
eIF4E_uAUG_Start_count$frame <- eIF4E_uAUG_Start_count$StartPos_to_Second_uAUG_start %% 3
head(eIF4E_uAUG_Start_count)

eIF4E_28_length_Start <- ggplot(data = eIF4E_uAUG_Start_count,
                                aes(x = StartPos_to_Second_uAUG_start, 
                                    y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('eIF4E')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.002),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
eIF4E_28_length_Start

eIF4G1_uAUG_Start_count <- read.table(paste0(path,'eIF4G1','/','eIF4G1','_StartCodonFrame_Length.txt'),header = TRUE)
eIF4G1_uAUG_Start_count <- eIF4G1_uAUG_Start_count[(eIF4G1_uAUG_Start_count$length == 28),]
eIF4G1_uAUG_Start_count$count <- eIF4G1_uAUG_Start_count$count / sum(eIF4G1_uAUG_Start_count$count)
eIF4G1_uAUG_Start_count <- eIF4G1_uAUG_Start_count[((eIF4G1_uAUG_Start_count$StartPos_to_Second_uAUG_start >= -25)&(eIF4G1_uAUG_Start_count$StartPos_to_Second_uAUG_start <= 91)),]
eIF4G1_uAUG_Start_count$frame <- eIF4G1_uAUG_Start_count$StartPos_to_Second_uAUG_start %% 3
head(eIF4G1_uAUG_Start_count)

eIF4G1_28_length_Start <- ggplot(data = eIF4G1_uAUG_Start_count,
                                 aes(x = StartPos_to_Second_uAUG_start, 
                                     y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('eIF4G1')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.002),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
eIF4G1_28_length_Start


PAB1_uAUG_Start_count <- read.table(paste0(path,'PAB1','/','PAB1','_StartCodonFrame_Length.txt'),header = TRUE)
PAB1_uAUG_Start_count <- PAB1_uAUG_Start_count[(PAB1_uAUG_Start_count$length == 28),]
PAB1_uAUG_Start_count$count <- PAB1_uAUG_Start_count$count / sum(PAB1_uAUG_Start_count$count)
PAB1_uAUG_Start_count <- PAB1_uAUG_Start_count[((PAB1_uAUG_Start_count$StartPos_to_Second_uAUG_start >= -25)&(PAB1_uAUG_Start_count$StartPos_to_Second_uAUG_start <= 91)),]
PAB1_uAUG_Start_count$frame <- PAB1_uAUG_Start_count$StartPos_to_Second_uAUG_start %% 3
head(PAB1_uAUG_Start_count)

PAB1_28_length_Start <- ggplot(data = PAB1_uAUG_Start_count,
                               aes(x = StartPos_to_Second_uAUG_start, 
                                   y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('PAB1')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.002),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
PAB1_28_length_Start

eIF1_uAUG_Start_count <- read.table(paste0(path,'eIF1','/','eIF1','_StartCodonFrame_Length.txt'),header = TRUE)
eIF1_uAUG_Start_count <- eIF1_uAUG_Start_count[(eIF1_uAUG_Start_count$length == 28),]
eIF1_uAUG_Start_count$count <- eIF1_uAUG_Start_count$count / sum(eIF1_uAUG_Start_count$count)
eIF1_uAUG_Start_count <- eIF1_uAUG_Start_count[((eIF1_uAUG_Start_count$StartPos_to_Second_uAUG_start >= -25)&(eIF1_uAUG_Start_count$StartPos_to_Second_uAUG_start <= 91)),]
eIF1_uAUG_Start_count$frame <- eIF1_uAUG_Start_count$StartPos_to_Second_uAUG_start %% 3
head(eIF1_uAUG_Start_count)

eIF1_28_length_Start <- ggplot(data = eIF1_uAUG_Start_count,
                               aes(x = StartPos_to_Second_uAUG_start, 
                                   y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('eIF1')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.002),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
eIF1_28_length_Start


eIF1A_uAUG_Start_count <- read.table(paste0(path,'eIF1A','/','eIF1A','_StartCodonFrame_Length.txt'),header = TRUE)
eIF1A_uAUG_Start_count <- eIF1A_uAUG_Start_count[(eIF1A_uAUG_Start_count$length == 28),]
eIF1A_uAUG_Start_count$count <- eIF1A_uAUG_Start_count$count / sum(eIF1A_uAUG_Start_count$count)
eIF1A_uAUG_Start_count <- eIF1A_uAUG_Start_count[((eIF1A_uAUG_Start_count$StartPos_to_Second_uAUG_start >= -25)&(eIF1A_uAUG_Start_count$StartPos_to_Second_uAUG_start <= 91)),]
eIF1A_uAUG_Start_count$frame <- eIF1A_uAUG_Start_count$StartPos_to_Second_uAUG_start %% 3
head(eIF1A_uAUG_Start_count)

eIF1A_28_length_Start <- ggplot(data = eIF1A_uAUG_Start_count,
                                aes(x = StartPos_to_Second_uAUG_start, 
                                    y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('eIF1A')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.002),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
eIF1A_28_length_Start


eIF4A_uAUG_Start_count <- read.table(paste0(path,'eIF4A','/','eIF4A','_StartCodonFrame_Length.txt'),header = TRUE)
eIF4A_uAUG_Start_count <- eIF4A_uAUG_Start_count[(eIF4A_uAUG_Start_count$length == 28),]
eIF4A_uAUG_Start_count$count <- eIF4A_uAUG_Start_count$count / sum(eIF4A_uAUG_Start_count$count)
eIF4A_uAUG_Start_count <- eIF4A_uAUG_Start_count[((eIF4A_uAUG_Start_count$StartPos_to_Second_uAUG_start >= -25)&(eIF4A_uAUG_Start_count$StartPos_to_Second_uAUG_start <= 91)),]
eIF4A_uAUG_Start_count$frame <- eIF4A_uAUG_Start_count$StartPos_to_Second_uAUG_start %% 3
head(eIF4A_uAUG_Start_count)

eIF4A_28_length_Start <- ggplot(data = eIF4A_uAUG_Start_count,
                                aes(x = StartPos_to_Second_uAUG_start, 
                                    y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('eIF4A')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.002),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
eIF4A_28_length_Start


eIF4B_uAUG_Start_count <- read.table(paste0(path,'eIF4B','/','eIF4B','_StartCodonFrame_Length.txt'),header = TRUE)
eIF4B_uAUG_Start_count <- eIF4B_uAUG_Start_count[(eIF4B_uAUG_Start_count$length == 28),]
eIF4B_uAUG_Start_count$count <- eIF4B_uAUG_Start_count$count / sum(eIF4B_uAUG_Start_count$count)
eIF4B_uAUG_Start_count <- eIF4B_uAUG_Start_count[((eIF4B_uAUG_Start_count$StartPos_to_Second_uAUG_start >= -25)&(eIF4B_uAUG_Start_count$StartPos_to_Second_uAUG_start <= 91)),]
eIF4B_uAUG_Start_count$frame <- eIF4B_uAUG_Start_count$StartPos_to_Second_uAUG_start %% 3
head(eIF4B_uAUG_Start_count)

eIF4B_28_length_Start <- ggplot(data = eIF4B_uAUG_Start_count,
                                aes(x = StartPos_to_Second_uAUG_start, 
                                    y = count, fill = factor(frame))) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_color_manual(values = c('red','lightgreen','skyblue')) +
  scale_fill_manual(values = c('red','lightgreen','skyblue')) +
  xlab('') + ylab('') + 
  ggtitle(paste0('eIF4B')) + 
  #geom_vline(xintercept = -12, linetype = "dashed", color = "red", lwd = 0.3) +
  scale_x_continuous(limits = c(-25, 91),breaks=seq(-24,90,6))+
  scale_y_continuous(limits = c(0, 0.002),breaks=seq(0,0.002,0.001))+
  mytheme_violin + 
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white"))
eIF4B_28_length_Start


library("grid")
pdf(paste0('/Dell/Dell15/zhanggy/project/4-Ribo-seq-data/tmp/Ribo-seq_data_final/figure/1-Periodicity/Periodicity_AlignAt_uAUG-2009-ingolia/Second_uAUG_Periodicity-RepMerge-28.pdf'),width=4,height=10,useDingbats = F)
print(HO_28_length_Start,vp=viewport(.9,.1,x=.5,y=.94))
print(CAF20_28_length_Start,vp=viewport(.9,.1,x=.5,y=.86))
print(EAP1_28_length_Start,vp=viewport(.9,.1,x=.5,y=.78))
print(eIF4E_28_length_Start,vp=viewport(.9,.1,x=.5,y=.7))
print(eIF4G1_28_length_Start,vp=viewport(.9,.1,x=.5,y=.62))
print(PAB1_28_length_Start,vp=viewport(.9,.1,x=.5,y=.54))
print(eIF1_28_length_Start,vp=viewport(.9,.1,x=.5,y=.46))
print(eIF1A_28_length_Start,vp=viewport(.9,.1,x=.5,y=.38))
print(eIF4A_28_length_Start,vp=viewport(.9,.1,x=.5,y=.30))
print(eIF4B_28_length_Start,vp=viewport(.9,.1,x=.5,y=.22))

dev.off()
#dev.new()