
library(ggthemes)
mytheme<- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#000000"),
                            
                            strip.text = element_text(size=8,hjust = 0.5),
                            plot.title = element_text(size=8,hjust = 0.5),
                            axis.text =element_text(size=8,color = "black"),
                            axis.title =element_text(size=8,color = "black"),
                            legend.text = element_text(size=8,color = "black"),
                            legend.title = element_text(size=8,color = "black"),
                            legend.background = element_blank(),
                            axis.line = element_line(color = "black",size=0.4))#移除整体的边???

rm(list = ls())
library(tidyverse)
library(readxl)
library(ggpubr)
lmce <- read_excel("D://R/banana data/Fig.2/2024-3-6/C.elegans线性相关.xlsx")
lmce$Year <- factor(lmce$Year, levels = c("Y1", "Y2", "Y4", "Y5", "Y7", "Y8", "Y10", "Y11"))

#C.elegans-Infect
p1 <- ggplot(lmce, aes(x = Celegans, y = Infect))+
    geom_point(aes(colour = Year),size = 1)+
    geom_smooth(method = "lm")+
    scale_color_manual(values = c('#C5B0D5FF','#9467BDFF','#FFBB78FF','#EEA236FF','#98DF8AFF', '#5CB85CFF','#9EDAE5FF', '#46B8DAFF'))+
    scale_fill_manual(values = c('#C5B0D5FF','#9467BDFF','#FFBB78FF','#EEA236FF','#98DF8AFF', '#5CB85CFF','#9EDAE5FF', '#46B8DAFF'))+
    
    stat_cor(method = "pearson", size = 3, label.x = 20)+
    theme_test(base_rect_size =1) +
    theme(text = element_text(colour='black'), axis.text=element_text(colour='black',size=12))+
    labs(y = "Disease incidence", x="Number of offsprings")+
  mytheme
p1

ggsave(paste("D://R/banana data/Fig.2/2024-3-6/C.elegans-Disease",".pdf",sep=""),
       device=cairo_pdf,width=80,height=50,dpi = 300,units = "mm")


#C.elegans-gall
p2 <- ggplot(lmce, aes(x = Celegans, y = Gall))+
  geom_point(aes(colour = Year),size = 1)+
  geom_smooth(method = "lm")+
  scale_color_manual(values = c('#C5B0D5FF','#9467BDFF','#FFBB78FF','#EEA236FF','#98DF8AFF', '#5CB85CFF','#9EDAE5FF', '#46B8DAFF'))+
  scale_fill_manual(values = c('#C5B0D5FF','#9467BDFF','#FFBB78FF','#EEA236FF','#98DF8AFF', '#5CB85CFF','#9EDAE5FF', '#46B8DAFF'))+
  
  stat_cor(method = "pearson", size = 3, label.x = 20)+
  theme_test(base_rect_size =1) +
  theme(text = element_text(colour='black'), axis.text=element_text(colour='black',size=12))+
  labs(y = "Gall index", x="Number of offsprings")+
  mytheme
p2

ggsave(paste("D://R/banana data/Fig.2/2024-3-6/C.elegans-Gall",".pdf",sep=""),
       device=cairo_pdf,width=80,height=50,dpi = 300,units = "mm")


#C.elegans-Abundance
p3 <- ggplot(lmce, aes(x = Celegans, y = Abundance))+
  geom_point(aes(colour = Year),size = 1)+
  geom_smooth(method = "lm")+
  scale_color_manual(values = c('#C5B0D5FF','#9467BDFF','#FFBB78FF','#EEA236FF','#98DF8AFF', '#5CB85CFF','#9EDAE5FF', '#46B8DAFF'))+
  scale_fill_manual(values = c('#C5B0D5FF','#9467BDFF','#FFBB78FF','#EEA236FF','#98DF8AFF', '#5CB85CFF','#9EDAE5FF', '#46B8DAFF'))+
  
  stat_cor(method = "pearson", size = 3, label.x = 20)+
  theme_test(base_rect_size =1) +
  theme(text = element_text(colour='black'), axis.text=element_text(colour='black',size=12))+
  labs(y = "Total soil nematode abundance \n per 100 g dry soil", x="Number of offsprings")+
  mytheme
p3

ggsave(paste("D://R/banana data/Fig.2/2024-3-6/C.elegans-Abundance",".pdf",sep=""),
       device=cairo_pdf,width=80,height=50,dpi = 300,units = "mm")


p <- ggarrange(p1,p2,p3,nrow = 1,labels = c("A","B","C"),common.legend = TRUE)
p
ggsave("lm.pdf",width = 12.58, height = 6.70, units = "in")
