
library(ggplot2)#作图 plot
library(ggpubr)#添加显著性标记, Add the significance marker
library(ggsignif)#添加显著性标记, Add the significance marker
library(dplyr)#数据清洗，Data cleaning
library(plyr)#数据清洗，Data cleaning
library(reshape2)#数据清洗，Data cleaning
library(ggthemes)#ggplot所用主题，Themes for ggplot2
library(grid)#分面和嵌合图，facet and Mosaic graph
library(agricolae)#多重比较，Multiple comparisons.
library(readxl)#读入 excel, read excel
library(ggsci)#配色，color scheme
library(showtext)#字体设置, font setting
library(car)#方差齐性检验，homogeneity test of variance, levene test
library(extrafont)#使用系统字体，Using the system fonts
library(sysfonts)#加载系统字体，loading the system fonts
library(Cairo)#抗锯齿,anti-aliasing
library(stringr)#字符串处理.string manipulation
library(graphics)#坐标轴表达式，expression for axis
library(vegan)
library(data.table)
library(forecast)#box-cox数据变换
library(PMCMRplus)
#### 2.

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



#Fig. 2b disease incidence
Disease_incidence<-read_excel("D://R/banana data/Fig.2/2024-3-6/Fig.2.xlsx", sheet="Fig 2b Infect")

Disease_incidence$Treatment <-factor(Disease_incidence$Treatment,levels=c("Y1", "Y2", "Y4", "Y5","Y7", "Y8","Y10", "Y11"))
#Fig 2b statistical analysis #
leveneTest(Infect ~ Treatment, data = Disease_incidence)#p>0.05，则满足方差齐性
shapiro.test(Disease_incidence$Infect)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_Disease_incidence<-aov(data=Disease_incidence,Infect~Treatment)
summary(aov_model_Disease_incidence)
LSD.test(aov_model_Disease_incidence,"Treatment",p.adj = "BH",console=T)
LSD.test(aov_model_Disease_incidence,"Treatment",p.adj = "BH",console=T,group = F)

p <- ggplot(data=Disease_incidence,mapping=aes(x=Treatment,y=Infect,fill=Treatment))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.4,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Disease incidence (%)", x="")+
  scale_fill_manual(values = c('#C5B0D5FF','#9467BDFF','#FFBB78FF','#EEA236FF','#98DF8AFF', '#5CB85CFF','#9EDAE5FF', '#46B8DAFF'))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,110))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  
  mytheme+
  theme(legend.position = "none")
p

ggsave(paste("D://R/banana data/Fig.2/2024-3-6/Disease_incidence-1",".pdf",sep=""),
       device=cairo_pdf,width=60,height=50,dpi = 300,units = "mm")


#Fig. 2c gall index
  Gall_index<-read_excel("D://R/banana data/Fig.2/2024-3-6/Fig.2.xlsx", sheet="Fig 2c Gall")
  Gall_index$Treatment <-factor(Gall_index$Treatment,levels=c("Y1", "Y2", "Y4", "Y5","Y7", "Y8","Y10", "Y11"))
  #Fig 2b statistical analysis #
  leveneTest(Gall ~ Treatment, data = Gall_index)#p>0.05，则满足方差齐性
  shapiro.test(Gall_index$Gall)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
  aov_model_Gall_index<-aov(data=Gall_index,Gall~Treatment)
  summary(aov_model_Gall_index)
  LSD.test(aov_model_Gall_index,"Treatment",p.adj = "BH",console=T)
  LSD.test(aov_model_Gall_index,"Treatment",p.adj = "BH",console=T,group = F)
  
  
  p <- ggplot(data=Gall_index,mapping=aes(x=Treatment,y=Gall,fill=Treatment))+
    geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
    geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.4,stroke = 0.2)+
    stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
                 geom='errorbar', width=0.2,size=0.2)+
    
    labs(y="Gall index", x="")+
    scale_fill_manual(values = c('#C5B0D5FF','#9467BDFF','#FFBB78FF','#EEA236FF','#98DF8AFF', '#5CB85CFF','#9EDAE5FF', '#46B8DAFF'))+
    scale_y_continuous(expand = c(0, 0),limits = c(0,80))+
    
    theme(axis.text=element_text(colour='black',size=8))+
    
    mytheme+
    theme(legend.position = "none")
  p
  ggsave(paste("D://R/banana data/Fig.2/2024-3-6/Gall_index",".pdf",sep=""),
         device=cairo_pdf,width=60,height=50,dpi = 300,units = "mm")
  
  
  
#Fig. 2d soil nematode abundance
Soil_nematode<-read_excel("D://R/banana data/Fig.2/2024-3-6/Fig.2.xlsx", sheet="Fig 2d Abundance")
Soil_nematode$Treatment <-factor(Soil_nematode$Treatment,levels=c("Y1", "Y2", "Y4", "Y5","Y7", "Y8","Y10", "Y11"))
  #Fig 2d statistical analysis #
leveneTest(Abundance ~ Treatment, data = Soil_nematode)#p>0.05，则满足方差齐性
shapiro.test(Soil_nematode$Abundance)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_Soil_nematode<-aov(data=Soil_nematode,Abundance~Treatment)
summary(aov_model_Soil_nematode)
LSD.test(aov_model_Soil_nematode,"Treatment",p.adj = "BH",console=T)
LSD.test(aov_model_Soil_nematode,"Treatment",p.adj = "BH",console=T,group = F)
 #bar -sd
p <- ggplot(data=Soil_nematode,mapping=aes(x=Treatment,y=Abundance,fill=Treatment))+
    geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
    geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.4,stroke = 0.2)+
    stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
                 geom='errorbar', width=0.2,size=0.2)+
    
    labs(y="Soil nematode abundance \n per 100 g dry soi", x="")+
    scale_fill_manual(values = c('#C5B0D5FF','#9467BDFF','#FFBB78FF','#EEA236FF','#98DF8AFF', '#5CB85CFF','#9EDAE5FF', '#46B8DAFF'))+
    scale_y_continuous(expand = c(0, 0),limits = c(0,16000))+
    
    theme(axis.text=element_text(colour='black',size=8))+
    
    mytheme+
    theme(legend.position = "none")
p
ggsave(paste("D://R/banana data/Fig.2/2024-3-6/Soil_nematode",".pdf",sep=""),
         device=cairo_pdf,width=60,height=50,dpi = 300,units = "mm")
  

#Fig. 2e C. elegans abundance
Elegans_abundance<-read_excel("D://R/banana data/Fig.2/2024-3-6/Fig.2.xlsx", sheet="Fig 2e Elegans")
Elegans_abundance$Treatment <-factor(Elegans_abundance$Treatment,levels=c("Y1", "Y2", "Y4", "Y5","Y7", "Y8","Y10", "Y11"))
#Fig 2d statistical analysis #
leveneTest(Celegans ~ Treatment, data = Elegans_abundance)#p>0.05，则满足方差齐性
shapiro.test(Elegans_abundance$Celegans)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_Elegans_abundance<-aov(data=Elegans_abundance,Celegans~Treatment)
summary(aov_model_Elegans_abundance)
LSD.test(aov_model_Elegans_abundance,"Treatment",p.adj = "BH",console=T)
LSD.test(aov_model_Elegans_abundance,"Treatment",p.adj = "BH",console=T,group = F)

#bar-sd
p <- ggplot(data=Elegans_abundance,mapping=aes(x=Treatment,y=Celegans,fill=Treatment))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.4,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Number of offsprings", x="")+
  scale_fill_manual(values = c('#C5B0D5FF','#9467BDFF','#FFBB78FF','#EEA236FF','#98DF8AFF', '#5CB85CFF','#9EDAE5FF', '#46B8DAFF'))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,80))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  
  mytheme+
  theme(legend.position = "none")
p
ggsave(paste("D://R/banana data/Fig.2/2024-3-6/Elegans_abundance-bar-sd",".pdf",sep=""),
       device=cairo_pdf,width=60,height=50,dpi = 300,units = "mm")
