
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
library(car)#方差齐性检验，homogeneity test of variance, levene test
library(Cairo)#抗锯齿,anti-aliasing
library(stringr)#字符串处理.string manipulation
library(graphics)#坐标轴表达式，expression for axis
library(vegan)
library(data.table)
library(forecast)#box-cox数据变换
library(PMCMRplus)

mytheme<- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#000000"),
                            strip.text = element_text(size=8,hjust = 0.5),
                            plot.title = element_text(size=8,hjust = 0.5),
                            axis.text =element_text(size=8,color = "black"),
                            axis.title =element_text(size=8,color = "black"),
                            legend.text = element_text(size=8,color = "black"),
                            legend.title = element_text(size=8,color = "black"),
                            legend.background = element_blank(),
                            axis.line = element_line(color = "black",size=0.4))#移除整体的边???
index <- read_excel("D://R/banana data/Fig.6/致死率/RKN-LB-2025-6-9/Death rate.xlsx")
index$Treatment <-factor(index$Treatment,levels=c("CK2", "0.1x", "0.5x","1x"))
leveneTest(death_rate~ Treatment, data =index)#p>0.05，则满足方差齐性
shapiro.test(index$death_rate)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_index<-aov(data=index,death_rate~Treatment)
summary(aov_model_index)
LSD.test(aov_model_index,"Treatment",p.adj = "none",console=T)
LSD.test(aov_model_index,"Treatment",p.adj = "none",console=T,group = F)

p <- ggplot(data=index,mapping=aes(x=Treatment,y=death_rate,fill=Treatment))+
  geom_bar(size = 0.25,color="black",fun="mean",position=position_dodge(0.8), stat="summary",width = 0.65)+  
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1),position=position_dodge(0.8), 
               geom='errorbar',width=0.2,size=0.2)+
  geom_jitter (data = index, aes(x=Treatment,y=death_rate,fill=Time),shape=21,
               position=position_jitterdodge(0.3), size=0.6, stroke = 0.2)+
  labs(y="Relative mortality (%)", x="")+
  scale_fill_manual(values = c("#9EDAE5FF","#9EDAE5FF","#9EDAE5FF", "#9EDAE5FF", "#9EDAE5FF","#9EDAE5FF","#9EDAE5FF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(-3,85))+
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p
ggsave(paste("D://R/banana data/Fig.6/致死率/RKN-LB-2025-6-9/Death rate-1",".pdf",sep=""),
       device=cairo_pdf,width=60,height=50,dpi = 300,units = "mm")

