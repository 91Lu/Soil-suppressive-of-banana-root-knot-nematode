
library(ggplot2)#作图 plot
library(ggpubr)#添加显著性标记, Add the significance marker
library(ggsignif)#添加显著性标记, Add the significance marker
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

#Fig s2a  PCoA
library(tidyverse)
library(vegan)
library(ggthemes)


metadata <- read.table("D://R/banana data/Fig.3/2024-3-7/PCoA/2017/metadata.csv",sep = ",",header=T)
otu <- read.table("D://R/banana data/Fig.3/2024-3-7/PCoA/2017/otu_table.csv",sep = ",",header=T, row.names = 1)
otut <- as.data.frame(t(otu))#转换数据集方向
set.seed(1)
#计算p值和R方
otut.div <- adonis2(otut~ Years, data = metadata, permutations = 999, method="bray")
otut.div
otut_adonis <- paste0("Adonis R2 = ",round(otut.div$R2,2), "; p < ", otut.div$`Pr(>F)`)
otut_dist <- vegdist(otut, method="bray", binary=F)#计算Bray-curtis距离
otut_pcoa <- cmdscale(otut_dist, k=3, eig=T)
otut_pcoa_points <- as.data.frame(otut_pcoa$points)
sum_eig <- sum(otut_pcoa$eig)
eig_percent <- round(otut_pcoa$eig/sum_eig*100,1)#获取维度的结实率
colnames(otut_pcoa_points) <- paste0("PCoA", 1:3)
otut_pcoa_result <- cbind(otut_pcoa_points, metadata)
p <- ggplot(otut_pcoa_result, aes(x=PCoA1, y=PCoA2, color=Years)) +
  scale_color_manual(values = c("#9467BDFF","#EEA236FF", "#5CB85CFF","#46B8DAFF"))+
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
       title = otut_adonis) +
  geom_point(size=2) + 
  theme_test() + geom_hline(yintercept = 0, linetype = "dashed")+geom_vline(xintercept = 0,linetype = "dashed")+
  mytheme

p

ggsave(paste("D://R/banana data/Fig.3/2024-3-7/PCoA/2017/PCoA-1",".pdf",sep=""),
       device=cairo_pdf,width=80,height=60,dpi = 300,units = "mm")



#2017 RKN abundance
RKN_abundance<-read_excel("D://R/banana data/Fig.3/2024-4-7-RKN/RKN.xlsx", sheet="2017")
RKN_abundance$Treatment <-factor(RKN_abundance$Treatment,levels=c("Y2", "Y5","Y8", "Y11"))
#Fig  statistical analysis #
leveneTest(Abundance ~ Treatment, data = RKN_abundance)#p>0.05，则满足方差齐性
shapiro.test(RKN_abundance$Abundance)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_RKN_abundance<-aov(data=RKN_abundance,Abundance~Treatment)
summary(aov_model_RKN_abundance)
LSD.test(aov_model_RKN_abundance,"Treatment",p.adj = "none",console=T)
LSD.test(aov_model_RKN_abundance,"Treatment",p.adj = "none",console=T,group = F)

#bar-sd
p <- ggplot(data=RKN_abundance,mapping=aes(x=Treatment,y=Abundance,fill=Treatment))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.4,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Relative abundance of \n Meloidogyne (%)", x="")+
  scale_fill_manual(values = c("#9467BDFF","#EEA236FF", "#5CB85CFF","#46B8DAFF"))+
  scale_y_continuous(limits = c(-8,40))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  
  mytheme+
  theme(legend.position = "none")
p


ggsave(paste("D://R/banana data/Fig.3/2024-4-7-RKN/2017-RKN-bar-sd-1",".pdf",sep=""),
       device=cairo_pdf,width=45,height=50,dpi = 300,units = "mm")

#bar-se
p <- ggplot(data=RKN_abundance,mapping=aes(x=Treatment,y=Abundance,fill=Treatment))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.4,stroke = 0.2)+
  stat_summary(fun.data=mean_se, 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Relative abundance of \n Meloidogyne (%)", x="")+
  scale_fill_manual(values = c("#9467BDFF","#EEA236FF", "#5CB85CFF","#46B8DAFF"))+
  scale_y_continuous(limits = c(0,40))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  
  mytheme+
  theme(legend.position = "none")
p


ggsave(paste("D://R/banana data/Fig.3/2024-4-7-RKN/2017-RKN-bar-se",".pdf",sep=""),
       device=cairo_pdf,width=45,height=50,dpi = 300,units = "mm")


#Fig. 2e RKN abundance boxplot
p <- ggplot(data=RKN_abundance,mapping=aes(x=Treatment,y=Abundance))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,size=0.2)+  
  geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.4,stroke = 0.2)+
  
  labs(y="Relative abundance of \n Meloidogyne (%)", x="",parse =T)+
  scale_fill_manual(values = c("#9467BDFF","#EEA236FF", "#5CB85CFF","#46B8DAFF"))+
  scale_y_continuous(limits = c(-1,40))+
  mytheme+
  guides(color="none",fill="none")

p 

ggsave(paste("D://R/banana data/Fig.3/2024-4-7-RKN/2017-RKN-boxplot",".pdf",sep=""),
       device=cairo_pdf,width=45,height=50,dpi = 300,units = "mm")
