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


#Fig. 3a PCoA
library(ggplot2)
library(vegan)
library(ape)
library(readxl)#读入 excel, read excel
library(tidyverse)
library(vegan)
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

abundance <- read.table("D://R/banana data/Fig.3/2024-3-7/PCoA/abundance.csv",sep = ",",header=T,row.names = 1)
group<-read.table("D://R/banana data/Fig.3/2024-3-7/PCoA/group.csv",sep = ",",header=T, row.names = 1)
group$Years <-factor(group$Years,levels=c("Y1", "Y4", "Y7", "Y10"))
tabundance <- as.data.frame(t(abundance))#转换数据集方向
set.seed(1)

#计算p值和R方
tabundance.div <- adonis2(tabundance~ Years, data = group, permutations = 999, method="bray")
tabundance.div
tabundance_adonis <- paste0("Adonis R2 = ",round(tabundance.div$R2,2), "; p < ", tabundance.div$`Pr(>F)`)

tabundance_dist <- vegdist(tabundance, method="bray", binary=F)#计算Bray-curtis距离
tabundance_pcoa <- cmdscale(tabundance_dist, k=3, eig=T)
tabundance_pcoa_points <- as.data.frame(tabundance_pcoa$points)
sum_eig <- sum(tabundance_pcoa$eig)
eig_percent <- round(tabundance_pcoa$eig/sum_eig*100,1)#获取维度的结实率
colnames(tabundance_pcoa_points) <- paste0("PCoA", 1:3)
tabundance_pcoa_result <- cbind(tabundance_pcoa_points, group)
p <- ggplot(tabundance_pcoa_result, aes(x=PCoA1, y=PCoA2, color=Years)) +
  scale_color_manual(values = c('#C5B0D5FF','#FFBB78FF','#98DF8AFF' ,'#9EDAE5FF'))+
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
       title = tabundance_adonis) +
  geom_point(size=2) + 
  theme_test() + geom_hline(yintercept = 0, linetype = "dashed")+geom_vline(xintercept = 0,linetype = "dashed")+
  mytheme

p

ggsave(paste("D://R/banana data/Fig.3/2024-3-7/PCoA/PCoA-1",".pdf",sep=""),
       device=cairo_pdf,width=80,height=60,dpi = 300,units = "mm")


#Fig. 3b trophic groups  代码有问题。。。
Trophic_levels <- read_excel("D://R/banana data/Fig.3/2024-3-7/Fig.3a Trophic levels.xlsx",sheet="Trophic1")
Trophic_levels<- Trophic_levels %>%
  group_by(Treatment, Trophic) %>%
  summarise(across(c(Abundance), mean))
#Fig 3a statistical analysis #
Trophic_levels2 <- read_excel("D://R/banana data/Fig.3/2024-3-7/Fig.3a Trophic levels.xlsx",sheet="Trophic2")
Trophic_levels2$Treatment <-factor(Trophic_levels2$Treatment,levels=c("Y1", "Y4", "Y7", "Y10"))

#Pp
leveneTest(Pp_abun ~ Treatment, data = Trophic_levels2)#p>0.05，则满足方差齐性
shapiro.test(Trophic_levels2$Pp_abun)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_Trophic_levels2<-aov(data=Trophic_levels2,Pp_abun~Treatment)
summary(aov_model_Trophic_levels2)
LSD.test(aov_model_Trophic_levels2,"Treatment",p.adj = "BH",console=T)
LSD.test(aov_model_Trophic_levels2,"Treatment",p.adj = "BH",console=T,group = F)
#Ba
leveneTest(Ba_abun ~ Treatment, data = Trophic_levels2)#p>0.05，则满足方差齐性
shapiro.test(Trophic_levels2$Ba_abun)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_Trophic_levels2<-aov(data=Trophic_levels2,Ba_abun~Treatment)
summary(aov_model_Trophic_levels2)
LSD.test(aov_model_Trophic_levels2,"Treatment",p.adj = "BH",console=T)
LSD.test(aov_model_Trophic_levels2,"Treatment",p.adj = "BH",console=T,group = F)
#Fu
leveneTest(Fu_abun ~ Treatment, data = Trophic_levels2)#p>0.05，则满足方差齐性
shapiro.test(Trophic_levels2$Fu_abun)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_Trophic_levels2<-aov(data=Trophic_levels2,Fu_abun~Treatment)
summary(aov_model_Trophic_levels2)
LSD.test(aov_model_Trophic_levels2,"Treatment",p.adj = "BH",console=T)
LSD.test(aov_model_Trophic_levels2,"Treatment",p.adj = "BH",console=T,group = F)
#OP
leveneTest(OP_abun ~ Treatment, data = Trophic_levels2)#p>0.05，则满足方差齐性
shapiro.test(Trophic_levels2$OP_abun)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_Trophic_levels2<-aov(data=Trophic_levels2,OP_abun~Treatment)
summary(aov_model_Trophic_levels2)
LSD.test(aov_model_Trophic_levels2,"Treatment",p.adj = "BH",console=T)
LSD.test(aov_model_Trophic_levels2,"Treatment",p.adj = "BH",console=T,group = F)


library(ggplot2)
library(dplyr) 
library(RColorBrewer)
Trophic_levels$Treatment <-factor(Trophic_levels$Treatment,levels=c("Y1", "Y4", "Y7", "Y10"))

Trophic_levels$Trophic <-factor(Trophic_levels$Trophic,levels=c("OP","Fu","Ba","Pp"))

p <- ggplot(data=Trophic_levels,mapping=aes(x=Treatment,y=Abundance,fill=Trophic))+
  geom_bar(size = 0.25,color="black", position='stack', stat="identity",width = 0.65)+  
  scale_fill_manual(values = c('#FDAE6BFF', '#BCBDDCFF', '#A1D99BFF', '#BDBDBDFF'),labels=c("OP","Fu","Ba","Pp"))+ 
  labs(y="Relative abundance of trophic groups (%)", x="")+
  scale_y_continuous(expand = c(0, 0),limits = c(0,100.5))+
  theme_test(base_line_size = 0.75,base_rect_size =0.75)+
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.text=element_text(size=8))+
  theme(legend.key.size = unit (8, "pt"))

p

ggsave(paste("D://R/banana data/Fig.3/2024-3-7/Trophic levels-1",".pdf",sep=""),
       device=cairo_pdf,width=70,height=60,dpi = 300,units = "mm")



#Fig. 3c RKN abundance
RKN_abundance<-read_excel("D://R/banana data/Fig.3/2024-4-7-RKN/RKN.xlsx", sheet="2016")
RKN_abundance$Treatment <-factor(RKN_abundance$Treatment,levels=c("Y1", "Y4","Y7", "Y10"))
#Fig  statistical analysis #
leveneTest(Abundance ~ Treatment, data = RKN_abundance)#p>0.05，则满足方差齐性
shapiro.test(RKN_abundance$Abundance)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_RKN_abundance<-aov(data=RKN_abundance,Abundance~Treatment)
summary(aov_model_RKN_abundance)
LSD.test(aov_model_RKN_abundance,"Treatment",p.adj = "none",console=T)
LSD.test(aov_model_RKN_abundance,"Treatment",p.adj = "none",console=T,group = F)

p <- ggplot(data=RKN_abundance,mapping=aes(x=Treatment,y=Abundance))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,size=0.2)+  
  geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.4,stroke = 0.2)+
  
  labs(y="Relative abundance of \n Meloidogyne (%)", x="",parse =T)+
  scale_fill_manual(values = c('#C5B0D5FF','#FFBB78FF','#98DF8AFF' ,'#9EDAE5FF'))+
  scale_y_continuous(limits = c(0,12))+
  mytheme+
  guides(color="none",fill="none")

p 

#Fig. 3d
rm(list = ls())
library(tidyverse)
library(readxl)
library(ggpubr)
library(ggthemes)

lmce <- read_excel("D://R/banana data/Fig.3/2024-5-13/根结线虫与自由生活线虫线性相关.xlsx")
lmce$Year <- factor(lmce$Year, levels = c("Y1", "Y4", "Y7", "Y10"))

#RKN-Pp
p1 <- ggplot(lmce, aes(x = RKN, y = Pp))+
  geom_point(aes(colour = Year),size = 2)+
  geom_smooth(method = "lm")+
  scale_color_manual(values = c('#C5B0D5FF','#FFBB78FF','#98DF8AFF' ,'#9EDAE5FF'))+
  scale_fill_manual(values = c('#C5B0D5FF','#FFBB78FF','#98DF8AFF' ,'#9EDAE5FF'))+
  
  stat_cor(method = "pearson", size = 3, label.x = 1)+
  theme_test(base_rect_size =1) +
  theme(text = element_text(colour='black'), axis.text=element_text(colour='black',size=12))+
  labs(y = "Pp abundance (%)", x="Meloidogyne abundance (%)")+
  mytheme
p1

ggsave(paste("D://R/banana data/Fig.3/2024-5-13/PP-RKN",".pdf",sep=""),
       device=cairo_pdf,width=80,height=50,dpi = 300,units = "mm")

#RKN-Ba
p2 <- ggplot(lmce, aes(x = RKN, y = Ba))+
  geom_point(aes(colour = Year),size = 2)+
  geom_smooth(method = "lm")+
  scale_color_manual(values = c('#C5B0D5FF','#FFBB78FF','#98DF8AFF' ,'#9EDAE5FF'))+
  scale_fill_manual(values = c('#C5B0D5FF','#FFBB78FF','#98DF8AFF' ,'#9EDAE5FF'))+
  
  stat_cor(method = "pearson", size = 3, label.x = 2)+
  theme_test(base_rect_size =1) +
  theme(text = element_text(colour='black'), axis.text=element_text(colour='black',size=12))+
  labs(y = "Ba abundance (%)", x="Meloidogyne abundance (%)")+
  mytheme
p2

ggsave(paste("D://R/banana data/Fig.3/2024-5-13/Ba-RKN",".pdf",sep=""),
       device=cairo_pdf,width=80,height=50,dpi = 300,units = "mm")

#RKN-Fu
p3 <- ggplot(lmce, aes(x = RKN, y = Fu))+
  geom_point(aes(colour = Year),size = 2)+
  geom_smooth(method = "lm")+
  scale_color_manual(values = c('#C5B0D5FF','#FFBB78FF','#98DF8AFF' ,'#9EDAE5FF'))+
  scale_fill_manual(values = c('#C5B0D5FF','#FFBB78FF','#98DF8AFF' ,'#9EDAE5FF'))+
  
  stat_cor(method = "pearson", size = 3, label.x = 1)+
  theme_test(base_rect_size =1) +
  theme(text = element_text(colour='black'), axis.text=element_text(colour='black',size=12))+
  labs(y = "Fu abundance (%)", x="Meloidogyne abundance (%)")+
  mytheme
p3

ggsave(paste("D://R/banana data/Fig.3/2024-5-13/Fu-RKN",".pdf",sep=""),
       device=cairo_pdf,width=80,height=50,dpi = 300,units = "mm")

#RKN-OP
p4 <- ggplot(lmce, aes(x = RKN, y = OP))+
  geom_point(aes(colour = Year),size = 2)+
  geom_smooth(method = "lm")+
  scale_color_manual(values = c('#C5B0D5FF','#FFBB78FF','#98DF8AFF' ,'#9EDAE5FF'))+
  scale_fill_manual(values = c('#C5B0D5FF','#FFBB78FF','#98DF8AFF' ,'#9EDAE5FF'))+
  
  stat_cor(method = "pearson", size = 3, label.x = 3)+
  theme_test(base_rect_size =1) +
  theme(text = element_text(colour='black'), axis.text=element_text(colour='black',size=12))+
  labs(y = "OP abundance (%)", x="Meloidogyne abundance (%)")+
  mytheme
p4

ggsave(paste("D://R/banana data/Fig.3/2024-5-13/OP-RKN",".pdf",sep=""),
       device=cairo_pdf,width=80,height=50,dpi = 300,units = "mm")

p <- ggarrange(p1,p2,p3,nrow = 1,labels = c("A","B","C"),common.legend = TRUE)
p
ggsave("lm.pdf",width = 12.58, height = 6.70, units = "in")



#Fig. 3e Shannon
Shannon<-read_excel("D://R/banana data/Fig.3/2024-7-20/Fig.3.xlsx", sheet="Fig 3c Shannon")

Shannon$Treatment <-factor(Shannon$Treatment,levels=c("Y1", "Y4", "Y7", "Y10"))
##Fig. 3c statistical analysis
leveneTest(H ~ Treatment, data = Shannon)#p>0.05，则满足方差齐性
shapiro.test(Shannon$H)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_Shannon<-aov(data=Shannon,H~Treatment)
summary(aov_model_Shannon)
LSD.test(aov_model_Shannon,"Treatment",p.adj = "none",console=T)
LSD.test(aov_model_Shannon,"Treatment",p.adj = "none",console=T,group = F)


#Fig. 3c Shannon boxplot
p <- ggplot(data=Shannon,mapping=aes(x=Treatment,y=H))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,size=0.2)+  
  geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.6,stroke = 0.2)+
  
  labs(y="Shannon index", x="",parse =T)+
  scale_fill_manual(values = c('#C5B0D5FF','#FFBB78FF','#98DF8AFF' ,'#9EDAE5FF'))+
  scale_y_continuous(limits = c(2,3.5))+
  mytheme+
  guides(color="none",fill="none")
p 

ggsave(paste("D://R/banana data/Fig.3/2024-7-20/Shannon-boxplot",".pdf",sep=""),
       device=cairo_pdf,width=40,height=45,dpi = 300,units = "mm")




#Fig. 3c MI
MI_Index<-read_excel("D://R/banana data/Fig.3/2024-7-20/Fig.3.xlsx", sheet="Fig 3c MI")

MI_Index$Treatment <-factor(MI_Index$Treatment,levels=c("Y1", "Y4", "Y7", "Y10"))
##Fig. 3c_MI_Index statistical analysis
leveneTest(MI ~ Treatment, data = MI_Index)#p>0.05，则满足方差齐性
shapiro.test(MI_Index$MI)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_MI_Index<-aov(data=MI_Index,MI~Treatment)
summary(aov_model_MI_Index)
LSD.test(aov_model_MI_Index,"Treatment",p.adj = "none",console=T)
LSD.test(aov_model_MI_Index,"Treatment",p.adj = "none",console=T,group = F)


#Fig. 3c MI_Index boxplot
p <- ggplot(data=MI_Index,mapping=aes(x=Treatment,y=MI))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,size=0.2)+  
  geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.6,stroke = 0.2)+
  
  labs(y="MI", x="",parse =T)+
  scale_fill_manual(values = c('#C5B0D5FF','#FFBB78FF','#98DF8AFF' ,'#9EDAE5FF'))+
  scale_y_continuous(limits = c(1,3))+
  mytheme+
  guides(color="none",fill="none")
p 

ggsave(paste("D://R/banana data/Fig.3/2024-7-20/MI-boxplot",".pdf",sep=""),
       device=cairo_pdf,width=40,height=45,dpi = 300,units = "mm")


#Fig. 3c WI
WI_Index<-read_excel("D://R/banana data/Fig.3/2024-7-20/Fig.3.xlsx", sheet="Fig 3c WI")

WI_Index$Treatment <-factor(WI_Index$Treatment,levels=c("Y1", "Y4", "Y7", "Y10"))
##Fig. 3c_WI_Index statistical analysis
leveneTest(WI ~ Treatment, data = WI_Index)#p>0.05，则满足方差齐性
shapiro.test(WI_Index$WI)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_WI_Index<-aov(data=WI_Index,WI~Treatment)
summary(aov_model_WI_Index)
LSD.test(aov_model_WI_Index,"Treatment",p.adj = "BH",console=T)
LSD.test(aov_model_WI_Index,"Treatment",p.adj = "BH",console=T,group = F)


#Fig. 3c WI_Index boxplot
p <- ggplot(data=WI_Index,mapping=aes(x=Treatment,y=WI))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,size=0.2)+  
  geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.6,stroke = 0.2)+
  
  labs(y="WI", x="",parse =T)+
  scale_fill_manual(values = c('#C5B0D5FF','#FFBB78FF','#98DF8AFF' ,'#9EDAE5FF'))+
  scale_y_continuous(limits = c(0,9))+
  mytheme+
  guides(color="none",fill="none")
p 

ggsave(paste("D://R/banana data/Fig.3/2024-7-20/WI-boxplot",".pdf",sep=""),
       device=cairo_pdf,width=40,height=45,dpi = 300,units = "mm")


#Fig. 3 PPI/MI
PPI_MI_Index<-read_excel("D://R/banana data/Fig.3/2024-7-20/Fig.3.xlsx", sheet="Fig 3 PPI-MI")

PPI_MI_Index$Treatment <-factor(PPI_MI_Index$Treatment,levels=c("Y1", "Y4", "Y7", "Y10"))
##Fig. S3_ppi_mi_Index statistical analysis
leveneTest(PPI_MI ~ Treatment, data = PPI_MI_Index)#p>0.05，则满足方差齐性
shapiro.test(PPI_MI_Index$PPI_MI)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_PPI_MI_Index<-aov(data=PPI_MI_Index,PPI_MI~Treatment)
summary(aov_model_PPI_MI_Index)
LSD.test(aov_model_PPI_MI_Index,"Treatment",p.adj = "none",console=T)
LSD.test(aov_model_PPI_MI_Index,"Treatment",p.adj = "none",console=T,group = F)


#Fig. S3 PPI_MI_Index boxplot
p <- ggplot(data=PPI_MI_Index,mapping=aes(x=Treatment,y=PPI_MI))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,size=0.2)+  
  geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.6,stroke = 0.2)+
  
  labs(y="PPI/MI", x="",parse =T)+
  scale_fill_manual(values = c('#C5B0D5FF','#FFBB78FF','#98DF8AFF' ,'#9EDAE5FF'))+
  scale_y_continuous(limits = c(0.5,2))+
  mytheme+
  guides(color="none",fill="none")
p 

ggsave(paste("D://R/banana data/Fig.3/2024-7-20/PPI_MI-boxplot",".pdf",sep=""),
       device=cairo_pdf,width=40,height=45,dpi = 300,units = "mm")
