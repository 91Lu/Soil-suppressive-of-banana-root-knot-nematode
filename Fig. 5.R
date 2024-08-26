library(ggthemes)#ggplot所用主题，Themes for ggplot2
library(ggplot2)
mytheme<- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#000000"),
                            
                            strip.text = element_text(size=8,hjust = 0.5),
                            plot.title = element_text(size=8,hjust = 0.5),
                            axis.text =element_text(size=8,color = "black"),
                            axis.title =element_text(size=8,color = "black"),
                            legend.text = element_text(size=8,color = "black"),
                            legend.title = element_text(size=8,color = "black"),
                            legend.background = element_blank(),
                            axis.line = element_line(color = "black",size=0.4))#移除整体的边???


#Fig.5a
rm(list=ls())
library(tidyverse)
library(psych)
# library(devtools)
# devtools::install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
#连作前后期差异OTU表
ind <- read.delim("ind.txt")
#线虫指标表
index <- read.delim("index.txt")
#实验设计表
metadata <- read.delim("metadata.txt")
#细菌丰度表
otutab <- read.delim("otutabB.txt")
#细菌分类表
tax <- read.delim("taxonomyB.txt")
#细菌网络模块表
module <- read.delim("ext_module.txt")
otutab[,-1] <- otutab[,-1]/colSums(otutab[,-1])
#提取ind中color为Y02Y05的数据,提取module中module为2的列并与indY2Y5按OTUID列合并，取两个数据框的交集
indY2Y5 <- ind %>% filter(color == "Y02Y05") %>% inner_join(module, by = "OTUID") %>% filter(module == 1)
#同理提取差异Y8Y11
indY8Y11 <- ind %>% filter(color == "Y08Y11") %>% inner_join(module, by = "OTUID") %>% filter(module == 2)
#对otutab除第一列外的所有列按行求和并删除多余列
otutab_sum <- otutab %>% mutate(sum = rowSums(.[,-1])) %>% select(OTUID, sum)
#将indY2Y5与otutab_sum按OTUID列合并，取交集，并按sum列降序排列
indY2Y5_sum <- indY2Y5 %>% inner_join(otutab_sum, by = "OTUID") 
indY2Y5 <- indY2Y5_sum %>% arrange(desc(sum))
indY2Y5 <- indY2Y5[1:10,]
indY8Y11_sum <- indY8Y11 %>% inner_join(otutab_sum, by = "OTUID")
indY8Y11 <- indY8Y11_sum %>% arrange(desc(sum))
indY8Y11 <- indY8Y11[1:10,]
#合并两个数据框
indY2Y5Y8Y11 <- rbind(indY2Y5, indY8Y11)
#提取otutab中在indY2Y5Y8Y11中的OTUID列存在的行
otutab_ind <- otutab %>% filter(OTUID %in% indY2Y5Y8Y11$OTUID)
#将otutab_ind中的行的顺序按indY2Y5Y8Y11中OTUID的顺序排列
otutab_ind <- otutab_ind[match(indY2Y5Y8Y11$OTUID, otutab_ind$OTUID),]
#将第一列变为行名并删除第一列
rownames(otutab_ind) <- otutab_ind[,1]
otutab_ind <- otutab_ind[,-1]
otutab_ind <- as.data.frame(t(otutab_ind))
otutab_ind$SampleID <- rownames(otutab_ind)
#将otutab_ind与index按SampleID列合并
index <- otutab_ind %>% inner_join(index, by = "SampleID")
#将index与metadata按SampleID列合并
index <- index %>% inner_join(metadata, by = "SampleID")

inf_rate <- corr.test(index$inf_rate, index[,1:20], method="spearman",adjust="BH",minlength=5)
inf_rate$r
inf_rate$p
meloidogyne <- corr.test(index$meloidogyne, index[,1:20], method="spearman",adjust="BH",minlength=5)
meloidogyne$r
meloidogyne$p
gall <- corr.test(index$gall, index[,1:20], method="spearman",adjust="BH",minlength=5)
gall$r
gall$p
abundance <- corr.test(index$abundance, index[,1:20], method="spearman",adjust="BH",minlength=5)
abundance$r
abundance$p
Free_living_nematode <- corr.test(index$Free_living_nematode, index[,1:20], method="spearman",adjust="BH",minlength=5)
Free_living_nematode$r
Free_living_nematode$p
Herbivorous_nematodes <- corr.test(index$Herbivorous_nematodes, index[,1:20], method="spearman",adjust="BH",minlength=5)
Herbivorous_nematodes$r
Herbivorous_nematodes$p
## 8.3.2 split ##
cor <- as.data.frame(t(rbind(inf_rate$r, inf_rate$p, meloidogyne$r,meloidogyne$p,abundance$r,abundance$p)))
colnames(cor) <- c("inf_rate_r","inf_rate_p","meloidogyne_r","meloidogyne_p","abundance_r","abundance_p")
cor$ASV <- rownames(cor)
# 绘制r矩阵热图
cor_r <- as.data.frame(t(rbind(inf_rate$r, meloidogyne$r,abundance$r)))
colnames(cor_r) <- c("Disease_incidence","meloidogyne_abundance","Soil_nematode_abundance")
max(cor_r)
min(cor_r)
col_cor = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#4d9221", "white", "#c51b7d"))
row_order <- c("ASV_6972","ASV_6738","ASV_5065","ASV_1467","ASV_5963","ASV_8041","ASV_4833","ASV_7704","ASV_8525","ASV_7026","ASV_7216","ASV_5729","ASV_6061","ASV_5244","ASV_4472","ASV_3920","ASV_6213","ASV_5425","ASV_4459","ASV_7438")
length(row_order)
row_order
mat_s = cor_r[row_order,]
mat_s = as.matrix(mat_s)
cor_r = cor_r[row_order,]
cor_r = as.matrix(cor_r)
Heatmap(cor_r,col = col_cor,
        cluster_columns = F,cluster_rows = F,
        row_order = row_order,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", mat_s[i, j]), x, y, gp = gpar(fontsize = 10))}
)
# 绘制p矩阵热图
cor_p <- as.data.frame(t(rbind(inf_rate$p, meloidogyne$p,abundance$p)))
colnames(cor_p) <- c("Disease_incidence","meloidogyne_abundance","Soil_nematode_abundance")
max(cor_r)
min(cor_r)
col_cor = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#4d9221", "white", "#c51b7d"))
row_order <- c("ASV_6972","ASV_6738","ASV_5065","ASV_1467","ASV_5963","ASV_8041","ASV_4833","ASV_7704","ASV_8525","ASV_7026","ASV_7216","ASV_5729","ASV_6061","ASV_5244","ASV_4472","ASV_3920","ASV_6213","ASV_5425","ASV_4459","ASV_7438")
length(row_order)
row_order
mat_s = cor_p[row_order,]
mat_s = as.matrix(mat_s)
cor_r = cor_r[row_order,]
cor_r = as.matrix(cor_r)
Heatmap(cor_r,col = col_cor,
        cluster_columns = F,cluster_rows = F,
        row_order = row_order,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.5f", mat_s[i, j]), x, y, gp = gpar(fontsize = 10))}
)



#Fig.5b
library(ggplot2)
ASV1467_Abundance<-read_excel("D://R/banana data/Fig.5/2024-3-12/Fig.5.xlsx",sheet="ASV1467 abundance")
ASV1467_Abundance$Treatment <-factor(ASV1467_Abundance$Treatment,levels=c("Y2", "Y5", "Y8", "Y11"))
##Fig. 5c ASV1467_Abundance statistical analysis
leveneTest(Abundance~ Treatment, data =ASV1467_Abundance)#p>0.05，则满足方差齐性
shapiro.test(ASV1467_Abundance$Abundance)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_ASV1467_Abundance<-aov(data=ASV1467_Abundance,Abundance~Treatment)
summary(aov_model_ASV1467_Abundance)
LSD.test(aov_model_ASV1467_Abundance,"Treatment",p.adj = "none",console=T)
LSD.test(aov_model_ASV1467_Abundance,"Treatment",p.adj = "none",console=T,group = F)


#Fig. 5c ASV1467_Abundance boxplot
p <- ggplot(data=ASV1467_Abundance,mapping=aes(x=Treatment,y=Abundance))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,size=0.2)+  
  geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.6,stroke = 0.2)+
  
  labs(y="Relative abundance\n of ASV1467 (%)", x="",parse =T)+
  scale_fill_manual(values = c("#9467BDFF","#EEA236FF", "#5CB85CFF","#46B8DAFF"))+
  scale_y_continuous(limits = c(0,0.8))+
  mytheme+
  guides(color="none",fill="none")
p 

ggsave(paste("D://R/banana data/Fig.5/2024-3-12/ASV1467_Abundance-boxplot",".pdf",sep=""),
       device=cairo_pdf,width=40,height=40,dpi = 300,units = "mm")


#Fig. 5e
library(reshape2) 
##正态 qq 图验证数据正态性
library(car)
library(ggplot2)    #ggplot2 作图
library(dplyr)
library(ggthemes)#ggplot所用主题，Themes for ggplot2

mytheme<- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#000000"),
                            
                            strip.text = element_text(size=8,hjust = 0.5),
                            plot.title = element_text(size=8,hjust = 0.5),
                            axis.text =element_text(size=8,color = "black"),
                            axis.title =element_text(size=8,color = "black"),
                            legend.text = element_text(size=8,color = "black"),
                            legend.title = element_text(size=8,color = "black"),
                            legend.background = element_blank(),
                            axis.line = element_line(color = "black",size=0.4))#移除整体的边???

#fresh_weight

library(ggplot2)
index<-read.csv("D://R/banana data/Fig.5/2025-4-22/灭菌/Sterilized soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Y11.1"))

m="fresh_weight"
df1<-index[,c(2,6)]
##正态 qq 图验证数据正态性
qqPlot(lm(index[[m]]~group, data = df1), simulate = TRUE, main = 'QQ Plot', labels = FALSE)

shapiro <- tapply(df1$fresh_weight, df1$group, shapiro.test)
shapiro
shapiro$'1'$p.value
shapiro$'2'$p.value

##独立样本的 t 检验
t_test <- t.test(index[[m]]~group, df1, paired = FALSE, alternative = 'two.sided')
t_test
t_test$p.value

data1<-df1%>%group_by(group)%>%summarise_at("fresh_weight",funs(mean,sd))
data1
df2<-as.data.frame(data1)

p <- ggplot(data=df1,mapping=aes(x=group,y=fresh_weight,fill=group))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=group),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Fresh weight of shoot (g/plant)", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#98DF8AFF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,50))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p


ggsave(paste("D://R/banana data/Fig.5/2025-4-22/灭菌/Fresh weight",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")


#soil nematodes
library(ggplot2)
index<-read.csv("D://R/banana data/Fig.5/2025-4-22/灭菌/Sterilized soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Y11.1"))

m="soil_nematode"
df1<-index[,c(2,7)]
##正态 qq 图验证数据正态性
qqPlot(lm(index[[m]]~group, data = df1), simulate = TRUE, main = 'QQ Plot', labels = FALSE)

shapiro <- tapply(df1$soil_nematode, df1$group, shapiro.test)
shapiro
shapiro$'1'$p.value
shapiro$'2'$p.value

##独立样本的 t 检验
t_test <- t.test(index[[m]]~group, df1, paired = FALSE, alternative = 'two.sided')
t_test
t_test$p.value

data1<-df1%>%group_by(group)%>%summarise_at("soil_nematode",funs(mean,sd))
data1
df2<-as.data.frame(data1)

p <- ggplot(data=df1,mapping=aes(x=group,y=soil_nematode,fill=group))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=group),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Nematode abundance per 100 g dry soil", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#98DF8AFF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,5200))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p

ggsave(paste("D://R/banana data/Fig.5/2025-4-22/灭菌/soil_nematode-1",".pdf",sep=""),
       device=cairo_pdf,width=43,height=55,dpi = 300,units = "mm")

#gall number
library(ggplot2)
index<-read.csv("D://R/banana data/Fig.5/2025-4-22/灭菌/Sterilized soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Y11.1"))

m="gall_number"
df1<-index[,c(2,8)]
##正态 qq 图验证数据正态性
qqPlot(lm(index[[m]]~group, data = df1), simulate = TRUE, main = 'QQ Plot', labels = FALSE)

shapiro <- tapply(df1$gall_number, df1$group, shapiro.test)
shapiro
shapiro$'1'$p.value
shapiro$'2'$p.value

##独立样本的 t 检验
t_test <- t.test(index[[m]]~group, df1, paired = FALSE, alternative = 'two.sided')
t_test
t_test$p.value

data1<-df1%>%group_by(group)%>%summarise_at("gall_number",funs(mean,sd))
data1
df2<-as.data.frame(data1)

p <- ggplot(data=df1,mapping=aes(x=group,y=gall_number,fill=group))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=group),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Gall numbers per 10 g root", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#98DF8AFF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,50))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p

ggsave(paste("D://R/banana data/Fig.5/2025-4-22/灭菌/Gall numbers-Y",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")



#Egg abundance
library(ggplot2)
index<-read.csv("D://R/banana data/Fig.5/2025-4-22/灭菌/Sterilized soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Y11.1"))

m="Egg_abundance"
df1<-index[,c(2,9)]
##正态 qq 图验证数据正态性
qqPlot(lm(index[[m]]~group, data = df1), simulate = TRUE, main = 'QQ Plot', labels = FALSE)

shapiro <- tapply(df1$Egg_abundance, df1$group, shapiro.test)
shapiro
shapiro$'1'$p.value
shapiro$'2'$p.value

##独立样本的 t 检验
t_test <- t.test(index[[m]]~group, df1, paired = FALSE, alternative = 'two.sided')
t_test
t_test$p.value

data1<-df1%>%group_by(group)%>%summarise_at("Egg_abundance",funs(mean,sd))
data1
df2<-as.data.frame(data1)

p <- ggplot(data=df1,mapping=aes(x=group,y=Egg_abundance,fill=group))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=group),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Egg mass numbers per 10 g root", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#98DF8AFF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,25))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p


ggsave(paste("D://R/banana data/Fig.5/2025-4-22/灭菌/Egg_abundance",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")



#Fig. 5f

#soil nematodes
library(ggplot2)
index<-read.csv("D://R/banana data/Fig.5/2025-4-22/不灭菌/natural soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Y11.1"))

m="soil_nematode"
df1<-index[,c(2,7)]
##正态 qq 图验证数据正态性
qqPlot(lm(index[[m]]~group, data = df1), simulate = TRUE, main = 'QQ Plot', labels = FALSE)

shapiro <- tapply(df1$soil_nematode, df1$group, shapiro.test)
shapiro
shapiro$'1'$p.value
shapiro$'2'$p.value
##独立样本的 t 检验
t_test <- t.test(index[[m]]~group, df1, paired = FALSE, alternative = 'two.sided')
t_test
t_test$p.value

data1<-df1%>%group_by(group)%>%summarise_at("soil_nematode",funs(mean,sd))
data1
df2<-as.data.frame(data1)

p <- ggplot(data=df1,mapping=aes(x=group,y=soil_nematode,fill=group))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=group),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Nematode abundance per 100 g dry soil", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#98DF8AFF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,6200))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p

ggsave(paste("D://R/banana data/Fig.5/2025-4-22/不灭菌/soil_nematode-1",".pdf",sep=""),
       device=cairo_pdf,width=44,height=55,dpi = 300,units = "mm")

#gall number
library(ggplot2)
index<-read.csv("D://R/banana data/Fig.5/2025-4-22/不灭菌/natural soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Y11.1"))

m="gall_number"
df1<-index[,c(2,8)]
##正态 qq 图验证数据正态性
qqPlot(lm(index[[m]]~group, data = df1), simulate = TRUE, main = 'QQ Plot', labels = FALSE)

shapiro <- tapply(df1$gall_number, df1$group, shapiro.test)
shapiro
shapiro$'1'$p.value
shapiro$'2'$p.value
##独立样本的 t 检验
t_test <- t.test(index[[m]]~group, df1, paired = FALSE, alternative = 'two.sided')
t_test
t_test$p.value

data1<-df1%>%group_by(group)%>%summarise_at("gall_number",funs(mean,sd))
data1
df2<-as.data.frame(data1)

p <- ggplot(data=df1,mapping=aes(x=group,y=gall_number,fill=group))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=group),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Gall numbers per 10 g root", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#98DF8AFF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,75))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")

p
ggsave(paste("D://R/banana data/Fig.5/2025-4-22/不灭菌/Gall numbers",".pdf",sep=""),
       device=cairo_pdf,width=45,height=55,dpi = 300,units = "mm")


#RKN abundance
library(ggplot2)
index<-read.csv("D://R/banana data/Fig.5/2025-4-22/不灭菌/natural soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Y11.1"))

m="RKN_abundance"
df1<-index[,c(2,9)]
##正态 qq 图验证数据正态性
qqPlot(lm(index[[m]]~group, data = df1), simulate = TRUE, main = 'QQ Plot', labels = FALSE)

shapiro <- tapply(df1$RKN_abundance, df1$group, shapiro.test)
shapiro
shapiro$'1'$p.value
shapiro$'2'$p.value
##独立样本的 t 检验
t_test <- t.test(index[[m]]~group, df1, paired = FALSE, alternative = 'two.sided')
t_test
t_test$p.value

data1<-df1%>%group_by(group)%>%summarise_at("RKN_abundance",funs(mean,sd))
data1
df2<-as.data.frame(data1)

p <- ggplot(data=df1,mapping=aes(x=group,y=RKN_abundance,fill=group))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=group),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Meloidogyne density (105 copies/g soil)", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#98DF8AFF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,3))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p

ggsave(paste("D://R/banana data/Fig.5/2025-4-22/不灭菌/RKN_abundance",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")

#Fresh weight
library(ggplot2)
index<-read.csv("D://R/banana data/Fig.5/2025-4-22/不灭菌/natural soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Y11.1"))

m="Fresh_weight"
df1<-index[,c(2,10)]
##正态 qq 图验证数据正态性
qqPlot(lm(index[[m]]~group, data = df1), simulate = TRUE, main = 'QQ Plot', labels = FALSE)

shapiro <- tapply(df1$Fresh_weight, df1$group, shapiro.test)
shapiro
shapiro$'1'$p.value
shapiro$'2'$p.value
##独立样本的 t 检验
t_test <- t.test(index[[m]]~group, df1, paired = FALSE, alternative = 'two.sided')
t_test
t_test$p.value

data1<-df1%>%group_by(group)%>%summarise_at("Fresh_weight",funs(mean,sd))
data1
df2<-as.data.frame(data1)

p <- ggplot(data=df1,mapping=aes(x=group,y=Fresh_weight,fill=group))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=group),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Fresh weight of shoot (g/plant)", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#98DF8AFF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,300))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p

ggsave(paste("D://R/banana data/Fig.5/2025-4-22/不灭菌/Fresh_weight",".pdf",sep=""),
       device=cairo_pdf,width=42,height=55,dpi = 300,units = "mm")

