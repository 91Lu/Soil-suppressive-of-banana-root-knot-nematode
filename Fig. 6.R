
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

#Fig. 6a MKB mortality
index <- read_excel("D://R/banana data/Fig.6/致死率/RKN-MKB-2025-6-9/致死率.xlsx")
index$Treatment <-factor(index$Treatment,levels=c("MKB", "Fe_0", "Fe_50"))
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
  labs(y="Relative mortality rate (%)", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF","#9EDAE5FF","#9EDAE5FF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(-7,60))+
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")

p
ggsave(paste("D://R/banana data/Fig.6/致死率/RKN-MKB-2025-6-9/Death rate-3",".pdf",sep=""),
       device=cairo_pdf,width=75,height=70,dpi = 300,units = "mm")

#Fig. 6b bacillibatin mortality
index <- read_excel("D://R/banana data/Fig.6/致死率/RKN-CELE - 2025-6-9/Death rate.xlsx")
index$Treatment <-factor(index$Treatment,levels=c( "R50", "R100", "R500", "50", "100", "500"))
leveneTest(death_rate~Treatment, data =index)#p>0.05，则满足方差齐性
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
  scale_fill_manual(values = c( "#9EDAE5FF","#9EDAE5FF","#9EDAE5FF","#9EDAE5FF",
                                "#9EDAE5FF","#9EDAE5FF","#9EDAE5FF"))+
  scale_y_continuous(limits = c(-3,22))+
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")

p

ggsave(paste("D://R/banana data/Fig.6/致死率/RKN-CELE - 2025-6-9/Death rate-2",".pdf",sep=""),
       device=cairo_pdf,width=80,height=70,dpi = 300,units = "mm")

#Fig. 6c

library(ggplot2)
CI<-read_excel("D://R/banana data/Fig.6/趋化性/RKN-CELE-2025-6-8/RKN-CELE.xlsx")
CI$Treatment <-factor(CI$Treatment,levels=c("CK1", "R10", "R100", "R500","CK2", "10", "100", "500"))
##Fig. 6D CI statistical analysis
leveneTest(Chemotaxis~ Treatment, data =CI)#p>0.05，则满足方差齐性
shapiro.test(CI$Chemotaxis)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_CI<-aov(data=CI,Chemotaxis~Treatment)
summary(aov_model_CI)
LSD.test(aov_model_CI,"Treatment",p.adj = "none",console=T)
LSD.test(aov_model_CI,"Treatment",p.adj = "none",console=T,group = F)

#bar-sd
p <- ggplot(data=CI,mapping=aes(x=Treatment,y=Chemotaxis,fill=Treatment))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.6,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  labs(y="Chemotaxis index", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF","#9EDAE5FF","#9EDAE5FF", "#C7C7C7FF", "#9EDAE5FF","#9EDAE5FF","#9EDAE5FF"))+
  
  scale_y_continuous(limits = c(-0.5,0.8))+
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p

ggsave(paste("D://R/banana data/Fig.6/趋化性/RKN-CELE-2025-6-8/趋化性",".pdf",sep=""),
       device=cairo_pdf,width=80,height=70,dpi = 300,units = "mm")



#Fig. 6e-Sterilized
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
index<-read.csv("D://R/banana data/2025-4铁载体盆栽实验/灭菌/Sterilized soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Bacillibactin"))

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
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,180))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p


ggsave(paste("D://R/banana data/2025-4铁载体盆栽实验/灭菌/Fresh weight",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")


#gall number
library(ggplot2)
index<-read.csv("D://R/banana data/2025-4铁载体盆栽实验/灭菌/Sterilized soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Bacillibactin"))

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
  
  labs(y="Gall numbers per plant", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,120))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p

ggsave(paste("D://R/banana data/2025-4铁载体盆栽实验/灭菌/Gall numbers",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")



#Egg abundance
library(ggplot2)
index<-read.csv("D://R/banana data/2025-4铁载体盆栽实验/灭菌/Sterilized soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Bacillibactin"))

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
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,27))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p


ggsave(paste("D://R/banana data/2025-4铁载体盆栽实验/灭菌/Egg_abundance-1",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")


#root_nematode
library(ggplot2)
index<-read.csv("D://R/banana data/2025-4铁载体盆栽实验/灭菌/Sterilized soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Bacillibactin"))

m="root_nematode"
df1<-index[,c(2,10)]
##正态 qq 图验证数据正态性
qqPlot(lm(index[[m]]~group, data = df1), simulate = TRUE, main = 'QQ Plot', labels = FALSE)

shapiro <- tapply(df1$root_nematode, df1$group, shapiro.test)
shapiro
shapiro$'1'$p.value
shapiro$'2'$p.value

##独立样本的 t 检验
t_test <- t.test(index[[m]]~group, df1, paired = FALSE, alternative = 'two.sided')
t_test
t_test$p.value

data1<-df1%>%group_by(group)%>%summarise_at("root_nematode",funs(mean,sd))
data1
df2<-as.data.frame(data1)

p <- ggplot(data=df1,mapping=aes(x=group,y=root_nematode,fill=group))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=group),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Nematode numbers per 10 g root", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,330))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p


ggsave(paste("D://R/banana data/2025-4铁载体盆栽实验/灭菌/root_nematode abundance",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")



#Fig. 6f-natural soil
#root nematodes abundance
library(ggplot2)
index<-read.csv("D://R/banana data/2025-4铁载体盆栽实验/不灭菌/natural soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Bacillibactin"))

m="root_nematode"
df1<-index[,c(2,6)]
##正态 qq 图验证数据正态性
qqPlot(lm(index[[m]]~group, data = df1), simulate = TRUE, main = 'QQ Plot', labels = FALSE)

shapiro <- tapply(df1$root_nematode, df1$group, shapiro.test)
shapiro
shapiro$'1'$p.value
shapiro$'2'$p.value
##独立样本的 t 检验
t_test <- t.test(index[[m]]~group, df1, paired = FALSE, alternative = 'two.sided')
t_test
t_test$p.value

data1<-df1%>%group_by(group)%>%summarise_at("root_nematode",funs(mean,sd))
data1
df2<-as.data.frame(data1)

p <- ggplot(data=df1,mapping=aes(x=group,y=root_nematode,fill=group))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=group),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Nematode numbers per 10 g root", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,350))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p
ggsave(paste("D://R/banana data/2025-4铁载体盆栽实验/不灭菌/root_nematode",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")




#gall number
library(ggplot2)
index<-read.csv("D://R/banana data/2025-4铁载体盆栽实验/不灭菌/natural soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Bacillibactin"))

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
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,125))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")

p
ggsave(paste("D://R/banana data/2025-4铁载体盆栽实验/不灭菌/Gall numbers",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")


#Egg abundance
library(ggplot2)
index<-read.csv("D://R/banana data/2025-4铁载体盆栽实验/不灭菌/natural soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Bacillibactin"))

m="Egg_abundance"
df1<-index[,c(2,11)]
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
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,30))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p


ggsave(paste("D://R/banana data/2025-4铁载体盆栽实验/不灭菌/Egg_abundance-1",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")

#fresh_weight

library(ggplot2)
index<-read.csv("D://R/banana data/2025-4铁载体盆栽实验/不灭菌/natural soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Bacillibactin"))

m="fresh_weight"
df1<-index[,c(2,10)]
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
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,160))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p


ggsave(paste("D://R/banana data/2025-4铁载体盆栽实验/不灭菌/Fresh weight",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")

