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


#Fig 9a. Root_dry_weight

library(ggplot2)
index<-read.csv("D://R/banana data/2025-4铁载体盆栽实验/灭菌/Sterilized soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Bacillibactin"))

m="root_dry_weight"
df1<-index[,c(2,11)]
##正态 qq 图验证数据正态性
qqPlot(lm(index[[m]]~group, data = df1), simulate = TRUE, main = 'QQ Plot', labels = FALSE)

shapiro <- tapply(df1$root_dry_weight, df1$group, shapiro.test)
shapiro
shapiro$'1'$p.value
shapiro$'2'$p.value

##独立样本的 t 检验
t_test <- t.test(index[[m]]~group, df1, paired = FALSE, alternative = 'two.sided')
t_test
t_test$p.value

data1<-df1%>%group_by(group)%>%summarise_at("root_dry_weight",funs(mean,sd))
data1
df2<-as.data.frame(data1)

p <- ggplot(data=df1,mapping=aes(x=group,y=root_dry_weight,fill=group))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=group),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Root dry weight (g/plant)", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,12))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p


ggsave(paste("D://R/banana data/2025-4铁载体盆栽实验/灭菌/Root dry weight",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")




#Fig 9a Dry weight of shoot
library(ggplot2)
index<-read.csv("D://R/banana data/2025-4铁载体盆栽实验/灭菌/Sterilized soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Bacillibactin"))

m="dry_weight"
df1<-index[,c(2,5)]
##正态 qq 图验证数据正态性
qqPlot(lm(index[[m]]~group, data = df1), simulate = TRUE, main = 'QQ Plot', labels = FALSE)

shapiro <- tapply(df1$dry_weight, df1$group, shapiro.test)
shapiro
shapiro$'1'$p.value
shapiro$'2'$p.value

##独立样本的 t 检验
t_test <- t.test(index[[m]]~group, df1, paired = FALSE, alternative = 'two.sided')
t_test
t_test$p.value

data1<-df1%>%group_by(group)%>%summarise_at("dry_weight",funs(mean,sd))
data1
df2<-as.data.frame(data1)

p <- ggplot(data=df1,mapping=aes(x=group,y=dry_weight,fill=group))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=group),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Dry weight of shoot (g/plant)", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,16))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p

ggsave(paste("D://R/banana data/2025-4铁载体盆栽实验/灭菌/ Shoot dry weight",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")



#Fig 9b Root-dry_weight

library(ggplot2)
index<-read.csv("D://R/banana data/2025-4铁载体盆栽实验/不灭菌/natural soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Bacillibactin"))

m="root_dry_weight"
df1<-index[,c(2,12)]
##正态 qq 图验证数据正态性
qqPlot(lm(index[[m]]~group, data = df1), simulate = TRUE, main = 'QQ Plot', labels = FALSE)

shapiro <- tapply(df1$root_dry_weight, df1$group, shapiro.test)
shapiro
shapiro$'1'$p.value
shapiro$'2'$p.value

##独立样本的 t 检验
t_test <- t.test(index[[m]]~group, df1, paired = FALSE, alternative = 'two.sided')
t_test
t_test$p.value

data1<-df1%>%group_by(group)%>%summarise_at("root_dry_weight",funs(mean,sd))
data1
df2<-as.data.frame(data1)

p <- ggplot(data=df1,mapping=aes(x=group,y=root_dry_weight,fill=group))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=group),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Root dry weight (g/plant)", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,12))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p


ggsave(paste("D://R/banana data/2025-4铁载体盆栽实验/不灭菌/root-dry- weight",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")




#Fig 9b shoot dry weight
library(ggplot2)
index<-read.csv("D://R/banana data/2025-4铁载体盆栽实验/不灭菌/natural soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Bacillibactin"))

m="dry_weight"
df1<-index[,c(2,5)]
##正态 qq 图验证数据正态性
qqPlot(lm(index[[m]]~group, data = df1), simulate = TRUE, main = 'QQ Plot', labels = FALSE)

shapiro <- tapply(df1$dry_weight, df1$group, shapiro.test)
shapiro

##独立样本的 t 检验
t_test <- t.test(index[[m]]~group, df1, paired = FALSE, alternative = 'two.sided')
t_test
t_test$p.value

data1<-df1%>%group_by(group)%>%summarise_at("dry_weight",funs(mean,sd))
data1
df2<-as.data.frame(data1)

p <- ggplot(data=df1,mapping=aes(x=group,y=dry_weight,fill=group))+
  geom_bar(size = 0.25,color="black",fun="mean", stat="summary",width = 0.65)+  
  geom_point(aes(fill=group),shape=21,position = position_jitterdodge(1),size=0.8,stroke = 0.2)+
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1), 
               geom='errorbar', width=0.2,size=0.2)+
  
  labs(y="Dry weight of shoot (g/plant)", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,15))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p
ggsave(paste("D://R/banana data/2025-4铁载体盆栽实验/不灭菌/Dry weight",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")

