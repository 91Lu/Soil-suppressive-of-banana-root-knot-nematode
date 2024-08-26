#Fig 6b sterilized soil
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




index<-read.csv("D://R/出站报告/功能菌盆栽实验/灭菌/Sterilized soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Y11.1"))


#dry weight
library(ggplot2)
index<-read.csv("D://R/出站报告/功能菌盆栽实验/灭菌/Sterilized soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Y11.1"))

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
  scale_fill_manual(values = c("#C7C7C7FF", "#98DF8AFF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,4))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p

ggsave(paste("D://R/出站报告/功能菌盆栽实验/灭菌/Dry weight",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")




#soil nematodes
library(ggplot2)
index<-read.csv("D://R/出站报告/功能菌盆栽实验/灭菌/Sterilized soil.csv")#,row.names=1
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
  
  labs(y="Soil nematode abundance per 100 g dry soil", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#98DF8AFF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,5000))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p

ggsave(paste("D://R/出站报告/功能菌盆栽实验/灭菌/soil_nematode",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")

#gall number
library(ggplot2)
index<-read.csv("D://R/出站报告/功能菌盆栽实验/灭菌/Sterilized soil.csv")#,row.names=1
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

ggsave(paste("D://R/出站报告/功能菌盆栽实验/灭菌/Gall numbers-Y",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")



#Egg abundance
library(ggplot2)
index<-read.csv("D://R/出站报告/功能菌盆栽实验/灭菌/Sterilized soil.csv")#,row.names=1
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


ggsave(paste("D://R/出站报告/功能菌盆栽实验/灭菌/Egg_abundance",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")



# Fig 6c natural soil

index<-read.csv("D://R/出站报告/功能菌盆栽实验/不灭菌/natural soil.csv",stringsAsFactors = FALSE)#row.names=1,,
index$group <-factor(index$group,levels=c("CK", "Y11.1"))



#dry weight
library(ggplot2)
index<-read.csv("D://R/出站报告/功能菌盆栽实验/不灭菌/natural soil.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "Y11.1"))

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
  scale_fill_manual(values = c("#C7C7C7FF", "#98DF8AFF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,30))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p
ggsave(paste("D://R/出站报告/功能菌盆栽实验/不灭菌/Dry weight",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")



#soil nematodes
library(ggplot2)
index<-read.csv("D://R/出站报告/功能菌盆栽实验/不灭菌/natural soil.csv")#,row.names=1
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
  
  labs(y="Soil nematode abundance per 100 g dry soil", x="")+
  scale_fill_manual(values = c("#C7C7C7FF", "#98DF8AFF"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,6000))+
  
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p

ggsave(paste("D://R/出站报告/功能菌盆栽实验/不灭菌/soil_nematode",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")

#gall number
library(ggplot2)
index<-read.csv("D://R/出站报告/功能菌盆栽实验/不灭菌/natural soil.csv")#,row.names=1
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
ggsave(paste("D://R/出站报告/功能菌盆栽实验/不灭菌/Gall numbers",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")


#RKN abundance
library(ggplot2)
index<-read.csv("D://R/出站报告/功能菌盆栽实验/不灭菌/natural soil.csv")#,row.names=1
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

ggsave(paste("D://R/出站报告/功能菌盆栽实验/不灭菌/RKN_abundance",".pdf",sep=""),
       device=cairo_pdf,width=40,height=55,dpi = 300,units = "mm")


#Fig 6d chemotaxis index

mytheme<- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#000000"),
                            
                            strip.text = element_text(size=8,hjust = 0.5),
                            plot.title = element_text(size=8,hjust = 0.5),
                            axis.text =element_text(size=8,color = "black"),
                            axis.title =element_text(size=8,color = "black"),
                            legend.text = element_text(size=8,color = "black"),
                            legend.title = element_text(size=8,color = "black"),
                            legend.background = element_blank(),
                            axis.line = element_line(color = "black",size=0.4))#移除整体的边???

library(ggplot2)
CI<-read_excel("D://R/banana data/Fig.6/趋化性/RKN-CELE.xlsx")
CI$Treatment <-factor(CI$Treatment,levels=c("CK1", "R10", "R100", "R500","R1000","CK2", "10", "100", "500","1000"))
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
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF","#9EDAE5FF","#9EDAE5FF","#9EDAE5FF",
                               "#C7C7C7FF", "#9EDAE5FF","#9EDAE5FF","#9EDAE5FF","#9EDAE5FF"))+
  
  scale_y_continuous(limits = c(-0.8,0.8))+
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = "none")
p

ggsave(paste("D://R/banana data/Fig.6/趋化性/趋化性",".pdf",sep=""),
       device=cairo_pdf,width=80,height=70,dpi = 300,units = "mm")


#Fig 6e morality
index <- read_excel("D://R/banana data/Fig.6/致死率/RKN-CELE/Death rate.xlsx")
index$Treatment <-factor(index$Treatment,levels=c("RCK", "R50", "R100", "R250","R500","R1000","CK", "50", "100", "250","500","1000"))
leveneTest(death_rate~ Treatment, data =index)#p>0.05，则满足方差齐性
shapiro.test(index$death_rate)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_index<-aov(data=index,death_rate~Treatment)
summary(aov_model_index)
LSD.test(aov_model_index,"Treatment",p.adj = "none",console=T)
LSD.test(aov_model_index,"Treatment",p.adj = "none",console=T,group = F)

p <- ggplot(data=index,mapping=aes(x=Treatment,y=death_rate,fill=Time))+
  geom_bar(size = 0.25,color="black",fun="mean",position=position_dodge(0.8), stat="summary",width = 0.65)+  
  stat_summary(fun.data=function(...) mean_sdl(..., mult=1),position=position_dodge(0.8), 
               geom='errorbar',width=0.2,size=0.2)+
  geom_jitter (data = index, aes(x=Treatment,y=death_rate,fill=Time),shape=21,
               position=position_jitterdodge(0.3), size=0.6, stroke = 0.2)+
  labs(y="Morality (%)", x="")+
  scale_y_continuous(expand = c(0, 0),limits = c(-0.5,25))+
  theme(axis.text=element_text(colour='black',size=8))+
  mytheme+
  theme(legend.position = c(0.2,0.9), legend.title = element_blank(),
        text = element_text(colour='black'), axis.text=element_text(colour='black'))

p
ggsave(paste("D://R/banana data/Fig.6/致死率/RKN-CELE/Death rate-1",".pdf",sep=""),
       device=cairo_pdf,width=100,height=70,dpi = 300,units = "mm")



