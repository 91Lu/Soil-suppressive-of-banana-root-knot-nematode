第一种算法
library(agricolae)
library(dplyr)
library(car)
library(ggplot2)

index<-read.csv("D://R/banana data/功能肥-功能菌/趋化实验/2023-11-20/缺铁培养基.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "10", "100","500", "1000"))

m="chemotaxis1"
model = aov(index[[m]] ~ group, data=index)

Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
Tukey_HSD_table = as.data.frame(Tukey_HSD$group)

out = LSD.test(model,"group", p.adj="none")
stat = out$groups

index$stat=stat[as.character(index$group),]$groups#
df1<-index[,c(2,3)]
data1<-df1%>%group_by(group)%>%summarise_at("chemotaxis1",funs(mean,sd))
data1
df2<-as.data.frame(data1)
df3<-cbind(df2,stat$groups)
colnames(df3)<-c("group","mean","sd","stat")

p <- ggplot(data=df2,mapping=aes(x=group,y=mean,fill=group))+
  geom_bar(size = 0.5, position=position_dodge(0.9), stat="identity",width = 0.7)+  
  geom_jitter(data=df1,mapping=aes(x=group,y=chemotaxis1),color="black",
              size = 0.71,height = 0.02,width = 0.1)+ 
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF","#9EDAE5FF","#9EDAE5FF","#9EDAE5FF"))+ 
  scale_color_manual(values = c("#C7C7C7FF", "#9EDAE5FF","#9EDAE5FF","#9EDAE5FF","#9EDAE5FF"))+
  geom_errorbar(data=df2,aes(x = group,ymin = mean-sd, ymax = mean+sd), width = 0.25,color="black",size=0.3)+
  labs(y="Chemotaxis index", x="")+
  scale_y_continuous(expand = c(0, 0),limits = c(-1,0.2))+
  theme_test(base_line_size = 0.75,base_rect_size =0.75)+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(colour='black',size=10),axis.text.x = element_text(vjust = 0.85,hjust = 0.75))+
  geom_text(data=df3,aes(x = group,y = mean-sd-0.05,label = stat), position = position_dodge(0.75),size=4)

p

ggsave(paste("D://R/banana data/功能肥-功能菌/趋化实验/2023-11-20/缺铁培养基",".pdf",sep=""),
       device=cairo_pdf,width=60,height=70,dpi = 300,units = "mm")


第二种算法
library(agricolae)
library(dplyr)
library(car)
library(ggplot2)

index<-read.csv("D://R/banana data/功能肥-功能菌/趋化实验/2023-11-20/缺铁培养基.csv")#,row.names=1
index$group <-factor(index$group,levels=c("CK", "10", "100","500", "1000"))

m="chemotaxis2"
model = aov(index[[m]] ~ group, data=index)

Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
Tukey_HSD_table = as.data.frame(Tukey_HSD$group)

out = LSD.test(model,"group", p.adj="none")
stat = out$groups

index$stat=stat[as.character(index$group),]$groups#
df1<-index[,c(2,4)]
data1<-df1%>%group_by(group)%>%summarise_at("chemotaxis2",funs(mean,sd))
data1
df2<-as.data.frame(data1)
df3<-cbind(df2,stat$groups)
colnames(df3)<-c("group","mean","sd","stat")

p <- ggplot(data=df2,mapping=aes(x=group,y=mean,fill=group))+
  geom_bar(size = 0.5, position=position_dodge(0.9), stat="identity",width = 0.7)+  
  geom_jitter(data=df1,mapping=aes(x=group,y=chemotaxis2),color="black",
              size = 0.71,height = 0.02,width = 0.1)+ 
  scale_fill_manual(values = c("#C7C7C7FF", "#9EDAE5FF","#9EDAE5FF","#9EDAE5FF","#9EDAE5FF"))+ 
  scale_color_manual(values = c("#C7C7C7FF", "#9EDAE5FF","#9EDAE5FF","#9EDAE5FF","#9EDAE5FF"))+
  geom_errorbar(data=df2,aes(x = group,ymin = mean-sd, ymax = mean+sd), width = 0.25,color="black",size=0.3)+
  labs(y="Chemotaxis index", x="")+
  scale_y_continuous(expand = c(0, 0),limits = c(-0.5,0.2))+
  theme_test(base_line_size = 0.75,base_rect_size =0.75)+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(colour='black',size=10),axis.text.x = element_text(vjust = 0.85,hjust = 0.75))+
  geom_text(data=df3,aes(x = group,y = mean-sd-0.02,label = stat), position = position_dodge(0.75),size=4)

p

ggsave(paste("D://R/banana data/功能肥-功能菌/趋化实验/2023-11-20/缺铁培养基-1",".pdf",sep=""),
       device=cairo_pdf,width=60,height=70,dpi = 300,units = "mm")

