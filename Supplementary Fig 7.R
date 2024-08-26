
library(ggplot2)
content<-read.csv("D://R/banana data/功能肥-功能菌/铁载体浓度/2023-11-7/bacillibactin.csv")

p<-ggplot()+geom_bar(data=content,mapping=aes(x=Treatment,y=mean,fill=Treatment),size = 0.75, 
                  position=position_dodge(0.9), stat="identity",width = 0.5)+  
  
  scale_fill_manual(values = c("#C7C7C7FF"))+ 
  geom_errorbar(data=content,mapping=aes(x = Treatment, ymin = mean-sd,ymax = mean+sd), width = 0.15, 
                color = 'black',size=0.5)+
  labs(y="Siderophore concentration in supernate \n (μM equivalents of DFOB) ", x="")+
  
  scale_y_continuous(expand = c(0, 0),limits = c(0,60))+
  theme_test(base_line_size = 0.75,base_rect_size =0.75)+
  theme(legend.position = "none")+
  theme(axis.text=element_text(colour='black',size=12))
p

ggsave(paste("D://R/banana data/功能肥-功能菌/铁载体浓度/2023-11-7/bacillibactin-3",".pdf",sep=""),
       device=cairo_pdf,width=60,height=80,dpi = 300,units = "mm")
