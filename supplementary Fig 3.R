rm(list=ls())
library(tidyverse)
library(RColorBrewer)
library(ggpmisc)
library(ggpubr)
module <- read.delim("ext_module.txt")
module <- filter(module,module %in% c(1,2))
otutab <- read.delim("otutab_rare.txt")
taxonomy <- read.delim("taxonomy.txt")
meloidogyne <- read.delim("meloidogyne.txt")
metadata <- read.delim("metadata.txt")
colnames(module)[1]="OTUID"
#依据OTUID列合并上述三个数据框
alltab <- filter(module, module %in% c(1,2)) %>% 
  left_join(.,otutab,by="OTUID") %>% 
  left_join(.,taxonomy,by="OTUID")
#alltab的module列前加"Module_"前缀
alltab$module <- paste("Module_",alltab$module,sep="")
alltab <- alltab %>% 
  select(-c("OTUID","Kingdom","Class","Order","Family","Genus","Species"))
alltab$mean <- rowMeans(alltab[,2:33])
alltab <- alltab[,c(1, 34,35)]
alltab <- alltab %>% group_by(module) %>% 
  reframe(mean2=mean/sum(mean)*100, Phylum = Phylum)
alltab2 <- alltab %>% group_by(module,Phylum) %>% 
  summarise(mean2 = sum(mean2))
write.csv(alltab2,"alltab2.csv",row.names = F)
p <- ggplot(alltab2, aes(x = module, y = mean2, fill = Phylum))+
  geom_col() + theme_bw() + 
  labs(y = "Prop. ASV numbers", x="")+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(colour='black',size=20), 
        axis.text=element_text(colour='black',size=17))
p
otutab[,-1] <- otutab[,-1]/11176
lmall <- t(otutab)
colnames(lmall) <- lmall[1,]
lmall <- lmall[-1,]
lmall <- as.data.frame(lmall)
#将每一列变为数字型
lmall <- lmall %>% mutate_all(as.numeric)
lmall$sample_name <- rownames(lmall)
lmall <- merge(lmall,meloidogyne,by="sample_name")
asvm1 <- filter(module, module == 1)$OTUID
lmm1 <- lmall %>% select(sample_name,meloidogyne,nematode,gall,disease,all_of(asvm1))
lmm1$ASVSUM <- rowSums(lmm1[,6:80])
lmm1$Years <- metadata$Years
lmm1$Modules <- "Module_1"
lmm1 <- lmm1 %>% select(ASVSUM,meloidogyne,nematode,gall,disease,Years,Modules)
asvm2 <- filter(module, module ==2)$OTUID
lmm2 <- lmall %>% select(sample_name,meloidogyne,nematode,gall,disease,all_of(asvm2))
lmm2$ASVSUM <- rowSums(lmm2[,6:104])
lmm2$Years <- metadata$Years
lmm2$Modules <- "Module_2"
lmm2 <- lmm2 %>% select(ASVSUM,meloidogyne,nematode,gall, disease,Years,Modules)

lmm <- rbind(lmm1,lmm2)

#画图
formula <- y ~ x
lmm$meloi <- log10(lmm$meloidogyne)
lmm$asv <- log10(lmm$ASVSUM)


p4 <- ggarrange(p2,p3,nrow = 1,common.legend = TRUE,legend="right")
p4
p5 <- ggplot(lmm, aes(x = asv, y = gall)) + geom_point(aes(colour = Modules,shape = Years), size = 3) +
  geom_smooth(method = "lm", aes(colour = Modules, fill = Modules)) +
  scale_colour_manual(values = c("#fba801","#87ceeb"))+
  scale_fill_manual(values = c("#fba801","#87ceeb"))+
  stat_cor(aes(color=Modules),method = "pearson",label.x = -1.7, size = 5) +
  theme_test() +
  theme(text = element_text(colour='black',size=17), axis.text=element_text(colour='black',size=13),
        legend.text = element_text(size = 15)) +
  labs(y = "Gall index", x="Relative abundance of module \n (log-transformed)")
p5
p6 <- ggplot(lmm, aes(x = asv, y = disease)) + geom_point(aes(colour = Modules,shape = Years), size = 3) +
  geom_smooth(method = "lm", aes(colour = Modules, fill = Modules)) +
  scale_colour_manual(values = c("#fba801","#87ceeb"))+
  scale_fill_manual(values = c("#fba801","#87ceeb"))+
  stat_cor(aes(color=Modules),method = "pearson",label.x = -1.7, size = 5) +
  theme_test() +
  theme(text = element_text(colour='black',size=17), axis.text=element_text(colour='black',size=13),
        legend.text = element_text(size = 15)) +
  labs(y = "Disease incidence (%)", x="Relative abundance of module \n (log-transformed)")
p6
p7 <-  ggarrange(p5,p6,nrow = 1,common.legend = TRUE,labels = c("a","b"),legend="right")
p7
