library(tidyverse)
library(agricolae)
library(dplyr)
library(car)
library(ggplot2)
library(MASS)
library(multcompView)
index <- read.csv("D://R/banana data/功能肥-功能菌/致死率/LB培养基2023-10-14/death rate1.csv")
mean <- index %>% group_by(treatment,Time) %>% 
  summarise(mean = mean(death_rate), se = sd(death_rate)/sqrt(n()))
i24h <- filter(index,Time == "24h")
i24h$treatment <- factor(i24h$treatment, levels = c("CK1","1x","0.5x","0.1x","CK2"))
i48h <- filter(index,Time == "48h")
i48h$treatment <- factor(i48h$treatment, levels = c("CK1","1x","0.5x","0.1x","CK2"))

model<- aov(death_rate ~ treatment, data=i24h)
Tukey_HSD <- TukeyHSD(model, conf.level = 0.95)
Tukey_HSD_table <- as.data.frame(Tukey_HSD$treatment)
significant <- Tukey_HSD_table[,4]
names(significant) <- rownames(Tukey_HSD_table)

dif <- multcompLetters(significant)
dif2 <- as.data.frame(dif$monospacedLetters)
colnames(dif2) <- "label"
dif2$label <- toupper(gsub(" ", "", dif2$label))
dif2$treatment <- row.names(dif2)
dif2$Time <- "24h"

model<- aov(death_rate ~ treatment, data=i48h)
Tukey_HSD <- TukeyHSD(model, conf.level = 0.95)
Tukey_HSD_table <- as.data.frame(Tukey_HSD$treatment)
significant <- Tukey_HSD_table[,4]
names(significant) <- rownames(Tukey_HSD_table)

dif <- multcompLetters(significant)
dif3 <- as.data.frame(dif$monospacedLetters)
colnames(dif3) <- "label"
dif3$label <- gsub(" ", "", dif3$label)
dif3$treatment <- row.names(dif3)
dif3$Time <- "48h"

dif <- rbind(dif2, dif3)
mean <- merge(mean,dif,by = c("treatment", "Time"))
mean$treatment <- factor(mean$treatment, levels = c("CK1","CK2","0.1x","0.5x","1x"))
p1 <- ggplot(mean, aes(x = treatment, y = mean, fill = Time)) +
  geom_bar(stat = "identity",position = position_dodge(0.8),width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width=.3, position = position_dodge(.8)) +
  labs(y = "Morality (%)", x="") + 
  
  scale_y_continuous(expand = c(0, 0),limits = c(0,80)) +
  theme_test(base_rect_size = 1)+
  theme_classic()+
  theme(legend.position = c(0.08,0.92), legend.title = element_blank(),
        text = element_text(colour='black',size=16), axis.text=element_text(colour='black',size=13)) +
  geom_text(data=mean,aes(x = treatment,y = mean + se+6,label = label), position = position_dodge(0.8),size=4)+
        geom_jitter (data = index, aes(x=treatment,y=death_rate),
                     position=position_jitterdodge(0.3), size=1.0, alpha=0.5,color="black")
  
p1  
ggsave(paste("D://R/banana data/功能肥-功能菌/致死率/LB培养基2023-10-14/death rate-7",".pdf",sep=""),
       device=cairo_pdf,width=100,height=80,dpi = 300,units = "mm")
