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
#以下热图绘制代码改自王男麒师兄github上nc代码
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
rm(list = ls())
library(tidyverse)
library(randomForest)
library(rfPermute)
library(A3)

meloi <- read.delim("meloidogyne.txt")
meloi$meloidogyne <- log10(meloi$meloidogyne)
otu <- read.delim("otutab_rare.txt")
ind <- read.delim("richY.txt")
# rich列的#1f78b4替换为Y08Y11,#228B22替换为Y02Y05
ind$rich <- ifelse(ind$rich == "#1f78b4", "Y08Y11", ind$rich)
ind$rich <- ifelse(ind$rich == "#228b22", "Y02Y05", ind$rich)

mod <- read.delim("ext_module.txt")
mod <- mod %>% filter(module %in% c(1,2))
ind <- ind %>% filter(rich %in% c("Y02Y05","Y08Y11"))
id <- intersect(mod$ASV, ind$OTUID)
mod <- mod %>% filter(ASV %in% id)
colnames(mod)[1] <- "OTUID"
ind <- ind %>% filter(OTUID %in% id)
otu <- otu  %>% filter(OTUID %in% id)
# join按OTUID合并otu和mod
otu_mod <- left_join(otu, mod, by = "OTUID")
otu_mod2 <- otu_mod %>% filter(module==2)
otu_mod1 <- otu_mod %>% filter(module==1)
# 删除otu_mod2中的module列
otu_mod2 <- otu_mod2[,-ncol(otu_mod2)]
rownames(otu_mod2) <- otu_mod2[,1]
otu_mod2 <- otu_mod2[,-1]
# 求行和
otu_mod2$sum <- rowSums(otu_mod2)
# 按sum列降序排列
otu_mod2 <- otu_mod2[order(-otu_mod2$sum),]
otu_mod2 <- otu_mod2[1:20,]
# 删除otu_mod1中的module列
otu_mod1 <- otu_mod1[,-ncol(otu_mod1)]
rownames(otu_mod1) <- otu_mod1[,1]
otu_mod1 <- otu_mod1[,-1]
otu_mod1$sum <- rowSums(otu_mod1)
otu_mod1 <- otu_mod1[order(-otu_mod1$sum),]
otu_mod1 <- otu_mod1[1:20,]
# 合并otu_mod1和otu_mod2
otu_mod <- rbind(otu_mod1,otu_mod2)
otu_mod <- otu_mod[,-ncol(otu_mod)]
totu <- as.data.frame(t(otu_mod))
totu$meloi <- meloi$meloidogyne
set.seed(1234)
otu_rfP <- rfPermute(meloi~., data = totu, importance = TRUE, ntree = 500, nrep = 1000)
otu_rfP
set.seed(1234)
otu_rfd <- randomForest(meloi~., data = totu, importance = TRUE, ntree = 500, nrep = 1000)
otu_rfd
#提取预测变量（细菌 OTU）的重要性得分（标准化后的得分）
importance_otu.scale <- data.frame(importance(otu_rfP, scale = TRUE), check.names = FALSE)
importance_otu.scale
totu <- totu %>% select(rownames(importance_otu.scale)[1:15],"meloi")
set.seed(1234)
otu_rfP <- rfPermute(meloi~., data = totu, importance = TRUE, ntree = 500, nrep = 1000)
otu_rfP
set.seed(1234)
otu_rfd <- randomForest(meloi~., data = totu, importance = TRUE, ntree = 500, nrep = 1000)
otu_rfd
#提取预测变量（细菌 OTU）的重要性得分（标准化后的得分）
importance_otu.scale <- data.frame(importance(otu_rfP, scale = TRUE), check.names = FALSE)
importance_otu.scale
#提取预测变量（细菌 OTU）的重要性得分的显著性（以标准化后的得分为例）
# summary(otu_rfP)
importance_otu.scale.pval <- (otu_rfP$pval)[ , , 2]
importance_otu.scale.pval
#对预测变量（细菌 OTU）按重要性得分排个序，例如根据“%IncMSE”
importance_otu.scale <- importance_otu.scale[order(importance_otu.scale$'%IncMSE', decreasing = TRUE), ]
importance_otu.scale
write.csv(importance_otu.scale, "importance_otu_scale.csv")

#简单地作图展示预测变量（细菌 OTU）的 %IncMSE 值
importance_otu.scale$OTU_name <- rownames(importance_otu.scale)
importance_otu.scale$OTU_name <- factor(importance_otu.scale$OTU_name, levels = importance_otu.scale$OTU_name)
#标记预测变量（细菌 OTU）的显著性信息
#默认以 p<0.05 为 *，p<0.01 为 **，p<0.001 为 ***
for (OTU in rownames(importance_otu.scale)) {
  importance_otu.scale[OTU,'%IncMSE.pval'] <- importance_otu.scale.pval[OTU,'%IncMSE']
  if (importance_otu.scale[OTU,'%IncMSE.pval'] >= 0.05) importance_otu.scale[OTU,'%IncMSE.sig'] <- ''
  else if (importance_otu.scale[OTU,'%IncMSE.pval'] >= 0.01 & importance_otu.scale[OTU,'%IncMSE.pval'] < 0.05) importance_otu.scale[OTU,'%IncMSE.sig'] <- '*'
  else if (importance_otu.scale[OTU,'%IncMSE.pval'] >= 0.001 & importance_otu.scale[OTU,'%IncMSE.pval'] < 0.01) importance_otu.scale[OTU,'%IncMSE.sig'] <- '**'
  else if (importance_otu.scale[OTU,'%IncMSE.pval'] < 0.001) importance_otu.scale[OTU,'%IncMSE.sig'] <- '***'
}
colnames(mod)[1] <- "OTU_name"
# join合并两个表格
importance_otu.scale <- left_join(importance_otu.scale, mod, by = "OTU_name")
# module列每个值前加Module
importance_otu.scale$module <- paste0("Module", importance_otu.scale$module)
importance_otu.scale <- importance_otu.scale[1:10,]
p <- ggplot(importance_otu.scale, aes(y = reorder(OTU_name,`%IncMSE`), x = `%IncMSE`, fill = module)) +
  geom_col(position = 'dodge', width = 0.7, color = "black") +
  labs(title = NULL, y = NULL, x = 'Increase in MSE (%)') +
  theme_test()+
  # theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  scale_x_continuous(expand = c(0, 0), limit = c(0, 10)) + 
  # geom_text(aes(y = reorder(OTU_name,`%IncMSE`), x = `%IncMSE`, label = `%IncMSE.sig`),nudge_x = 0.1)+
  theme(legend.position = c(0.9,0.1),legend.title = element_blank()) + scale_fill_manual(values = c("#fba801","#87ceeb"))
p
set.seed(1234)
otu_forest.pval <- a3(meloi~., data = totu, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
write.csv(print(otu_forest.pval),"oru_forest.pval.csv")
p <- p +
  # annotate('text', label = 'Plant Age', x = 9, y = 15, size = 4) +
  annotate('text', label = sprintf('italic(R^2) == %.2f', 26.1), x = 7.5, y = 5, size = 4, parse = TRUE) +
  annotate('text', label = sprintf('italic(p) == %.2f', 0.009), x = 7.5, y = 4.5, size = 4, parse = TRUE)

p
save.image()


#Fig.5c 
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

