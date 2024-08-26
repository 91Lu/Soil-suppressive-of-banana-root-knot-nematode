#Fig.S3a
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
