rm(list=ls())
library(tidyverse)
library(psych)
# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
#连作前后期差异OTU表
ind <- read.delim("richY.txt")
ind$rich <- ifelse(ind$rich == "#1f78b4", "Y08Y11", ind$rich)
ind$rich <- ifelse(ind$rich == "#228b22", "Y02Y05", ind$rich)
#线虫指标表
index <- read.delim("meloidogyne.txt")
index$meloidogyne <- log10(index$meloidogyne)
#实验设计表
metadata <- read.delim("metadata.txt")
#细菌丰度表
otutab <- read.delim("otutab_rare.txt")
otutab[,-1] <- otutab[,-1]/colSums(otutab[,-1])
#细菌分类表
tax <- read.delim("taxonomyB.txt")
#细菌网络模块表
module <- read.delim("ext_module.txt")
colnames(module)[1] <- "OTUID"
#提取ind中color为Y02Y05的数据,提取module中module为2的列并与indY2Y5按OTUID列合并，取两个数据框的交集
indY2Y5 <- ind %>% filter(rich == "Y02Y05") %>% inner_join(module, by = "OTUID") %>% filter(module == 1)
#同理提取差异Y8Y11
indY8Y11 <- ind %>% filter(rich == "Y08Y11") %>% inner_join(module, by = "OTUID") %>% filter(module == 2)
#对otutab除第一列外的所有列按行求和并删除多余列
otutab_sum <- otutab %>% mutate(sum = rowSums(.[,-1])) %>% select(OTUID, sum)
#将indY2Y5与otutab_sum按OTUID列合并，取交集，并按sum列降序排列
indY2Y5_sum <- indY2Y5 %>% inner_join(otutab_sum, by = "OTUID") 
indY2Y5 <- indY2Y5_sum %>% arrange(desc(sum))
indY2Y5 <- indY2Y5[1:20,]
indY8Y11_sum <- indY8Y11 %>% inner_join(otutab_sum, by = "OTUID")
indY8Y11 <- indY8Y11_sum %>% arrange(desc(sum))
indY8Y11 <- indY8Y11[1:20,]
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
colnames(index)[1] <- "SampleID"
#将otutab_ind与index按SampleID列合并
index <- otutab_ind %>% inner_join(index, by = "SampleID")
#将index与metadata按SampleID列合并
index <- index %>% inner_join(metadata, by = "SampleID")
#以下热图绘制代码改自王男麒师兄github上nc代码
#### 7. Fig. 2b-Heatmap_genus_biomarker_MPvsIP####
### 7.1 Import and process data ###
index_mean <- index %>% group_by(Years,Group) %>%
  summarise_at(vars(ASV_6213:ASV_774),funs(mean))
index_relativeabundance <- index_mean[,-c(1:2)]
index_relativeabundance_norm <- t(scale(index_relativeabundance,center = T))
colnames(index_relativeabundance_norm) <- index_mean$Years
max(index_relativeabundance_norm)
min(index_relativeabundance_norm)
col_fun = circlize::colorRamp2(c(-1.5, 0, 1.5), c("#2166ac", "white", "#b2182b"))
col_fun(seq(-3, 3))
column_split_order<-c(rep("Y02",1),rep("Y05",1),rep("Y08",1),rep("Y11",1))
column_split_order<-factor(column_split_order, levels = c("Y02","Y05","Y08","Y11"))
index_Y08Y11 <- filter(index, Group == "Y08Y11")
index_Y08Y11_value <- index_Y08Y11[,1:40]
index_Y08Y11_value <- t(index_Y08Y11_value)
ha = rowAnnotation(RA_mean = anno_boxplot(index_Y08Y11_value,height = unit(4, "cm"),
                                          box_width = 0.5))
ha
## 7.2.4 plots ##
Heatmap(index_relativeabundance_norm,col = col_fun,
        cluster_columns = F,column_split = column_split_order,
        name="Relative abundance",row_km = 4,row_gap = unit(c(2,4,2), "mm"),
        column_title_gp = gpar(fill = c("#E31A1C", "#1F78B4", "#A6CEE3")),
        right_annotation = ha, row_dend_reorder = TRUE,
        clustering_distance_rows = "spearman")
inf_rate <- corr.test(index$disease, index[,1:40], method="spearman",adjust="BH",minlength=5)
inf_rate$r
inf_rate$p
inf_rate$p.adj
hehe1 <- as.data.frame(rbind(inf_rate$r,inf_rate$p))
hehe1 <- hehe1[,c("ASV_8459","ASV_6972","ASV_8453","ASV_8525","ASV_4693",
                  "ASV_4365","ASV_7026","ASV_8041","ASV_7704","ASV_774",
                  "ASV_3638","ASV_7123","ASV_6008","ASV_5963","ASV_1467",
                  "ASV_4833","ASV_5065","ASV_6738","ASV_447","ASV_6459",
                  "ASV_7438","ASV_7216","ASV_5729","ASV_6522","ASV_7537",
                  "ASV_5223","ASV_6061","ASV_7898","ASV_3999","ASV_44",
                  "ASV_5244","ASV_4472","ASV_5425","ASV_7676","ASV_5305",
                  "ASV_4585","ASV_3920","ASV_6213","ASV_4459","ASV_8058")]
meloidogyne <- corr.test(index$meloidogyne, index[,1:40], method="spearman",adjust="BH",minlength=5)
meloidogyne$r
meloidogyne$p
meloidogyne$p.adj
hehe2 <- as.data.frame(rbind(meloidogyne$r,meloidogyne$p))
hehe2 <- hehe2[,c("ASV_8459","ASV_6972","ASV_8453","ASV_8525","ASV_4693",
                  "ASV_4365","ASV_7026","ASV_8041","ASV_7704","ASV_774",
                  "ASV_3638","ASV_7123","ASV_6008","ASV_5963","ASV_1467",
                  "ASV_4833","ASV_5065","ASV_6738","ASV_447","ASV_6459",
                  "ASV_7438","ASV_7216","ASV_5729","ASV_6522","ASV_7537",
                  "ASV_5223","ASV_6061","ASV_7898","ASV_3999","ASV_44",
                  "ASV_5244","ASV_4472","ASV_5425","ASV_7676","ASV_5305",
                  "ASV_4585","ASV_3920","ASV_6213","ASV_4459","ASV_8058")]
gall <- corr.test(index$gall, index[,1:40], method="spearman",adjust="BH",minlength=5)
gall$r
gall$p
abundance <- corr.test(index$nematode, index[,1:40], method="spearman",adjust="BH",minlength=5)
abundance$r
abundance$p
hehe3 <- as.data.frame(rbind(abundance$r,abundance$p))
hehe3 <- hehe3[,c("ASV_8459","ASV_6972","ASV_8453","ASV_8525","ASV_4693",
                  "ASV_4365","ASV_7026","ASV_8041","ASV_7704","ASV_774",
                  "ASV_3638","ASV_7123","ASV_6008","ASV_5963","ASV_1467",
                  "ASV_4833","ASV_5065","ASV_6738","ASV_447","ASV_6459",
                  "ASV_7438","ASV_7216","ASV_5729","ASV_6522","ASV_7537",
                  "ASV_5223","ASV_6061","ASV_7898","ASV_3999","ASV_44",
                  "ASV_5244","ASV_4472","ASV_5425","ASV_7676","ASV_5305",
                  "ASV_4585","ASV_3920","ASV_6213","ASV_4459","ASV_8058")]
# Free_living_nematode <- corr.test(index$Free_living_nematode, index[,1:40], method="spearman",adjust="BH",minlength=5)
# Free_living_nematode$r
# Free_living_nematode$p
# Herbivorous_nematodes <- corr.test(index$Herbivorous_nematodes, index[,1:40], method="spearman",adjust="BH",minlength=5)
# Herbivorous_nematodes$r
# Herbivorous_nematodes$p
cor <- as.data.frame(t(rbind(inf_rate$r, inf_rate$p, meloidogyne$r,meloidogyne$p,abundance$r,abundance$p)))
colnames(cor) <- c("inf_rate_r","inf_rate_p","meloidogyne_r","meloidogyne_p","abundance_r","abundance_p")
cor$ASV <- rownames(cor)
cor_r <- as.data.frame(t(rbind(inf_rate$r, meloidogyne$r,abundance$r)))
colnames(cor_r) <- c("inf_rate_r","meloidogyne_r","abundance_r")
max(cor_r)
min(cor_r)
col_cor = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#4d9221", "white", "#c51b7d"))
row_order <- c("ASV_8459","ASV_6972","ASV_8453","ASV_8525","ASV_4693",
               "ASV_4365","ASV_7026","ASV_8041","ASV_7704","ASV_774",
               "ASV_3638","ASV_7123","ASV_6008","ASV_5963","ASV_1467",
               "ASV_4833","ASV_5065","ASV_6738","ASV_447","ASV_6459",
               "ASV_7438","ASV_7216","ASV_5729","ASV_6522","ASV_7537",
               "ASV_5223","ASV_6061","ASV_7898","ASV_3999","ASV_44",
               "ASV_5244","ASV_4472","ASV_5425","ASV_7676","ASV_5305",
               "ASV_4585","ASV_3920","ASV_6213","ASV_4459","ASV_8058")
length(row_order)
row_order
Heatmap(cor_r,col = col_cor,
        cluster_columns = F,cluster_rows = F,
        row_order = row_order
)
