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
library(showtext)#字体设置, font setting
library(car)#方差齐性检验，homogeneity test of variance, levene test
library(extrafont)#使用系统字体，Using the system fonts
library(sysfonts)#加载系统字体，loading the system fonts
library(Cairo)#抗锯齿,anti-aliasing
library(stringr)#字符串处理.string manipulation
library(graphics)#坐标轴表达式，expression for axis
library(vegan)
library(data.table)
library(rcompanion)
library(forecast)#box-cox数据变换
library(PMCMRplus)
#### 2.

library(ggthemes)
mytheme<- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#000000"),
                                
                                strip.text = element_text(size=8,hjust = 0.5),
                                plot.title = element_text(size=8,hjust = 0.5),
                                axis.text =element_text(size=8,color = "black"),
                                axis.title =element_text(size=8,color = "black"),
                                legend.text = element_text(size=8,color = "black"),
                                legend.title = element_text(size=8,color = "black"),
                                legend.background = element_blank(),
                                axis.line = element_line(color = "black",size=0.4))#移除整体的边???



#Fig. 4a Shannon
Shannon_Index<-read_excel("D://R/banana data/Fig.4/2024-3-11/Fig.4.xlsx", sheet="Bacterial_diversity")

Shannon_Index$Treatment <-factor(Shannon_Index$Treatment,levels=c("Y2", "Y5", "Y8", "Y11"))
##Fig. 4_Shannon_Index statistical analysis
leveneTest(Shannon ~ Treatment, data = Shannon_Index)#p>0.05，则满足方差齐性
shapiro.test(Shannon_Index$Shannon)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_Shannon_Index<-aov(data=Shannon_Index,Shannon~Treatment)
summary(aov_model_Shannon_Index)
LSD.test(aov_model_Shannon_Index,"Treatment",console=T)
LSD.test(aov_model_Shannon_Index,"Treatment",console=T,group = F)


#Fig. 4a Shannon_Index boxplot
p <- ggplot(data=Shannon_Index,mapping=aes(x=Treatment,y=Shannon))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,size=0.2)+  
  geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.6,stroke = 0.2)+
  
  labs(y="Shannon Index", x="",parse =T)+
  scale_fill_manual(values = c("#9467BDFF","#EEA236FF", "#5CB85CFF","#46B8DAFF"))+
  scale_y_continuous(limits = c(7.0,8.5))+
  mytheme+
  guides(color="none",fill="none")
p 

ggsave(paste("D://R/banana data/Fig.4/2024-3-11/Shannong-boxplot",".pdf",sep=""),
       device=cairo_pdf,width=35,height=40,dpi = 300,units = "mm")


#Fig. 4a Chao1
Chao1_Index<-read_excel("D://R/banana data/Fig.4/2024-3-11/Fig.4.xlsx", sheet="Bacterial_diversity")

Chao1_Index$Treatment <-factor(Chao1_Index$Treatment,levels=c("Y2", "Y5", "Y8", "Y11"))
##Fig. 4_Chao1_Index statistical analysis
leveneTest(Chao1 ~ Treatment, data = Chao1_Index)#p>0.05，则满足方差齐性
shapiro.test(Chao1_Index$Chao1)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_Chao1_Index<-aov(data=Chao1_Index,Chao1~Treatment)
summary(aov_model_Chao1_Index)
LSD.test(aov_model_Chao1_Index,"Treatment",console=T)
LSD.test(aov_model_Chao1_Index,"Treatment",console=T,group = F)


#Fig. 4a Chao1_Index boxplot
p <- ggplot(data=Chao1_Index,mapping=aes(x=Treatment,y=Chao1))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,size=0.2)+  
  geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.6,stroke = 0.2)+
  
  labs(y="Chao1", x="",parse =T)+
  scale_fill_manual(values = c("#9467BDFF","#EEA236FF", "#5CB85CFF","#46B8DAFF"))+
  scale_y_continuous(limits = c(3000,5000))+
  mytheme+
  guides(color="none",fill="none")
p 

ggsave(paste("D://R/banana data/Fig.4/2024-3-11/Chao1-boxplot",".pdf",sep=""),
       device=cairo_pdf,width=40,height=40,dpi = 300,units = "mm")


#Fig. 4b Abundance
Bacterial_Abundance<-read_excel("D://R/banana data/Fig.4/2024-3-11/Fig.4.xlsx", sheet="Bacterial_abundance")

Bacterial_Abundance$Treatment <-factor(Bacterial_Abundance$Treatment,levels=c("Y2", "Y5", "Y8", "Y11"))
##Fig. 4_Shannon_Index statistical analysis
leveneTest(Abundance ~ Treatment, data = Bacterial_Abundance)#p>0.05，则满足方差齐性
shapiro.test(Bacterial_Abundance$Abundance)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
#BoxCox.lambda(Bacterial_Abundance$YL_ActiveFe,method="loglik")#a number indicating the Box-Cox transformation parameter
Bacterial_Abundance<-Bacterial_Abundance%>%mutate(boxcox_Abundance=BoxCox(Bacterial_Abundance$Abundance,BoxCox.lambda(Bacterial_Abundance$Abundance)))
leveneTest(Abundance ~ Treatment, data = Bacterial_Abundance)#p>0.05，则满足方差齐性
shapiro.test(Bacterial_Abundance$Abundance)
aov_model_Bacterial_Abundance<-aov(data=Bacterial_Abundance,Abundance~Treatment)
summary(aov_model_Bacterial_Abundance)
LSD.test(aov_model_Bacterial_Abundance,"Treatment",console=T)
LSD.test(aov_model_Bacterial_Abundance,"Treatment",console=T,group = F)

#Fig. 4b Bacterial_Abundance boxplot
p <- ggplot(data=Bacterial_Abundance,mapping=aes(x=Treatment,y=Abundance))+
  geom_boxplot(aes(fill=Treatment),width = 0.65,outlier.alpha = 0.5,outlier.size=0.5,size=0.2)+  
  geom_point(aes(fill=Treatment),shape=21,position = position_jitterdodge(1),size=0.6,stroke = 0.2)+
  
  labs(y="16S rRNA gene copies \n (109 copies/g soil)", x="",parse =T)+
  scale_fill_manual(values = c("#9467BDFF","#EEA236FF", "#5CB85CFF","#46B8DAFF"))+
  scale_y_continuous(limits = c(0,45))+
  mytheme+
  guides(color="none",fill="none")
p 

ggsave(paste("D://R/banana data/Fig.4/2024-3-11/Bacterial_Abundance-boxplot",".pdf",sep=""),
       device=cairo_pdf,width=40,height=40,dpi = 300,units = "mm")


#Fig. 4c PCoA
library(tidyverse)
library(vegan)
library(ggthemes)
mytheme<- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#000000"),
                            
                            strip.text = element_text(size=8,hjust = 0.5),
                            plot.title = element_text(size=8,hjust = 0.5),
                            axis.text =element_text(size=8,color = "black"),
                            axis.title =element_text(size=8,color = "black"),
                            legend.text = element_text(size=8,color = "black"),
                            legend.title = element_text(size=8,color = "black"),
                            legend.background = element_blank(),
                            axis.line = element_line(color = "black",size=0.4))#移除整体的边???


metadata <- read.table("D://R/banana data/Fig.4/2024-3-11/PCoA/metadata.csv",sep = ",",header=T)
otu <- read.table("D://R/banana data/Fig.4/2024-3-11/PCoA/otu_table.csv",sep = ",",header=T, row.names = 1)
otut <- as.data.frame(t(otu))#转换数据集方向
set.seed(1)
#计算p值和R方
otut.div <- adonis2(otut~ Years, data = metadata, permutations = 999, method="bray")
otut.div
otut_adonis <- paste0("Adonis R2 = ",round(otut.div$R2,2), "; p < ", otut.div$`Pr(>F)`)

otut_dist <- vegdist(otut, method="bray", binary=F)#计算Bray-curtis距离
otut_pcoa <- cmdscale(otut_dist, k=3, eig=T)
otut_pcoa_points <- as.data.frame(otut_pcoa$points)
sum_eig <- sum(otut_pcoa$eig)
eig_percent <- round(otut_pcoa$eig/sum_eig*100,1)#获取维度的结实率
colnames(otut_pcoa_points) <- paste0("PCoA", 1:3)
otut_pcoa_result <- cbind(otut_pcoa_points, metadata)
p <- ggplot(otut_pcoa_result, aes(x=PCoA1, y=PCoA2, color=Years)) +
  scale_color_manual(values = c("#9467BDFF","#EEA236FF", "#5CB85CFF","#46B8DAFF"))+
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
       title = otut_adonis) +
  geom_point(size=1.5) + 
  theme_test() + geom_hline(yintercept = 0, linetype = "dashed")+geom_vline(xintercept = 0,linetype = "dashed")+
  mytheme

p

ggsave(paste("D://R/banana data/Fig.4/2024-3-11/PCoA/PCoA-3",".pdf",sep=""),
       device=cairo_pdf,width=70,height=50,dpi = 300,units = "mm")


#Fig. 4d-f network
rm(list=ls())
library(edgeR)#用于基于edgeR的组间差异OTU识别
library(indicspecies)#用于基于indicator species的组间差异OTU识别
library(igraph)#用于共现性网络的绘制
library(Hmisc)#进行共现性网络构建前的OTU两两相关性计算
library(sciplot)#用于基于模块丰度的绘图
library(reshape2)#用于长宽数据转换
library(ggpmisc)#用于其它观测指标与组间差异模块丰度关系分析
library(tidyverse)

#####输入OTU表格数据#####
otu_16s <- read.table("otutab.txt",row.names=1,sep="\t",header=T, blank.lines.skip=F,check.names=F)
otu_16s <- as.matrix(otu_16s)

#####物种分类表#####
tax_16s <- read.table("taxonomy.txt",row.names=1, sep="\t", header=T,stringsAsFactors=F,quote="")


#####输入metadata数据 #####
design_16s <- read.table("metadata.txt", header=T, row.names=1, na.strings="NA")
design_16s$Group <- factor(design_16s$Group, c("Y02Y05", "Y08Y11"))

#去除一些低发现率的otu
otu_16s <- otu_16s[rowSums(otu_16s != 0) / ncol(otu_16s) > 0.85, ]
tax_16s <- tax_16s[rownames(otu_16s),]

#####寻找基于indicator species的组间显著差异OTU并保存结果#####
edgeR_16s <- DGEList(counts=otu_16s,
                     group=design_16s$Group,
                     genes=tax_16s)
#CPM标准化
otu_norm_16s <- cpm(edgeR_16s, normalized.lib.sizes=T, log=F)
# 准备数据：T时期下组间的indicator species
indic_16s <- as.data.frame(t(otu_norm_16s))
indic_groups_16s <- design_16s$Group
# 设置随机数种子，保证置换检验可重复
set.seed(951001)
# 鉴定各组指示种
indicatorsp_16s <- multipatt(indic_16s,indic_groups_16s,func = "r.g",control=how(nperm=1000))
indic_df_16s <- indicatorsp_16s$sign

#按照阈值p.value < 0.05筛选各组显著的指示OTU
Y02Y05_16s <- as.matrix(indic_df_16s[which(indic_df_16s$s.Y02Y05 == 1 & indic_df_16s$p.value < 0.05),])
Y08Y11_16s <- as.matrix(indic_df_16s[which(indic_df_16s$s.Y08Y11 == 1 & indic_df_16s$p.value < 0.05),])
#合并
r_values_16s <- rbind(Y02Y05_16s,Y08Y11_16s)
# 组名修正，删除多余的"s."
# colnames(t_r_values_16s)[1:3] <- c("CPT","PT","ZT")
colnames(r_values_16s)[1:3] <- gsub("s.","",colnames(r_values_16s)[1:3])

#####寻找基于edgeR的组间显著差异OTU#####
model_mat_16s <- model.matrix(~Group, data=design_16s)
edgeR_16s_Year <- DGEList(counts=otu_16s, group=design_16s$Group, genes=tax_16s)
edgeR_16s_Year <- calcNormFactors(edgeR_16s_Year)
dge_Year_16s <- estimateGLMRobustDisp(edgeR_16s_Year, design=model_mat_16s)
fit_Year_16s <- glmFit(dge_Year_16s, design=model_mat_16s)
# 2，3组分别与1组比较
lrt_Year_16s <- glmLRT(fit_Year_16s, coef=2)
Year_16s <- topTags(lrt_Year_16s, n=Inf, p.value=0.05)
Year_16s <- Year_16s$table

##### (可选)limma分析，取indicator species、edgeR、limma的交集，但是这样可能太严格导致OTU极少
#limma_voom_t_16s <- voom(edgeR_16s_t_tillage, model_matt_16s)
#fit_t_16s <- lmFit(limma_voom_t_16s, model_matt_16s)
#fit_t_16s <- eBayes(fit_t_16s)
#limma_t_16s<-topTable(fit, coef = 2:3)
#indic_edge_16s_t <- Reduce(intersect,list(rownames(t_r_values_16s),rownames(tillage_t_16s),rownames(limma_t_16s)))

######绘制能够体现组件丰度差异OTU/模块的共现性网络图co-occurence network #####
#取基于indicator species分析和edgeR分析得到的显著组间差异OTU的交集
indic_edge_16s <- intersect(rownames(r_values_16s),rownames(Year_16s))
#基于TMM标准化的OTU表格进行OTU的两两Spearman相关计算
b16s_otu_cor <- rcorr(t(otu_norm_16s),type=c("spearman"))
## 邻接矩阵转化为边表
CorrDF <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    from = rownames(cormat)[col(cormat)[ut]],
    to = rownames(cormat)[row(cormat)[ut]],
    cor =(cormat)[ut],
    p = pmat[ut]
  )
}
b16s_cor_df <- CorrDF(b16s_otu_cor$r,b16s_otu_cor$P)
# p值校正
b16s_cor_df$padj <- p.adjust(b16s_cor_df$p, method = "none") #method可选c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#取Spearman's rho > 0.7且p-value < 0.001的关系作为入选共现性网络co-occurence network的边
b16s_cor_df_padj <- b16s_cor_df[which(b16s_cor_df$cor > 0.6),]
b16s_cor_df_padj <- b16s_cor_df_padj[which(b16s_cor_df_padj$padj < 0.001),]
#生成node属性表
# 边两列合并为结点
nodeattrib_16s <- data.frame(node = union(b16s_cor_df_padj$from,b16s_cor_df_padj$to))
# 显著的添加标签
nodeattrib_16s$indicgroup <- 0
for (i in as.character(nodeattrib_16s$node))
{
  if (i %in% indic_edge_16s == TRUE)
  {nodeattrib_16s[nodeattrib_16s$node==i,"indicgroup"] <- paste(colnames(r_values_16s)[which(r_values_16s[i,1:2]==1)],collapse = "_")}
  else
  { nodeattrib_16s[nodeattrib_16s$node==i,"indicgroup"] <- "NA"}
}
#将OTU，即节点分类信息添加到node属性表
nodeattrib_16s <- cbind(nodeattrib_16s,tax_16s[as.character(nodeattrib_16s$node),])
#用igraph绘制共现性网络图co-occurence network
net_16s <- graph_from_data_frame(b16s_cor_df_padj,direct=F, vertices = nodeattrib_16s)
# comps <- components(t_net_16s)
# comps
# low_deg <- which(comps$csize[comps$membership] < 3)
# low_deg
# t_net_16s <- delete_vertices(t_net_16s, low_deg)
## 网络中的节点相对丰度
ra_16s <- apply(otu_norm_16s,1,mean)
ra_16s <- ra_16s[V(net_16s)$name]
#将上述显著组间差异的OTU着色，这里有CPT,ZT,PT三种处理，因此理论上有6种可能的差异丰度分布情况
cs <- c("Y02Y05","Y02Y05_Y08Y11","Y08Y11")
unique(V(net_16s)$indicgroup)
V(net_16s)$color <- V(net_16s)$indicgroup
V(net_16s)$color[!V(net_16s)$color %in% cs] <- "#adadad"
V(net_16s)$color[V(net_16s)$color == "Y02Y05"] <- "#228b22"
V(net_16s)$color[V(net_16s)$color == "Y02Y05_Y08Y11"] <- "#00FFFF"
V(net_16s)$color[V(net_16s)$color == "Y08Y11"] <- "#1f78b4"
V(net_16s)$frame.color <- V(net_16s)$color
#提取富集信息
richY <- as.data.frame(V(net_16s)$name)
colnames(richY) <- "OTUID"
richY$rich <- V(net_16s)$color
write.table(richY,"richY.txt",sep="\t",row.names=F,col.names=T,quote=F)
#上述着色信息映射到node属性表中
b16s_nodes <- rownames(nodeattrib_16s[nodeattrib_16s$indicgroup %in% cs,])
#设置节点形状
V(net_16s)$shape <- "circle"
##设置节点大小，非显著组间差异OTU设置为"3"，显著的为"5"
V(net_16s)$size <- V(net_16s)$name
V(net_16s)$size[!V(net_16s)$size %in% b16s_nodes] <- 3
V(net_16s)$size[V(net_16s)$size %in% b16s_nodes] <- 5
b16s_nodesizes <- as.numeric(V(net_16s)$size)
#向量化各类型的显著组间差异OTU以便后续计算
Y02Y05_nodes <- rownames(nodeattrib_16s[nodeattrib_16s$indicgroup=="Y02Y05",])
Y02Y05_Y08Y11_nodes <- rownames(nodeattrib_16s[nodeattrib_16s$indicgroup=="Y02Y05_Y08Y11",])
Y08Y11_nodes <- rownames(nodeattrib_16s[nodeattrib_16s$indicgroup=="Y08Y11",])
cs_nodes_all <- c(Y02Y05_nodes,Y02Y05_Y08Y11_nodes, Y08Y11_nodes)
#将网络中的节点/OTU进行聚类，这里采用fast greedy法
cfg <- cluster_fast_greedy(as.undirected(net_16s))
# cfg <- cluster_walktrap(as.undirected(net_16s))
#查看包含OTU数量最多的10个模块，以进行后续的着色
modules <- sort(table(membership(cfg)),decr=T)
modules_10 <- modules[1:10]
sm10_plot <- modules_10
names(sm10_plot) <- as.factor(1:10)
#寻找包含OTU数量最多的10个模块中具有显著组间差异OTU的模块
modules_cs <- table(factor(membership(cfg)[cs_nodes_all],levels=names(modules)))
modules_cs_10 <- modules_cs[names(modules_10)]
smcs10_plot <- modules_cs_10
names(smcs10_plot) <- as.factor(1:10)
#将OTU数量最多的10个模块中的OTU向量化
modules_points <- membership(cfg)[membership(cfg) %in% names(modules_10)]
points <- NULL
for(i in modules_points){
  tx <- which(names(modules_10)==i)
  points <- c(points, tx)
}
names(points) <- names(modules_points)
#按照组间差异OTU的类型着色这些OTU
all_cols <- sort(points)
all_cols[!names(all_cols) %in% cs] <- "#adadad"
all_cols[names(all_cols) %in% Y02Y05_nodes] <- "#228b22"
all_cols[names(all_cols) %in% Y02Y05_Y08Y11_nodes] <- "#00FFFF"
all_cols[names(all_cols) %in% Y08Y11_nodes] <- "#1f78b4"
#设置节点样式，1为空心圆代表普通节点，16为实心圆代表组间差异OTU
all_pch <- sort(points)
all_pch[names(all_pch) %in% rownames(otu_norm_16s)] <- 1
all_pch[names(all_pch) %in% intersect(rownames(otu_norm_16s),cs_nodes_all)] <- 16
#设置节点缩放倍数，1为普通节点1倍，2为组间差异OTU2倍
all_cex <- sort(points)
all_cex[!names(all_cex) %in% cs_nodes_all] <- 1
all_cex[names(all_cex) %in% cs_nodes_all] <- 2
#哪些模块包含有组间差异OTU
mods_list_cs <- list()
for (i in names(modules_cs_10)){
  x1 <- names(membership(cfg)[membership(cfg)==i])
  x2 <- x1[x1 %in% cs_nodes_all]
  mods_list_cs[[i]] <- as.numeric(V(net_16s)[x2])
}
# t_mods_list_cs
#设定layout出图，这里选择Fruchterman & Reingold
set.seed(951001)
coords_16s <- layout_(net_16s,with_fr(niter=9999, grid="nogrid"))
#每次运行耗费时间，可以把文件储存
#write.table(coords_t_16s,"coords_t_16s.txt",sep="\t",row.names=F,col.names=F,quote=F)
#coords_t_16s <- as.matrix(read.table("coords_t_16s.txt"))
#dimnames(coords_t_16s) <- NULL
#出图
# 设置节点的标签属性
V(net_16s)$label <- ifelse(V(net_16s)$name == "ASV_1467", V(net_16s)$name, NA)
write_graph(net_16s, file = "mygraph.graphml", format = "graphml")
pdf(paste0("p1_network.pdf"),width=7,height=5)
par(mfrow=c(1,1), mar=c(0,0,0,0))
cols <- c("#fba801","#87ceeb")
plot(net_16s,vertex.label=V(net_16s)$label,vertex.size=b16s_nodesizes, layout=coords_16s,
     mark.groups=list(mods_list_cs$`1`,mods_list_cs$`2`),
     mark.col=cols, mark.border=cols)
legend("right",legend=c("Module 1","Module 2"),col=cols,
       bty="n",fill=cols,border=cols)
dev.off()

#模块中组内OTU相对丰度比较
pdf(paste0("p2_module_abundance.pdf"),width=7,height=3)
par(mfrow=c(1,6), mar=c(0.5,3.5,2,0))#根据模块数调整mfrow
CS_cols <- c("#F8766D","#619CFF")
names(CS_cols) <- c("Y02Y05","Y08Y11")
#T module 1
bargraph.CI(design_16s$Group, colSums(otu_norm_16s[cfg[[1]],])/1000,
            las=2, ylab="cumulative relative abundance", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 1", col=CS_cols, border=F)
#T module 2
bargraph.CI(design_16s$Group, colSums(otu_norm_16s[cfg[[2]],])/1000,
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 2", col=CS_cols, border=F)
#T module 3
# bargraph.CI(design_16s$Group, colSums(otu_norm_16s[cfg[[3]],])/1000,
#             las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
#             err.width=.025, main="Module 3", col=CS_cols, border=F)
plot.new()
par(mar=c(0.5,0,2,0))
legend("left", bty="n", cex=1, #x.intersp=0.1, y.intersp=1,
       legend=names(CS_cols),
       fill=CS_cols,
       border=CS_cols , xpd=T)
dev.off()
m1 = as.data.frame(colSums(otu_norm_16s[cfg[[1]],])/1000)
colnames(m1) = "Abundance"
m1$Module = "Module_1"
m1$Group = design_16s$Group
m2 = as.data.frame(colSums(otu_norm_16s[cfg[[2]],])/1000)
colnames(m2) = "Abundance"
m2$Module = "Module_2"
m2$Group = design_16s$Group
# m3 = as.data.frame(colSums(otu_norm_16s[cfg[[3]],])/1000)
# colnames(m3) = "Abundance"
# m3$Module = "Module_3"
# m3$Group = design_16s$Group
#将三个数据框按行合并
module = rbind(m1,m2)
write.csv(module, "module.csv")
module <- module %>% group_by(Module,Group) %>% 
  summarise(mean = mean(Abundance), se = sd(Abundance)/sqrt(n()))
#绘制分面柱形图
p <- ggplot(module, aes(x = Group, y = mean, fill = Group))+
  geom_col(colour = "black",width = 0.6)+
  scale_y_continuous(expand = c(0, 0),limits = c(0,370))+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width=.2, position = position_dodge(.9)) +
  labs(y = "Cumulative relative abundance", x="") +
  facet_wrap(.~Module)+
  theme_test()+
  theme(strip.background = element_rect(fill = "lightblue"), 
        strip.text = element_text(size = 20),
        legend.title = element_blank(),
        text = element_text(colour='black',size=20), 
        axis.text=element_text(colour='black',size=17),legend.position = "none")+
  scale_fill_manual(values = c("#228b22", "#1f78b4"))

p
# 提取每个ASV对应的模块
ext_module <- as.data.frame(cfg$names)
colnames(ext_module) <- c("ASV")
ext_module$module <- cfg$membership
write.table(ext_module,"ext_module.txt",sep="\t",row.names=F,col.names=T,quote=F)
#分析zipi
set.seed(12)
#'net_16s'是网络对象
communities <- cluster_fast_greedy(net_16s)
comm_membership <- membership(communities)
# 计算每个节点的度
node_degree <- degree(net_16s)
# 计算每个模块的平均度和标准差
module_mean <- tapply(node_degree, comm_membership, mean)
module_sd <- tapply(node_degree, comm_membership, sd)
# 计算Zi，并处理标准差为0的情况
Zi <- ifelse(module_sd[comm_membership] != 0,
             (node_degree - module_mean[comm_membership]) / module_sd[comm_membership],
             0)
# 创建一个矩阵来存储每个节点与每个模块的连接数
module_connections <- matrix(0, nrow=length(V(net_16s)), ncol=max(comm_membership))
# 计算每个节点与每个模块的连接数
for (v in V(net_16s)) {
  neighbors_of_v <- neighbors(net_16s, v)
  modules_of_neighbors <- comm_membership[neighbors_of_v]
  module_connections[v, ] <- table(factor(modules_of_neighbors, levels=1:max(comm_membership)))
}
# 计算Pi
Pi <- 1 - rowSums((module_connections / node_degree)^2)
# 将结果整合到一个数据框
zi_pi_metrics <- data.frame(
  name = V(net_16s)$name,
  Zi = Zi,
  Pi = Pi
)
# 节点分类
zi_pi_metrics$type <- ifelse(zi_pi_metrics$Zi > 2.5 & zi_pi_metrics$Pi < 0.62, "Module hubs",
                             ifelse(zi_pi_metrics$Zi < 2.5 & zi_pi_metrics$Pi > 0.62, "Connectors",
                                    ifelse(zi_pi_metrics$Zi > 2.5 & zi_pi_metrics$Pi > 0.62, "Network hubs",
                                           "Peripherals")))
# 导出CSV
write.csv(zi_pi_metrics, "zi_pi_metrics.csv", row.names = FALSE)
p <- ggplot(zi_pi_metrics, aes(x = Pi, y = Zi, color = type)) +
  geom_point(alpha = 0.7, size = 3) +
  theme_minimal() +
  labs(
    x = "Among-module conectivities (Pi)",
    y = "Within-module conectivities (Zi)",
    color = "Node Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  geom_vline(xintercept = 0.62, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 2.5, linetype = "dashed", alpha = 0.5)

p


# Fig 4g line-
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
p2 <- ggplot(lmm, aes(x = asv, y = meloi)) + geom_point(aes(colour = Modules,shape = Years), size = 3) +
  geom_smooth(method = "lm", aes(colour = Modules, fill = Modules)) +
  scale_colour_manual(values = c("#fba801","#87ceeb"))+
  scale_fill_manual(values = c("#fba801","#87ceeb"))+
  stat_cor(aes(color=Modules),method = "pearson",label.x = -1.7, size = 5) +
  theme_test() +
  theme(text = element_text(colour='black',size=17), axis.text=element_text(colour='black',size=13),
        legend.text = element_text(size = 15)) +
  labs(y = "Relative abundance of meloidogyne \n (log-transformed)", x="Relative abundance of module \n (log-transformed)")
p2  
p3 <- ggplot(lmm, aes(x = asv, y = nematode)) + geom_point(aes(colour = Modules,shape = Years), size = 3) +
  geom_smooth(method = "lm", aes(colour = Modules, fill = Modules)) +
  scale_colour_manual(values = c("#fba801","#87ceeb"))+
  scale_fill_manual(values = c("#fba801","#87ceeb"))+
  stat_cor(aes(color=Modules),method = "pearson",label.x = -1.7, size = 5) +
  theme_test() +
  theme(text = element_text(colour='black',size=17), axis.text=element_text(colour='black',size=13),
        legend.text = element_text(size = 15)) +
  labs(y = "Soil nematode abundance \n (per 100 g dry soil)", x="Relative abundance of module \n (log-transformed)")
p3
p4 <- ggarrange(p2,p3,nrow = 1,common.legend = TRUE,legend="right")
p4



