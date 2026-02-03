# 加载必要的包
library(ggplot2)

# 读取数据
data <- read.csv("bacillibactin.csv", header = TRUE, strip.white = TRUE)

# 计算统计量
mean_val <- mean(data$mean)
sd_val <- sd(data$mean)

# 创建汇总数据框
summary_data <- data.frame(
  Group = "Y11.1",
  Mean = mean_val,
  SD = sd_val
)

# 为原始数据添加分组
data$Group <- "Y11.1"

# 绑图
p <- ggplot() +
  # 柱状图
  geom_bar(data = summary_data, aes(x = Group, y = Mean), 
           stat = "identity", fill = "grey75", color = "black", width = 0.5) +
  # 误差棒
  geom_errorbar(data = summary_data, aes(x = Group, ymin = Mean - SD, ymax = Mean + SD),
                width = 0.15, linewidth = 0.8) +
  # 离散点
  geom_jitter(data = data, aes(x = Group, y = mean), 
              width = 0.1, size = 2, alpha = 0.7) +
  # Y轴标签
  labs(y = "Siderophore concentration in supernate\n(μM equivalents of DFOB)",
       x = "") +
  # Y轴范围
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 20), expand = c(0, 0)) +
  # 主题设置
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 11),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )

# 显示图形
print(p)