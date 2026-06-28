library(readxl)
library(dplyr)
library(ggplot2)
library(forcats)

# 读文件，替换成你的文件名
file_cc <- "C:/Users/20379/Desktop/新建文件夹/1-1-CC.xlsx"
file_mf <- "C:/Users/20379/Desktop/新建文件夹/1-1-MF.xlsx"
file_bp <- "C:/Users/20379/Desktop/新建文件夹/1-1-BP.xlsx"

df_cc <- read_excel(file_cc)
df_mf <- read_excel(file_mf)
df_bp <- read_excel(file_bp)

# 给类别列
df_cc$Category <- "CC"
df_mf$Category <- "MF"
df_bp$Category <- "BP"

# 合并
df_all <- bind_rows(df_cc, df_mf, df_bp)
  
# 筛选每个类别PValue值最小的7个
topn <- 10
df_top <- df_all %>%
  group_by(Category) %>%
  arrange(PValue) %>%
  slice_head(n = topn) %>%
  ungroup()

# 计算-log10(PValue)
df_top <- df_top %>%
  mutate(logBen = -log10(PValue))

# 保证Term顺序合理（按类别和logBen排序）
df_top <- df_top %>%
  arrange(Category, logBen) %>%
  mutate(Term = factor(Term, levels = unique(Term)))

df_top$Category <- factor(df_top$Category, levels = c("MF", "CC", "BP")) 
# 画图
# p<- ggplot(df_top, aes(x = Term, y = logBen, fill = Category)) +
#   geom_col(position = position_dodge(width = 0.8)) +
#   coord_flip() +
#   labs(x = "GO Term", y = "-log10(PValue)", fill = "Category") +
#   theme_minimal() +
#   theme(axis.text.y = element_text(size = 15))



p <- ggplot(df_top) +
  # 线段
  geom_segment(
    aes(
      x = Term, 
      xend = Term, 
      y = 0, 
      yend = logBen, 
      color = Category
    ),
    linewidth = 2, lineend = "round"
  ) +
  # 圆点（在线段末端）
  geom_point(
    aes(
      x = Term, 
      y = logBen, 
      color = Category
    ),
    size = 4
  ) +
  coord_flip(expand = TRUE) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),  # 保留原有的扩展规则（可根据需要调整）
    limits = c(0, 8)  # 新增：将y轴范围设为 0 到 6
  )+
  scale_x_discrete(
    expand = expansion(add = c(3, 1))  # 关键：顶部强制留白3个单位
  ) +
  scale_color_manual(
    values = c(
      "MF" = hcl(h = 210, c = 40, l = 75), # 蓝
      "CC" = hcl(h = 55, c = 55, l = 85) ,  # 马卡龙淡黄色
      "BP" = hcl(h = 130, c = 40, l = 75)  # 绿
    )
    
  ) +
  scale_x_discrete(expand = expansion(add = c(1, 1))) +
  labs(
    x = "GO Term", 
    y = "-log10(PValue)", 
    color = "Category"
  ) +
  theme_minimal(base_size = 30) +
  theme(
    axis.text.y = element_text(size = 33),
    axis.text.x = element_text(size = 25),
    axis.title  = element_text(size = 25),
    plot.title  = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 30),
    legend.title= element_text(size = 30),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, b = 20, l = 10, r = 20)
  )

# 保存为tiff格式
ggsave("C:/Users/20379/Desktop/新建文件夹/预实验GO_enrichment_barplot.tiff", plot = p, width = 25, height = 15, dpi = 300)



