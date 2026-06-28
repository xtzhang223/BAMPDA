library(readxl) 
library(dplyr) 
library(tidyr) 
library(stringr) 
library(writexl) 
library(circlize)
library(viridis)
library(ggsci)
library(colorspace)
# ====================== Step 1：读取 KEGG 结果 =======================
input_excel <- "C:/Users/20379/Desktop/新建文件夹/1-2-KEGG.xlsx"
output_table <- "C:/Users/20379/Desktop/新建文件夹/C_input.xlsx"
output_pdf   <- "C:/Users/20379/Desktop/新建文件夹/KEGG_chord.pdf"

# ====================== Step 2：筛选 Top N 通路 =======================
top_n <- 20
kegg_raw <- read_excel(input_excel, sheet = 1)

kegg_top <- kegg_raw %>%
  arrange(PValue) %>%
  slice_head(n = top_n)

# ====================== Step 3：拆分 Genes 字段 =======================
kegg_long <- kegg_top %>%
  select(Term, Genes) %>%
  mutate(Genes = str_replace_all(Genes, " ", "")) %>%
  separate_rows(Genes, sep = ",")

write_xlsx(kegg_long, output_table)


# ====================== Step 4：绘制 KEGG 弦图为 PDF =======================
terms <- unique(kegg_long$Term)


# 为每个色系定义色相（h），中等饱和度（c），中亮度（l）
blues       <- hcl(h = 210, c = 60, l = seq(65, 75, length.out = 2)) # 蓝
cyan_green  <- hcl(h = 170, c = 60, l = seq(65, 75, length.out = 2)) # 青绿
grass_green <- hcl(h = 130, c = 60, l = seq(65, 75, length.out = 2)) # 草绿
beige       <- hcl(h = 50,  c = 35, l = seq(80, 85, length.out = 2)) # 米黄
orange      <- hcl(h = 35,  c = 60, l = seq(68, 78, length.out = 2)) # 橙
light_brown <- hcl(h = 30,  c = 40, l = seq(70, 78, length.out = 2)) # 浅棕
mint_green  <- hcl(h = 150, c = 50, l = seq(70, 80, length.out = 2)) # 薄荷绿
pink        <- hcl(h = 340, c = 45, l = seq(75, 82, length.out = 2)) # 粉
light_purple<- hcl(h = 280, c = 45, l = seq(75, 82, length.out = 2)) # 浅紫
sky_blue    <- hcl(h = 200, c = 50, l = seq(75, 82, length.out = 2)) # 天蓝
lemon_yellow<- hcl(h = 60,  c = 60, l = seq(78, 85, length.out = 2)) # 柠檬黄
honey_orange<- hcl(h = 40,  c = 55, l = seq(75, 82, length.out = 2)) # 蜜橙

# 合并所有色系（每个取 1-2 个色），并截取前 20 个
palette_20 <- c(
  blues, cyan_green, grass_green, beige, orange, light_brown,
  mint_green, pink, light_purple, sky_blue, lemon_yellow, honey_orange
)[1:20]

# 假设 palette_20 是之前生成的 20 种颜色
set.seed(123)  # 固定随机顺序，方便复现
palette_20 <- sample(palette_20)  # 打乱颜色顺序

# 绑定到 terms
term_colors <- setNames(
  rep(palette_20, length.out = length(terms)),
  terms
)








# term_colors <- c( install.packages("viridis")
# "Human cytomegalovirus infection" = "#E18182", 
# "Endocrine resistance" = "#CADDE3", 
# "ErbB signaling pathway" = "#BBCFC3",
# "PD-L1 expression and PD-1 checkpoint pathway in cancer" = "#D9B9D4", 
# "Proteoglycans in cancer" = "#7C9895", 
# "Chemical carcinogenesis - receptor activation" = "#F4EEAC", 
# "Breast cancer" = "#DAA87C", # "Shigellosis" = "#FBDFE2", 
# "EGFR tyrosine kinase inhibitor resistance" = "#B2D3A4", 
# "Pathways in cancer" = "#B696B6", 
# "Choline metabolism in cancer" = "#FC8D62", 
# "Human papillomavirus infection" = "#8DA0CB", 
# "Leukocyte transendothelial migration" = "#599CB4",
# "Growth hormone synthesis, secretion and action" = "#B2DF8A", 
# "Adipocytokine signaling pathway" = "#FB9A99", 
# "JAK-STAT signaling pathway" = "#FDBF6F", 
# "Pancreatic cancer" = "#CAB2D6", 
# "Longevity regulating pathway" = "#FFFF99", 
# "Focal adhesion" = "#B15928", 
# "Human immunodeficiency virus 1 infection" = "#81B5D5" # ) 
# 打开 PDF 设备 
pdf(file = output_pdf, 
    width = 65, 
    height = 55,  
    bg = "transparent") 
chordDiagram(kegg_long, 
             grid.col = term_colors, 
             transparency = 0.25, 
             annotationTrack = "grid", 
             preAllocateTracks = 1) 
circos.trackPlotRegion(track.index = 1, 
                       panel.fun = function(x, y) { 
                         circos.text(CELL_META$xcenter, 
                                     CELL_META$ylim[1], 
                                     CELL_META$sector.index, 
                                     facing = "clockwise", 
                                     niceFacing = TRUE, 
                                     adj = c(0, 0.5), 
                                     cex = 7.8) 
                         }, bg.border = NA) 
dev.off()