import matplotlib.pyplot as plt

# ==========================================
# 全局字体设置 (Times New Roman)
# ==========================================
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']

# 1. 数据准备
# 顺序依次为：Phosphorylation, Acetylation, Methylation, Glycosylation, Ubiquitination

labels = ["Phosphorylation", "Acetylation", "Methylation", "Glycosylation", "Ubiquitination"]
counts = [600, 42, 131, 100, 127]

# 2. 设定低饱和度科研配色 (莫兰迪高级色系)
color_dict = {
    "Phosphorylation": "#4DBBD5", # 莫兰迪蓝
    "Acetylation":     "#D2C6A9", # 莫兰迪沙/黄
    "Methylation":     "#00A087", # 莫兰迪红/灰粉
    "Glycosylation":   "#3C5488", # 莫兰迪紫
    "Ubiquitination":  "#F39B7F"  # 莫兰迪绿/豆绿
}
colors = [color_dict[label] for label in labels]

# 3. 创建图表
fig, ax = plt.subplots(figsize=(8, 6)) # 设置画布大小

# 绘制饼图
wedges, texts, autotexts = ax.pie(
    counts, 
    colors=colors,
    autopct='%1.1f%%',       # 自动计算并显示百分比，保留1位小数
    startangle=90,           # 从90度（正上方）开始绘制
    pctdistance=0.65,        # 百分比数值距离圆心的距离
    wedgeprops={
        'edgecolor': 'white', # 扇区之间的白色极简描边
        'linewidth': 1.5      # 描边粗细
    },
    textprops={
        'fontsize': 24,       # 【修改这里：调大百分比数值的字体大小】
        'fontweight': 'bold', # 数值加粗
        'color': '#333333'    # 字体颜色为深灰/黑色
    }
)

# 4. 设置图例 (Legend)
legend = ax.legend(
    wedges, labels,
    title="PTM Types",
    loc="center left",
    bbox_to_anchor=(1, 0, 0.5, 1), # 将图例移到饼图的右侧外面
    fontsize=30                   # 【修改这里：调大图例具体 PTM 标签的字体】
)

# 【调大图例标题的字体】
legend.get_title().set_fontsize(30)
legend.get_title().set_fontweight('bold')

# 自动调整布局，防止右侧图例被截断
plt.tight_layout()

# 5. 保存为高质量出版级图片 (300 DPI)
plt.savefig("PTM_Distribution_Morandi.pdf", dpi=300, bbox_inches='tight')
plt.savefig("PTM_Distribution_Morandi.png", dpi=300, bbox_inches='tight')
print("图像已成功生成，并保存为 PDF 和 PNG 格式！")

# 6. 在屏幕上显示图表
plt.show()