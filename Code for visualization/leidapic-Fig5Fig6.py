import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# 模型名称和指标
models = ['BAMPDA', 'BAMPDA-MDHGI', 'BAMPDA-BNNR']
categories = ['F1-score', 'Precision', 'Recall', 'AUC']
data = np.array([
    [0.9092, 0.8829, 0.8800, 0.8866],  # BAMPDA
    [0.8815, 0.7237, 0.8138, 0.8815],  # BAMPDA-MDHGI
    [0.6766, 0.4920, 0.6906, 0.6766],  # BAMPDA-BNNR
])

# 设置字体为 Times New Roman（全局）
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 20  # 全局字体变大

# 设置雷达图角度和闭合循环
angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False).tolist()
angles += angles[:1]

# 初始化雷达图
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'polar': True})

# 绘制每个模型的数据
for i in range(len(models)):
    values = data[i].tolist()
    values += values[:1]
    ax.plot(angles, values, linewidth=2, label=models[i])
    ax.fill(angles, values, alpha=0.25)

# 设置刻度标签
ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories, fontsize=20, fontname='Times New Roman')
ax.tick_params(axis='x', which='major', pad=20)

# 设置径向坐标轴
ax.set_rlabel_position(30)
plt.yticks([0.45, 0.55, 0.65, 0.75, 0.85, 0.95],
           color="grey", size=20, fontname='Times New Roman')
plt.ylim(0.45, 0.95)

# 图例
plt.legend(loc='upper right', bbox_to_anchor=(1.2, 1.17), fontsize=18)

# 桌面路径
desktop_path = Path.home() / "Desktop"
desktop_path.mkdir(parents=True, exist_ok=True)

# 保存为 JPG
file_path = desktop_path / "radar_plot.jpg"
plt.savefig(file_path, dpi=400, format='jpg', bbox_inches='tight')
print(f"图片已保存到: {file_path}")

plt.show()
