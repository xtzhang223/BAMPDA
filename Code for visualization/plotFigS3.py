import numpy as np
import matplotlib.pyplot as plt

# Data
methods = [
    "Phosphorylation",
    "Ubiquitination",
    "Methylation",
    "Acetylation",
    "Glycosylation",
]

auc = [0.9862, 0.9885, 0.9868, 0.9650, 0.8595]
f1 = [0.8993, 0.8679, 0.8904, 0.8908, 0.8222]

x = np.arange(len(methods))
width = 0.35

# Figure
plt.figure(figsize=(8, 5.2), dpi=120)
ax = plt.gca()

bars1 = ax.bar(
    x - width / 2,
    auc,
    width,
    label="AUC",
    color="#4C90D9",
    edgecolor="black",
    linewidth=0.6,
)

bars2 = ax.bar(
    x + width / 2,
    f1,
    width,
    label="F1-score",
    color="#F39C12",
    edgecolor="black",
    linewidth=0.6,
)

# Axis
ax.set_ylabel("Performance Score", fontweight="bold")
ax.set_ylim(0.68, 1.17)
ax.set_yticks(np.arange(0.7, 1.2, 0.1))

ax.set_xticks(x)
ax.set_xticklabels(methods, rotation=35, ha="right", fontweight="bold")

# Legend
ax.legend(frameon=False, loc="upper right")

# Value labels
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height + 0.008,
            f"{height:.4f}",
            ha="center",
            va="bottom",
            rotation=90,
            fontsize=8,
        )

# Style
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(direction="out")
ax.grid(False)

plt.tight_layout()

# Save and show
plt.savefig("performance_barplot.png", dpi=600, bbox_inches="tight")
plt.show()