import pandas as pd
import numpy as np
import xgboost as xgb
import shap
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import FuncFormatter

# --- 1. 路径配置 ---
base_path = r"C:\Users\zxt\Desktop\run_MDHGIBNNR\1036-390\five_fold\dataset"
predict_score_path = r"C:\Users\zxt\Desktop\run_MDHGIBNNR\1036-390\paint-roc-data2\Predict_Score_Matrix3_3.xlsx"

files = {
    'A': os.path.join(base_path, "Protein_Disease_Associations.xlsx"),
    'S_seq': os.path.join(base_path, "Protein_Sequence_Similarity_Matrix.xlsx"),
    'S_pg': os.path.join(base_path, "Protein_Gaussian_Similarity_Matrix.xlsx"),
    'S_sem': os.path.join(base_path, "Disease_Semantic_Similarity_Matrix.xlsx"),
    'S_dg': os.path.join(base_path, "Disease_Gaussian_Similarity_Matrix.xlsx")
}

# --- 2. 定义特征名称 ---
feature_display_names = [
    'Known Protein-Disease Associations',
    'Protein Sequence Profile',
    'Protein Gaussian Interaction Profile',
    'Disease Semantic Profile',
    'Disease Gaussian Interaction Profile'
]

# --- 3. 柱状图专用的科学计数法格式化函数 ---
def bar_sci_formatter(x, pos):
    if x == 0: return "0"
    # 将数值放大1000倍，保留1位小数，拼接上 \times 10^-3
    return f"{x*1000:.1f}$\\times$10$^{{-3}}$"

# --- 4. 数据加载与处理函数 ---
def load_and_clean_data():
    print("正在加载数据...")
    mats = {}
    for k, v in files.items():
        df = pd.read_excel(v, header=0, index_col=0)
        mats[k] = df.values.astype(float)
    df_score = pd.read_excel(predict_score_path, header=0, index_col=0)
    p_score = df_score.values.astype(float)
    expected_shape = mats['A'].shape
    if p_score.shape != expected_shape:
        p_score = p_score[:expected_shape[0], :expected_shape[1]]
    return mats, p_score

def extract_features(p_idx, d_idx, mats):
    A = mats['A']
    assoc_p = np.where(A[:, d_idx] > 0)[0] 
    assoc_d = np.where(A[p_idx, :] > 0)[0] 
    f1 = len(assoc_d)                                                          
    f2 = np.mean(mats['S_seq'][p_idx, assoc_p]) if len(assoc_p) > 0 else 0     
    f3 = np.mean(mats['S_pg'][p_idx, :])                                       
    f4 = np.mean(mats['S_sem'][d_idx, assoc_d]) if len(assoc_d) > 0 else 0     
    f5 = np.mean(mats['S_dg'][d_idx, :])                                       
    return [f1, f2, f3, f4, f5]

# --- 5. 主程序 ---
try:
    matrices, P_score = load_and_clean_data()
    
    print("正在采样 100,000 个样本进行分析...")
    np.random.seed(42)
    rows, cols = matrices['A'].shape
    sample_indices = np.random.choice(rows * cols, 100000, replace=False)
    
    X_data, y_target = [], []
    for idx in sample_indices:
        r, c = divmod(idx, cols)
        X_data.append(extract_features(r, c, matrices))
        y_target.append(P_score[r, c])
        
    X_data = np.array(X_data)
    y_target = np.array(y_target)
    
    model = xgb.XGBRegressor(n_estimators=100, max_depth=4, learning_rate=0.05)
    model.fit(X_data, y_target)
    
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_data)

    # --- 绘图 1: 柱状图 (采用科学计数法) ---
    plt.figure(figsize=(11, 6))
    shap.summary_plot(shap_values, X_data, feature_names=feature_display_names, plot_type="bar", show=False)
    
    plt.xlabel("mean(|SHAP value|)", fontsize=14, fontname='Times New Roman')
    plt.title("Feature Importance Ranking for Dataset1", fontsize=16, fontname='Times New Roman')
    
    # 应用科学计数法格式化器
    plt.gca().xaxis.set_major_formatter(FuncFormatter(bar_sci_formatter))
    
    plt.yticks(fontsize=12, fontname='Times New Roman')
    plt.xticks(fontsize=10, fontname='Times New Roman')
    plt.tight_layout()
    plt.savefig("Feature_Bar_Scientific.png", dpi=300)
    print("柱状图已保存。")
    plt.show()

    # --- 绘图 2: 蜂群图 (不采用科学计数法，保持默认) ---
    plt.figure(figsize=(11, 6))
    shap.summary_plot(shap_values, X_data, feature_names=feature_display_names, show=False)
    
    plt.xlabel("SHAP value", fontsize=14, fontname='Times New Roman')
    plt.title("Detailed Feature Impact on Prediction Score for Dataset1", fontsize=16, fontname='Times New Roman')
    
    # 蜂群图不设置 set_major_formatter，保持默认的小数显示
    
    plt.yticks(fontsize=12, fontname='Times New Roman')
    plt.xticks(fontsize=12, fontname='Times New Roman')
    plt.tight_layout()
    plt.savefig("Feature_Beeswarm_Standard.png", dpi=300)
    print("蜂群图已保存。")
    plt.show()

except Exception as e:
    print(f"发生错误: {e}")