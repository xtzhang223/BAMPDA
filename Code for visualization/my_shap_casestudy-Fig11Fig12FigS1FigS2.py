import pandas as pd
import numpy as np
import xgboost as xgb
import shap
import matplotlib.pyplot as plt
import os

# --- 0. 全局绘图设置 ---
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['font.size'] = 18
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

# --- 1. 路径配置 ---
base_path = r"C:\Users\zxt\Desktop\run_MDHGIBNNR\1036-390\five_fold\dataset"
predict_score_path = r"C:\Users\zxt\Desktop\run_MDHGIBNNR\1036-390\paint-roc-data1\Predict_Score_Matrix3_3.xlsx"

files = {
    'A': os.path.join(base_path, "Protein_Disease_Associations.xlsx"),
    'S_seq': os.path.join(base_path, "Protein_Sequence_Similarity_Matrix.xlsx"),
    'S_pg': os.path.join(base_path, "Protein_Gaussian_Similarity_Matrix.xlsx"),
    'S_sem': os.path.join(base_path, "Disease_Semantic_Similarity_Matrix.xlsx"),
    'S_dg': os.path.join(base_path, "Disease_Gaussian_Similarity_Matrix.xlsx")
}

feature_display_names = [
    'Known Protein-Disease Associations',      # f1
    'Protein Sequence Profile',                # f2
    'Protein Gaussian Interaction Profile',    # f3
    'Disease Semantic Profile',                # f4
    'Disease Gaussian Interaction Profile'     # f5
]

# --- 2. 数据处理函数 ---
def load_and_clean_data():
    print("正在加载数据...")
    mats = {}
    for k, v in files.items():
        df = pd.read_excel(v, header=0, index_col=0)
        mats[k] = df.values.astype(float)
    df_score = pd.read_excel(predict_score_path, header=0, index_col=0)
    p_score = df_score.values.astype(float)
    if p_score.shape != mats['A'].shape:
        p_score = p_score[:mats['A'].shape[0], :mats['A'].shape[1]]
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

# --- 3. 主程序 ---
try:
    # A. 准备数据和全局模型
    matrices, P_score = load_and_clean_data()
    
    print("1. 正在训练全局解释器 (约需 30 秒)...")
    np.random.seed(42)
    rows, cols = matrices['A'].shape
    sample_indices = np.random.choice(rows * cols, 100000, replace=False)
    X_train, y_train = [], []
    for idx in sample_indices:
        r, c = divmod(idx, cols)
        X_train.append(extract_features(r, c, matrices))
        y_train.append(P_score[r, c])
    
    model = xgb.XGBRegressor(n_estimators=100, max_depth=4, learning_rate=0.05)
    model.fit(np.array(X_train), np.array(y_train))
    explainer = shap.Explainer(model, np.array(X_train))
    
    # ======================================================
    #      B. Case Study: 智能选择并绘图
    # ======================================================
    
    # 设定你要分析的疾病列索引 (0 代表第一列 CML)
    TARGET_COL_INDEX = 4  
    DISEASE_NAME = "HUNTINGTON’S DISEASE"
    
    print(f"\n2. 开始 Case Study 分析: {DISEASE_NAME} ...")
    
    # 提取该列预测分并找到 Top 20
    scores = P_score[:, TARGET_COL_INDEX]
    top20_indices = np.argsort(scores)[-20:][::-1] # 从高到低排序
    
    # 提取特征并计算 SHAP
    X_target = []
    for p_idx in top20_indices:
        X_target.append(extract_features(p_idx, TARGET_COL_INDEX, matrices))
    X_target = np.array(X_target)
    
    shap_vals = explainer(X_target)
    shap_vals.feature_names = feature_display_names
    
    # --- 扫描前 5 名，寻找最佳绘图对象 ---
    print("\n--- 正在检查前 5 名候选者的特征贡献详情 ---")
    
    best_candidate_idx = 0       # 默认选第1名
    max_active_features = 0      # 记录哪个候选者的非零特征最多
    
    for i in range(5):
        vals = shap_vals[i].values
        active_count = np.sum(np.abs(vals) > 0.001)
        
        print(f"Rank {i+1} (Protein Index {top20_indices[i]}): 有 {active_count} 个显著特征")
        for name, val in zip(feature_display_names, vals):
            print(f"   - {name[:20]}... : {val:.4f}")
            
        if active_count > max_active_features:
            max_active_features = active_count
            best_candidate_idx = i
            
    print("-" * 40)
    print(f"最终决定绘制: Rank {best_candidate_idx+1} (Protein Index: {top20_indices[best_candidate_idx]})")
    

    # ======================================================
    #   绘图 2: 热力图 (Top 20 整体) - 严格锁定字体格式版
    # ======================================================
    print("\n3. 正在绘制热力图并调整格式...")
    fig = plt.figure(figsize=(16, 9)) 
    
    # 画图，但不要立刻显示 (show=False)
    shap.plots.heatmap(shap_vals, instance_order=shap_vals.sum(1), show=False, plot_width=16)

    # 获取图表中所有的坐标轴对象 (包括主图、顶部曲线、颜色条)
    all_axes = plt.gcf().get_axes()
    
    # 1. 设置主图坐标轴 (通常是 axes[0])
    ax_main = all_axes[0]
    
    # 强制修改 X 轴 (横轴)
    ax_main.tick_params(axis='x', labelsize=18)
    for tick in ax_main.get_xticklabels():
        tick.set_fontname("Times New Roman")
    ax_main.set_xlabel("Top 20 Candidate Proteins", fontsize=18, fontfamily='serif', name='Times New Roman', fontweight='bold', labelpad=15)
    
    # 强制修改 Y 轴 (纵轴：特征名称)
    ax_main.tick_params(axis='y', labelsize=18)
    for tick in ax_main.get_yticklabels():
        tick.set_fontname("Times New Roman")
    

    # 2. 遍历并修改所有子图的刻度字体 
    for temp_ax in all_axes:
        temp_ax.tick_params(labelsize=18)
        
        # [核心修复]：把 SHAP value 的标签往外推
        current_ylabel = temp_ax.get_ylabel()
        if current_ylabel == "SHAP value":
            # labelpad=30 强行增加间距
            temp_ax.set_ylabel("SHAP value", labelpad=30, fontsize=18, fontfamily='serif', name='Times New Roman')
        elif current_ylabel and temp_ax != ax_main:
            temp_ax.yaxis.labelpad = 15 

        # 获取所有可能的文本对象并强制修改为 Times New Roman
        for item in ([temp_ax.xaxis.label, temp_ax.yaxis.label] +
                     temp_ax.get_xticklabels() + temp_ax.get_yticklabels()):
            item.set_fontsize(18)
            item.set_fontname("Times New Roman")
            
        # 如果子图自带标题(例如颜色条旁边的说明)，也改掉
        if temp_ax.get_title():
            temp_ax.title.set_fontsize(18)
            temp_ax.title.set_fontname("Times New Roman")

    # 3. 添加总标题 (使用 suptitle 防止被 SHAP 的内部排版挤掉，y=1.05 稍微往上抬一点)
    plt.suptitle("Feature Contribution Patterns for HUNTINGTON’S DISEASE Candidates", 
                 fontsize=18, fontfamily='serif', name='Times New Roman', fontweight='bold', y=1.05)

    # 4. 调整边界
    # left=0.4 给左侧特征名留空间，right=0.85 给右侧被推开的 SHAP value 留空间
    plt.subplots_adjust(left=0.4, bottom=0.2, top=0.9, right=0.85)

    # 5. 高清保存 (推荐 600 DPI，符合顶刊标准)
    output_filename = "CaseStudy_Heatmap.png"
    plt.savefig(output_filename, dpi=600, bbox_inches='tight')
    plt.show()
    print(f"✅ 已成功保存高分辨率格式要求图片: {output_filename}")

except Exception as e:
    import traceback
    print(f"发生错误: {e}")
    traceback.print_exc()