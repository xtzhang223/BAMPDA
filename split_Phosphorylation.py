import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, recall_score, precision_score, f1_score

# ==========================================
# 1. 配置文件路径和目标蛋白质名单
# ==========================================
# 使用 r"" 防止 Windows 路径转义报错
true_matrix_path = r"C:\Users\zxt\Desktop\run_MDHGIBNNR\749-276\five_fold\dataset\Protein_Disease_Associations.xlsx"
pred_matrix_path = r"C:\Users\zxt\Desktop\run_MDHGIBNNR\749-276\paint-roc-data1\Predict_Score_Matrix3_3.xlsx"

# 磷酸化相关的蛋白质名单 (非常长，包含多项重复)

target_proteins = [
"ogt", "cd55", "mapt", "hnrnph1", "rasgef1b", "fgg", "shbg", "rnase1", "f5",
"cftr", "prph2", "gabrb3", "fbn1", "ifngr2", "pdia4", "b4gat1", "sgce", "f8",
"serpinc1", "fga", "orm1", "opn1mw", "app", "itgb2", "serping1", "ctsb", "sp1",
"igf1r", "gusb", "met", "ctsa", "mapt", "prf1", "muc1", "epcam", "apex1",
"cd40lg", "il2rg", "psen1", "tusc3", "dag1", "dag1", "erap2", "stt3b", "muc16",
"gpa33", "chordc1", "mok", "ryr1", "zdhhc14", "app", "bdnf", "pafah1b3", "slc6a4"
]

# 将目标蛋白质名字全部转换为大写，并使用 set() 自动去重，提高匹配效率
target_proteins = list(set([p.upper() for p in target_proteins]))

print("正在加载数据，请稍候...")

# ==========================================
# 2. 读取 Excel 文件
# ==========================================
true_df = pd.read_excel(true_matrix_path, index_col=0)
pred_df = pd.read_excel(pred_matrix_path, index_col=0)

# 将矩阵的行索引（蛋白质名称）也全部转换为大写，并去除可能存在的首尾空格
true_df.index = true_df.index.astype(str).str.strip().str.upper()
pred_df.index = pred_df.index.astype(str).str.strip().str.upper()

print(f"✅ 成功加载数据！真实矩阵大小: {true_df.shape}, 预测矩阵大小: {pred_df.shape}")

# 检查矩阵的形状和索引是否对齐
if not true_df.index.equals(pred_df.index) or not true_df.columns.equals(pred_df.columns):
    print("⚠️ 警告：真实矩阵和预测矩阵的行或列可能没有完全对齐！正在自动按真实矩阵的行列顺序对齐...")
    pred_df = pred_df.loc[true_df.index, true_df.columns]

# ==========================================
# 3. 筛选并提取子矩阵
# ==========================================
valid_proteins = [p for p in target_proteins if p in true_df.index]
print(f"🔍 目标名单去重后共有 {len(target_proteins)} 个独立蛋白，在矩阵中成功匹配到 {len(valid_proteins)} 个。")

if len(valid_proteins) == 0:
    print("❌ 错误：没有找到任何匹配的蛋白质！请检查原始表格中的蛋白质名称格式。")
else:
    # 提取只包含这些蛋白质的行
    true_submatrix = true_df.loc[valid_proteins]
    pred_submatrix = pred_df.loc[valid_proteins]
    
    # ==========================================
    # 4. 展平矩阵
    # ==========================================
    y_true_flat = true_submatrix.values.flatten()
    y_pred_flat = pred_submatrix.values.flatten()
    
    num_positive = np.sum(y_true_flat == 1)
    print(f"📊 该子矩阵总计包含 {len(y_true_flat)} 个“蛋白-疾病”对，其中真实存在的关联有 {num_positive} 个。\n")
    
    # ==========================================
    # 5. 计算指标 (AUC + 动态寻找最佳阈值)
    # ==========================================
    try:
        # 1. 计算整体 AUC
        auc_score = roc_auc_score(y_true_flat, y_pred_flat)
        print(f"🏆 【总体性能】糖基化类别预测 AUC 为: {auc_score:.4f}\n")
        
        # 2. 动态寻找最佳阈值
        \
        best_f1 = 0.0
        best_threshold = 0.5
        best_recall = 0.0
        best_precision = 0.0

        for t in np.arange(0.01, 1.0, 0.01):
            y_pred_binary_temp = (y_pred_flat >= t).astype(int)
            current_f1 = f1_score(y_true_flat, y_pred_binary_temp, zero_division=0)
            
            if current_f1 > best_f1:
                best_f1 = current_f1
                best_threshold = t
                best_recall = recall_score(y_true_flat, y_pred_binary_temp, zero_division=0)
                best_precision = precision_score(y_true_flat, y_pred_binary_temp, zero_division=0)

     
        print(f"🎉 Recall: {best_recall:.4f}")
        print(f"🎉 Precision:  {best_precision:.4f}")
        print(f"🎉  F1-score:  {best_f1:.4f}")
        
    except ValueError:
        print("❌ 计算出错：提取出的真实标签中可能全是 0 或全是 1，无法计算指标。")