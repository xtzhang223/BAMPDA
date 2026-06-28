import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, recall_score, precision_score, f1_score

# ==========================================
# 1. 配置文件路径和目标蛋白质名单
# ==========================================
# 使用 r"" 防止 Windows 路径转义报错
true_matrix_path = r"C:\Users\zxt\Desktop\run_MDHGIBNNR\749-276\five_fold\dataset\Protein_Disease_Associations.xlsx"
pred_matrix_path = r"C:\Users\zxt\Desktop\run_MDHGIBNNR\749-276\paint-roc-data1\Predict_Score_Matrix3_3.xlsx"

# 乙酰化相关的蛋白质名单
target_proteins = [
    "got2", "esr1", "tp53", "eno1", "ar", "gfap", "tal1", 
    "calr", "cdkn1a", "acly", "hist1h4a", "hist1h3a", 
    "hist1h3b", "h3f3a", "insm1", "pou5f1", "rela", 
    "khdrbs1", "rbl2", "klf6", "bub1b"
]

# 将目标蛋白质名字全部转换为大写，确保匹配万无一失
target_proteins = [p.upper() for p in target_proteins]

print("正在加载数据，请稍候...")

# ==========================================
# 2. 读取 Excel 文件 (index_col=0 将第一列设为行索引)
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
# 检查我们的名单里有多少个蛋白质真正在矩阵中
valid_proteins = [p for p in target_proteins if p in true_df.index]
print(f"🔍 目标名单共有 {len(target_proteins)} 个蛋白，在矩阵中成功匹配到 {len(valid_proteins)} 个。")

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
    
    # 统计一下这个子矩阵里一共有多少个正样本 (即真实存在的关联)
    num_positive = np.sum(y_true_flat == 1)
    print(f"📊 该子矩阵总计包含 {len(y_true_flat)} 个“蛋白-疾病”对，其中真实存在的关联有 {num_positive} 个。\n")
    
    # ==========================================
    # 5. 计算指标 (AUC + 动态寻找最佳阈值)
    # ==========================================
    try:
        # 1. 计算整体 AUC (AUC不受阈值影响，最能代表全局排序能力)
        auc_score = roc_auc_score(y_true_flat, y_pred_flat)
        print(f"🏆 【总体性能】乙酰化类别预测 AUC 为: {auc_score:.4f}\n")
        
        # 2. 动态寻找最佳阈值 (使 F1-score 最大的阈值)
       
        best_f1 = 0.0
        best_threshold = 0.5
        best_recall = 0.0
        best_precision = 0.0

        # 遍历从 0.01 到 0.99 的所有可能阈值，步长为 0.01
        for t in np.arange(0.01, 1.0, 0.01):
            # 根据当前阈值 t 将预测分数转化为 0 和 1
            y_pred_binary_temp = (y_pred_flat >= t).astype(int)
            
            # 计算当前的 F1-score
            current_f1 = f1_score(y_true_flat, y_pred_binary_temp, zero_division=0)
            
            # 如果发现更高的 F1，就更新我们的记录
            if current_f1 > best_f1:
                best_f1 = current_f1
                best_threshold = t
                best_recall = recall_score(y_true_flat, y_pred_binary_temp, zero_division=0)
                best_precision = precision_score(y_true_flat, y_pred_binary_temp, zero_division=0)

        
        print(f"🎉  Recall:    {best_recall:.4f}")
        print(f"🎉  Precision: {best_precision:.4f}")
        print(f"🎉  F1-score:  {best_f1:.4f}")
        
    except ValueError:
        print("❌ 计算出错：提取出的真实标签中可能全是 0 或全是 1，无法计算指标。")