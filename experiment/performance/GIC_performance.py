import pandas as pd
import numpy as np
from sklearn.metrics import (
    confusion_matrix, f1_score, matthews_corrcoef, accuracy_score,
    precision_score, recall_score
)

# === Step 1: 读取数据 ===

# 替换为你自己的文件路径
gic_file = '../../results/mouse/mouse_GIC_score.csv'
essential_file = '../../data/benchMarking/mouse/ess_lpi.csv'
nonessential_file = '../../data/benchMarking/mouse/noness_lpi.csv'

# 读取数据
gic_df = pd.read_csv(gic_file)
ess_ids = pd.read_csv(essential_file, header=None)[0].tolist()
noness_ids = pd.read_csv(nonessential_file, header=None)[0].tolist()

# 只保留在必需和非必需列表中的样本
gic_df = gic_df[gic_df['lncRNA_ID'].isin(ess_ids + noness_ids)]

# 加标签：1 = 必需，0 = 非必需
gic_df['Label'] = gic_df['lncRNA_ID'].apply(lambda x: 1 if x in ess_ids else (0 if x in noness_ids else np.nan))
gic_df.dropna(inplace=True)

# 准备数据
y_true = gic_df['Label'].astype(int).values
y_scores = gic_df['Score'].values

# === Step 2: 指标计算函数 ===
def compute_metrics(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    sen = tp / (tp + fn) if (tp + fn) > 0 else 0
    spe = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
    acc = accuracy_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)
    mcc = matthews_corrcoef(y_true, y_pred)
    return {
        'Sensitivity': sen,
        'Specificity': spe,
        'PPV': ppv,
        'Accuracy': acc,
        'F1': f1,
        'MCC': mcc
    }

# === Step 3: 使用固定阈值 0.5 计算指标 ===
threshold = 0.5
y_pred = (y_scores >= threshold).astype(int)
metrics = compute_metrics(y_true, y_pred)

# === Step 4: 输出结果 ===
print("📌 使用固定阈值 0.5 的结果")
print(f"Threshold: {threshold}")
print(metrics)
