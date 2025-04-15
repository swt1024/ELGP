import pandas as pd

# === 配置文件路径 ===
mapping_path = "right_protein.csv"         # 无表头：旧symbol\t新symbol
protein_path = "../../data/LPPI/mouse/protein.csv"         # symbol\tp+symbol
lppi_path = "../../data/LPPI/mouse/LPPI.csv"               # 两列，可能混合lncRNA和protein

# === 读取 mapping 映射表 ===
df_map = pd.read_csv(mapping_path, header=None, names=["old_symbol", "new_symbol"])
symbol_map = dict(zip(df_map["old_symbol"], df_map["new_symbol"]))

# === 处理 protein 文件 ===
df_protein = pd.read_csv(protein_path)
df_protein["NewSymbol"] = df_protein["protein"].map(symbol_map).fillna(df_protein["protein"])
df_protein["NewProteinID"] = "p" + df_protein["NewSymbol"]

# 保存修正后的 protein 文件
update_protein = df_protein[["NewSymbol", "NewProteinID"]]
update_protein.columns=['protein','protein_ID']
update_protein = update_protein.drop_duplicates()
update_protein.to_csv("../../data/LPPI/mouse/protein_updated.csv", index=False)
print("✅ 已保存：protein_updated.csv")

# === 处理 LPPI 文件 ===
df_lppi = pd.read_csv(lppi_path)

# 构建 protein ID 映射表：p+old_symbol → p+new_symbol
protein_id_map = {
    "p" + old: "p" + new
    for old, new in zip(df_map["old_symbol"], df_map["new_symbol"])
}

# 替换 LPPI 中两个节点中，前缀为 p 的 ID
def replace_protein_id(value):
    return protein_id_map.get(value, value) if value.startswith("p") else value

df_lppi["Node_i"] = df_lppi["Node_i"].apply(replace_protein_id)
df_lppi["Node_j"] = df_lppi["Node_j"].apply(replace_protein_id)

# 保存修正后的 LPPI 文件
df_lppi.to_csv("../../data/LPPI/mouse/LPPI_updated.csv", index=False)
print("✅ 已保存：LPPI_updated.csv（含 LPI + PPI）")
