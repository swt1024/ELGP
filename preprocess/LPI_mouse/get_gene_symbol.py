import pandas as pd

# ====== 配置文件路径 ======
input_synonym_list = "invalid_proteins.txt"          # 每行一个 synonym
mapping_file = "synonym_to_symbol.csv"          # 已生成的映射文件
output_mapped = "mapped_synonyms.csv"           # 输出成功结果
output_unmapped = "unmapped_synonyms.txt"       # 输出未匹配结果

# ====== 读取 synonym 映射表 ======
df_map = pd.read_csv(mapping_file)
df_map["synonym"] = df_map["synonym"].astype(str).str.strip()

# 创建查找表（小写匹配更稳）
syn2symbol = dict(zip(df_map["synonym"].str.lower(), df_map["GeneSymbol"]))

# ====== 读取用户 synonym 列表 ======
with open(input_synonym_list, "r", encoding="utf-8") as f:
    query_synonyms = [line.strip() for line in f if line.strip()]

# ====== 执行映射 ======
mapped = []
unmapped = []

for syn in query_synonyms:
    key = syn.lower()
    if key in syn2symbol:
        mapped.append((syn, syn2symbol[key]))
    else:
        unmapped.append(syn)

# ====== 保存映射成功结果 ======
df_out = pd.DataFrame(mapped, columns=["InputSynonym", "GeneSymbol"])
df_out.to_csv(output_mapped, index=False)

# ====== 保存映射失败列表 ======
with open(output_unmapped, "w", encoding="utf-8") as f:
    for syn in unmapped:
        f.write(syn + "\n")

# ====== 完成提示 ======
print(f"✅ 映射完成！共成功 {len(mapped)} 条，失败 {len(unmapped)} 条")
print(f"📝 成功结果：{output_mapped}")
print(f"⚠️ 未匹配结果：{output_unmapped}")
