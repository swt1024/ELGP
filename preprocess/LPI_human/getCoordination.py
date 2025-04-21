import pandas as pd
import os
import re

# Set the directory containing gencode data
gencode_dir = "../../data/reference_lncRNA/human/bed/gencode/"

# Read lncRNA ID and Symbol information
inter_lncbook = pd.read_csv('../../data/LncBook_LPI/LncBook_LPI.csv')
inter_npinter = pd.read_csv('../../data/NPInter_LPI/human/correct_NPInter_LPI.csv')

lncRNA_lncbook = inter_lncbook[['Gene ID', 'Symbol']].drop_duplicates()
lncRNA_lncbook.columns = ['gene_id', 'symbol']
lncRNA_npinter = inter_npinter[['gene_id', 'symbol']].drop_duplicates()

# Read NONCODE and LncBook BED files
lncbook_bed = pd.read_csv('../../data/reference_lncRNA/human/bed/lncRNA_gene_LncBookv2.0_GRCh38.bed', sep='\t', header=None, names=['chr', 'start', 'end', 'gene_id', 'score', 'strand'])
noncodev5_bed = pd.read_csv('../../data/reference_lncRNA/human/bed/NONCODEv5_hg38.lncRNAGene.bed', sep='\t', header=None, names=['chr', 'start', 'end', 'gene_id', 'score', 'strand'])
noncodev6_bed = pd.read_csv('../../data/reference_lncRNA/human/bed/NONCODEv6_hg38.lncRNAGene.bed', sep='\t', header=None, names=['chr', 'start', 'end', 'gene_id', 'score', 'strand'])

noncodev5_bed['gene_id'] = noncodev5_bed['gene_id'].str.split('.').str[0]
noncodev6_bed['gene_id'] = noncodev6_bed['gene_id'].str.split('.').str[0]

# Get genomic position by lncbook_id
pos_lnc_lncbook = pd.merge(lncRNA_lncbook, lncbook_bed[['chr', 'start', 'end', 'strand', 'gene_id']], on='gene_id', how='inner')
pos_lnc_lncbook['identifier'] = pos_lnc_lncbook['gene_id']

# Initialize remaining lncRNA list
remained_lncRNA_npinter = lncRNA_npinter[['gene_id', 'symbol']].copy()
remained_lncRNA_npinter['identifier'] = remained_lncRNA_npinter['gene_id']
results = []

# Get genomic position by noncode_id (noncodev6)
pos_lnc_noncodev6 = pd.merge(
    remained_lncRNA_npinter,
    noncodev6_bed[['chr', 'start', 'end', 'strand', 'gene_id']],
    on='gene_id',
    how='inner'
)
results.append(pos_lnc_noncodev6)
remained_lncRNA_npinter = remained_lncRNA_npinter[~remained_lncRNA_npinter['gene_id'].isin(pos_lnc_noncodev6['gene_id'])]

# Get genomic position by noncode_id (noncodev5)
pos_lnc_noncodev5 = pd.merge(
    remained_lncRNA_npinter,
    noncodev5_bed[['chr', 'start', 'end', 'strand', 'gene_id']],
    on='gene_id',
    how='inner'
)
results.append(pos_lnc_noncodev5)
remained_lncRNA_npinter = remained_lncRNA_npinter[~remained_lncRNA_npinter['gene_id'].isin(pos_lnc_noncodev5['gene_id'])]

# **STEP 1: 提取 BED 文件版本号**
def extract_version(filename):
    #match = re.search(r'GRCh38\.(\d+)', filename)  # ensembl
    match = re.search(r"v(\d+)", filename) # gencode
    return int(match.group(1)) if match else -1  # 

# **STEP 2: 获取所有 BED 文件，并按版本号排序**
bed_files = [f for f in os.listdir(gencode_dir) if f.endswith(".bed")]
bed_files_sorted = sorted(bed_files, key=extract_version, reverse=True)  # 按版本号降序排列

# **STEP 3: 按排序后的顺序遍历 BED 文件**
for bed_file in bed_files_sorted:
    bed_path = os.path.join(gencode_dir, bed_file)
    print(f"Processing {bed_file}...")

    # 读取 gencode BED 文件
    gencode_bed = pd.read_csv(bed_path, sep='\t', header=None, names=['chr', 'start', 'end', 'symbol', 'gene_id', 'strand'])

    # Match by gene_id
    pos_lnc_gencode = pd.merge(remained_lncRNA_npinter, gencode_bed[['chr', 'start', 'end', 'strand', 'gene_id']], on='gene_id', how='inner')
    results.append(pos_lnc_gencode)
    remained_lncRNA_npinter = remained_lncRNA_npinter[~remained_lncRNA_npinter['gene_id'].isin(pos_lnc_gencode['gene_id'])]

    # Match by symbol
    pos_lnc_symbol = pd.merge(remained_lncRNA_npinter, gencode_bed[['chr', 'start', 'end', 'strand', 'symbol']], on='symbol', how='inner')
    remained_lncRNA_npinter = remained_lncRNA_npinter[~remained_lncRNA_npinter['symbol'].isin(pos_lnc_symbol['symbol'])]
    pos_lnc_symbol['identifier'] = pos_lnc_symbol['symbol']
    results.append(pos_lnc_symbol)

# Combine all results
pos_lnc_npinter = pd.concat(results, ignore_index=True).drop_duplicates(subset=['identifier'])

# Save remaining lncRNAs without genomic positions
remained_lncRNA_npinter.drop_duplicates().to_csv('gencode/lnc_no_pos.csv', index=False)


# Generate BED files
pos_lnc_npinter[['chr', 'start', 'end', 'gene_id', 'symbol', 'strand', 'identifier']].to_csv('gencode/npinter_lnc_0based.bed', index=False, sep='\t', header=None)
pos_lnc_lncbook[['chr', 'start', 'end', 'gene_id', 'symbol', 'strand', 'identifier']].to_csv('lncbook_lnc_0based.bed', index=False, sep='\t', header=None)

print("Processing complete. Results saved.")
