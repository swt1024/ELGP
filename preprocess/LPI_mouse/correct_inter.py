import pandas as pd

# Load the LPI file that needs the corrections
inter_path = '../../data/NPInter_LPI/mouse/NPInter_LPI.csv'  # Replace with the actual file path
inter = pd.read_csv(inter_path)

# Delete inter with lncRNA that without position
lncRNA_pos = pd.read_csv("../../data/LPI/mouse/lncRNA.csv")

noncodev5_trans = pd.read_csv('../../data/reference_lncRNA/mouse/transcript/NONCODEv5_mouse_mm10_lncRNA_trans.txt', sep='\t', header=None, names=['gene_id', 'transcript_id'])

noncodev5_trans['transcript_id'] = noncodev5_trans['transcript_id'].str.split('.').str[0]

lnc_trans_dict_v5 = dict(zip(noncodev5_trans['transcript_id'], noncodev5_trans['gene_id']))

# Replace transcript IDs with corresponding gene IDs
inter['gene_id'] = inter['gene_id'].map(lnc_trans_dict_v5).fillna(inter['gene_id'])
print("Replace transcript IDs with corresponding gene IDs")

inter.loc[(inter['gene_id'] == '-') & (inter['symbol'].str.startswith('ENSMUSG')), 'gene_id'] = inter['symbol']

inter.loc[(inter['tarName'] == '-'), 'tarName' ] = inter['tarID']

# Save the modified CSV file
inter = inter[['gene_id', 'symbol', 'tarName']]
inter.to_csv("../../data/NPInter_LPI/mouse/correct_NPInter_LPI.csv", index=False)