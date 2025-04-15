import pandas as pd

# Load the LPI file that needs the corrections
inter_path = '../../data/NPInter_LPI/human/NPInter_LPI.csv'  # Replace with the actual file path
inter = pd.read_csv(inter_path)

noncodev5_trans = pd.read_csv('../../data/reference_lncRNA/human/transcript/NONCODEv5_human_hg38_lncRNA_trans.txt', sep='\t', header=None, names=['gene_id', 'transcript_id'])
noncodev6_trans = pd.read_csv('../../data/reference_lncRNA/human/transcript/NONCODEv6_human_hg38_lncRNA_trans.txt', sep='\t', header=None, names=['gene_id', 'transcript_id'])

noncodev5_trans['transcript_id'] = noncodev5_trans['transcript_id'].str.split('.').str[0]
noncodev6_trans['transcript_id'] = noncodev6_trans['transcript_id'].str.split('.').str[0]

lnc_trans_dict_v5 = dict(zip(noncodev5_trans['transcript_id'], noncodev5_trans['gene_id']))
lnc_trans_dict_v6 = dict(zip(noncodev6_trans['transcript_id'], noncodev6_trans['gene_id']))

# Replace transcript IDs with corresponding gene IDs
inter['gene_id'] = inter['gene_id'].map(lnc_trans_dict_v5).fillna(inter['gene_id'])
inter['gene_id'] = inter['gene_id'].map(lnc_trans_dict_v6).fillna(inter['gene_id'])
print("Replace transcript IDs with corresponding gene IDs")

# Load the CSV file that contains the mapping of incorrect to correct proteins
protein_mapping_path = 'right_protein.csv'  # Replace with the actual file path
protein_mapping = pd.read_csv(protein_mapping_path, header=None, names=['wrong', 'correct'])
invalid_protein = pd.read_csv('exact_invalid_protein.csv', header=None, names=['protein'])

# Create a dictionary mapping from incorrect names to correct names
protein_dict = dict(zip(protein_mapping['wrong'], protein_mapping['correct']))

# Replace incorrect proteins in the 'tarName' column using the dictionary
inter['tarName'] = inter['tarName'].replace(protein_dict)
inter = inter[~inter['tarName'].isin(invalid_protein['protein'])]
print("Correct invalid protein names")

# Load the file that contains right cell name
mapping_file = "correct_cell_tissue.csv"  # The mapping file containing Original → Standardized names
output_csv = "../../data/NPInter_LPI/human/correct_NPInter_LPI.csv"  # The output CSV file with updated values

# Read the mapping file (standardized names)
mapping_df = pd.read_csv(mapping_file)
mapping_dict = dict(zip(mapping_df['Original Name'], mapping_df['Standardized Name']))

# Replace values in the 'tissueOrCell' column based on the mapping dictionary
inter['tissueOrCell'] = inter['tissueOrCell'].replace(mapping_dict)
print("Fix incorrect cell line and tissue names")

# Save the modified CSV file
inter.to_csv(output_csv, index=False)
print(f"Output saved to {output_csv}")
