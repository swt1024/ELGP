import pandas as pd

# Load example DataFrames
df1 = pd.read_csv("HinSAGE/mouse/lncRNA_embeddings_brain.csv",header=None)
df2 = pd.read_csv("HGVAE/mouse/lncRNA_embeddings_brain.csv",header=None)

# Set the 'lncRNA_ID' column as index
df1.set_index(0, inplace=True)
df2.set_index(0, inplace=True)

# Concatenate based on the index (original 'lncRNA_ID' column)
result = pd.concat([df1, df2], axis=1)

# Save the merged result to CSV without the index
result.to_csv('mouse_merged_brain.csv',header=None)

