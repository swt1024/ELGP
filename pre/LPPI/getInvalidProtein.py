import pandas as pd

protein = pd.read_csv("../../data/PPI/human/protein_in_ppi.csv")
pLoF = pd.read_csv('../../protein/human/pLoF.txt', sep='\t')
valid_protein = pLoF[['gene']]

# Use a regular expression to filter out rows containing only uppercase letters, digits, and hyphens
# We use `~` to negate the condition, so it selects rows that do NOT match the pattern
invalid_protein = protein[~protein['protein'].isin(valid_protein['gene'])]
# Save the result to a CSV file
invalid_protein.to_csv("invalid_protein.csv", index=False)
