import pandas as pd

# Adjust the path according to the actual file location
# Since you uploaded a file, let's use that
NPInter_inter = pd.read_csv("../../data/all_LPI/human/NPInter_LPI.csv")
pLoF = pd.read_csv('../../protein/pLoF.txt', sep='\t')
valid_protein = pLoF[['gene']]

# Extract the 'Protein' column 
NPinter_protein = NPInter_inter[['tarName', 'tarID', 'reference']]

# Extract invalid proteins in npinter_inter
# Drop duplicates
NPinter_protein = NPinter_protein.drop_duplicates(subset=['tarName', 'tarID'])

# Use a regular expression to filter out rows containing only uppercase letters, digits, and hyphens
# We use `~` to negate the condition, so it selects rows that do NOT match the pattern
NPinter_protein = NPinter_protein[~NPinter_protein['tarName'].isin(valid_protein['gene'])]
NPinter_protein = NPinter_protein.drop_duplicates(['tarName'])
# Save the result to a CSV file
NPinter_protein.to_csv("invalid_npinter_protein.csv", index=False)
