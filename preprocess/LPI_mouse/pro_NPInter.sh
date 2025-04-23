#!/bin/bash

# Set the input file
input_file="../../data/NPInter_LPI/lncRNA_interaction.txt"

# Output files
mouse_output="../../data/NPInter_LPI/mouse/NPInter_LPI.csv"

# Process the input file and filter records based on specific conditions
awk -F'\t' '($11 == "Mus musculus") && $4 == "lncRNA" && $7 == "protein" && $14 == "binding" {

    print (($1 ~ /^".*"$/ ? $1 : "\"" $1 "\"")) "," \
            (($2 ~ /^".*"$/ ? $2 : "\"" $2 "\"")) "," \
            (($3 ~ /^".*"$/ ? $3 : "\"" $3 "\"")) "," \
            (($5 ~ /^".*"$/ ? $5 : "\"" $5 "\"")) "," \
            (($6 ~ /^".*"$/ ? $6 : "\"" $6 "\"")) "," \
            (($10 ~ /^".*"$/ ? $10 : "\"" $10 "\"")) "," \
            (($12 ~ /^".*"$/ ? $12 : "\"" $12 "\"")) > "'$mouse_output'"

}' "$input_file"

# Add column headers to the output files
echo -e "inter_id,symbol,gene_id,tarName,tarID,reference,tissueOrCell" | cat - "$mouse_output" > temp && mv temp "$mouse_output"

# Display completion message
echo "Filtering and splitting completed. Data has been saved to '$mouse_output'."
