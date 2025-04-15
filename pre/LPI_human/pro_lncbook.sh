#!/bin/bash

# Set the input file
input_file="../data/lncrna_rbp_LncBook2.0.csv"

# Set the output file
output_file="../data/LncBook_LPI.csv"

# Use gawk for advanced CSV parsing to handle fields enclosed in double quotes that may contain commas
gawk -v FPAT='("([^"]*)")|([^,]+)' -v OFS=',' '{
    # Create a unique identifier based on the first, second, fourth, and fifth columns
    unique_key = $1"\t"$2"\t"$4"\t"$5

    # Remove double quotes from fields
    gsub(/"/, "", $1)
    gsub(/"/, "", $2)
    gsub(/"/, "", $4)
    gsub(/"/, "", $5)

    # Check if fields contain commas, and if so, enclose them in double quotes
    if (index($1, ",") != 0) $1 = "\"" $1 "\""
    if (index($2, ",") != 0) $2 = "\"" $2 "\""
    if (index($4, ",") != 0) $4 = "\"" $4 "\""
    if (index($5, ",") != 0) $5 = "\"" $5 "\""

    # Only keep the first, second, fourth, and fifth columns, and output unique entries
    if (!seen[unique_key]++) {
        print $1, $2, $4, $5 > "'$output_file'"
    }
}' "$input_file"

echo "Filtering and splitting completed. Data has been saved to '$output_file'"