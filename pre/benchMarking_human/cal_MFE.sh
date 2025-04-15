#!/bin/bash

Check the number of input arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_fasta_file>"
    exit 1
fi

input_fasta="$1"
output_file="trans_MFE.csv"

Check if the input FASTA file exists
if [ ! -f "$input_fasta" ]; then
    echo "Error: Input FASTA file does not exist."
    exit 1
fi

Split the FASTA file so each transcript is stored in a separate file (used for parallel RNAfold computation)
mkdir -p temp_fasta
awk '/^>/ {
    if (seq) print seq > filename;  
    gsub(/\..*/, "", $1);           
    id=substr($1,2);                
    filename="temp_fasta/"id".fasta";  
    print $0 > filename;           
    seq="";                         
    next;
} 
{ 
    seq = seq ? seq "\n" $0 : $0;   
} 
END { 
    if (seq) print seq > filename   
}' "$input_fasta"

echo "transcript_id,MFE" > trans_MFE.csv

# Perform parallel RNAfold calculations of MFE
find "$fasta_dir" -name "*.fasta" | parallel --jobs 16 --bar '
    id=$(basename {} .fasta);
    mfe=$(RNAfold --noPS < {} | awk "/)/ {print \$NF}" | tr -d "()")
    echo "$id,$mfe"
' > trans_MFE.csv

# Clean up temporary files
rm -rf temp_fasta 

echo "Final MFE calculations have been saved to $output_file"
