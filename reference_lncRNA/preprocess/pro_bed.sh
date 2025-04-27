# Retain rows categorized as 'gene' and the required columns
awk '$4 ~ /^NONHSAG/ {sub(/\..*/, "", $4); print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ../human/bed/NONCODEv6_hg38.lncAndGene.bed > ../human/bed/NONCODEv6_hg38.lncRNAGene.bed
awk '$4 ~ /^NONHSAG/ {sub(/\..*/, "", $4); print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ../human/bed/NONCODEv5_hg38.lncAndGene.bed > ../human/bed/NONCODEv5_hg38.lncRNAGene.bed
awk '$4 ~ /^NONMMUG/ {sub(/\..*/, "", $4); print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ../mouse/bed/NONCODEv5_mm10.lncAndGene.bed > ../mouse/bed/NONCODEv5_mm10.lncRNAGene.bed
awk '$4 ~ /^NONMMUG/ {sub(/\..*/, "", $4); print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ../mouse/bed/NONCODEv6_mm10.lncAndGene.bed > ../mouse/bed/NONCODEv6_mm10.lncRNAGene.bed

# Retain rows categorized as 'gene' and the required columns from lncRNA_LncBookv2.0_GRCh38.gtf file, and convert to bed format
awk -F'\t' '$3 == "gene" {split($9, a, "\""); print $1 "\t" $4-1 "\t" $5 "\t" a[2] "\t" $6 "\t" $7}' ../human/bed/lncRNA_LncBookv2.0_GRCh38.gtf > ../human/bed/lncRNA_gene_LncBookv2.0_GRCh38.bed
