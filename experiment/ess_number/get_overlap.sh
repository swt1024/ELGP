bedtools intersect \
  -a mouse_essential_genes_union.bed \
  -b mouse_essential_genes_union.bed \
  -s \
  -wo | \
awk '{
    len_a = $3 - $2;
    len_b = $9 - $8;
    shorter = (len_a < len_b) ? len_a : len_b;
    if ($NF / shorter >= 0.5 && $4 != $10)
        print $4, $10;
}' > mouse_overlapping_genes.txt./