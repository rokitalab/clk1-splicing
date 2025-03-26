#!/bin/bash

cd results
# Output file
echo -e "Filename\tGeneset\tLogFC\tadj.P.Val" > kegg_spliceosome_summary.tsv

# Loop through all *pathway_stranded.tsv files in the current directory
for file in *pathway_stranded.tsv; do
  # Skip if no matching files
  [ -e "$file" ] || continue

  # Extract the header to get column numbers
  header=$(head -n 1 "$file")
  geneset_col=$(echo "$header" | tr '\t' '\n' | grep -n '^Geneset$' | cut -d: -f1)
  logfc_col=$(echo "$header" | tr '\t' '\n' | grep -n '^logFC$' | cut -d: -f1)
  adjp_col=$(echo "$header" | tr '\t' '\n' | grep -n '^adj.P.Val$' | cut -d: -f1)

  # If both columns are found
  if [[ -n "$geneset_col" && -n "$logfc_col" && -n "$adjp_col" ]]; then
    awk -v file="$file" -v gcol="$geneset_col" -v lcol="$logfc_col" -v apcol="$adjp_col" 'BEGIN{FS=OFS="\t"}
    NR > 1 && $gcol == "KEGG_SPLICEOSOME" {
      print file, $gcol, $lcol, $apcol
    }' "$file" >> kegg_spliceosome_summary.tsv
  fi
done
