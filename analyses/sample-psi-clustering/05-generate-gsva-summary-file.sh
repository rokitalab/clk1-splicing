#!/bin/bash

cd results

# Output file
output_file="all_gsva_de_results_stranded.tsv"
first=1

# Loop through all *pathway_stranded.tsv files in the current directory
for file in *pathway_stranded.tsv; do
  [ -e "$file" ] || continue

  if [[ $first -eq 1 ]]; then
    # Add header from first file, with "Filename" column
    echo -e "Filename\t$(head -n 1 "$file")" > "$output_file"
    first=0
  fi

  # Append data with filename as first column
  tail -n +2 "$file" | awk -v fname="$file" 'BEGIN{OFS="\t"} {print fname, $0}' >> "$output_file"
done