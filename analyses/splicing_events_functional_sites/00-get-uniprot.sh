#!/bin/bash

# Target directory
target_dir="input"

# List of URLs
urls=(
  "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/genome_annotation_tracks/UP000005640_9606_beds/UP000005640_9606_signal.bed"
  "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/genome_annotation_tracks/UP000005640_9606_beds/UP000005640_9606_disulfid.bed"
  "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/genome_annotation_tracks/UP000005640_9606_beds/UP000005640_9606_domain.bed"
  "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/genome_annotation_tracks/UP000005640_9606_beds/UP000005640_9606_mod_res.bed"
)

# Loop over URLs
for url in "${urls[@]}"; do
  filename=$(basename "$url")
  filepath="$target_dir/$filename"
  if [ -f "$filepath" ]; then
    echo "$filename already exists in $target_dir. Skipping download."
  else
    echo "Downloading $filename to $target_dir..."
    curl -o "$filepath" "$url"
  fi
done

