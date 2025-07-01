  #!/bin/sh

## Cross reference splicing events with Unipro database (*hg38.col.bed)
#!/bin/sh

# Define the input and output directories
input_dir="input/"
output_dir="results/"

# Define the bed files
bed_files=("UP000005640_9606_mod_res.bed" "UP000005640_9606_disulfid.bed" "UP000005640_9606_signal.bed" "UP000005640_9606_domain.bed")

# Define the splicing event files
splicing_event_files=("splicing_events.SE.total.pos.bed" "splicing_events.SE.total.neg.bed")

# Loop through each splicing event file and bed file combination
for splicing_file in "${splicing_event_files[@]}"; do
    for bed_file in "${bed_files[@]}"; do
        output_file="${output_dir}${splicing_file%.*}.intersect_${bed_file%.*}.wo.txt"
        bedtools intersect -wo -a "${output_dir}${splicing_file}" -b "${input_dir}${bed_file}" | sort -u > "$output_file"
    done
done
