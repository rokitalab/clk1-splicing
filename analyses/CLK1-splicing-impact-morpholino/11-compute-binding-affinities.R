# 11-compute-binding-affinities.R
# This script computes the binding affinity of a primer to two target sequences.
# The script performs the following:
# - Calculates a basic alignment score based on mismatches
# - Computes the melting temperature (Tm) and Gibbs free energy (ΔG) for primer-target binding
#
# Dependencies:
# - Biostrings: for sequence handling and pairwise alignment

library(Biostrings)

# Define the primer sequence
primer_seq <- DNAString("ACTCTTCTGGAAACGTCAAGTGGGC")

# Define the target sequences
seq1_CLK1 <- DNAString(paste0(
  "ATCTTGCACTTAGTACGTCTCCACTCTGACAGATCAGGTGACCCTCCTCATCATCCTCTACACTCCTGGTTCTTTTCCTTCGGTGACTCTTCTGGAAACGTCAAGTGGGC"
))

seq2_CLK4 <- DNAString(paste0(
  "ATCTTGCCTTAGACGTCTCCACTTGACAGATCAGGTGACCCTCCTCATCATCCTCTAACTCCTGGTCTTTTCCTTCGGTGCTCTTCTGGAAGAGG"
))

# Function to compute binding affinity (basic score based on mismatches)
compute_binding_affinity <- function(primer, target_seq) {
  alignments <- pairwiseAlignment(primer, target_seq, type = "local")
  mismatch_count <- sum(mismatchTable(alignments)$mismatch == TRUE)
  alignment_score <- score(alignments)  # Higher score means better binding
  
  return(list(
    alignment_score = alignment_score,
    mismatches = mismatch_count
  ))
}

# Compute binding affinity for each sequence
binding_seq1 <- compute_binding_affinity(primer_seq, seq1_CLK1)
binding_seq2 <- compute_binding_affinity(primer_seq, seq2_CLK4)

# Print results for binding affinity
print(binding_seq1)
print(binding_seq2)

## TM and Gibbs Free Energy Calculation
# Define nearest neighbor values (simplified)
nn_values <- list(
  AA = list(H = -7.6, S = -21.3), AC = list(H = -8.4, S = -22.4), AG = list(H = -8.4, S = -22.4), AT = list(H = -7.8, S = -21.0),
  CA = list(H = -8.4, S = -22.4), CC = list(H = -8.0, S = -22.0), CG = list(H = -10.6, S = -27.2), CT = list(H = -7.8, S = -21.0),
  GA = list(H = -8.4, S = -22.4), GC = list(H = -10.6, S = -27.2), GG = list(H = -9.0, S = -24.0), GT = list(H = -8.4, S = -22.4),
  TA = list(H = -7.6, S = -21.3), TC = list(H = -7.8, S = -21.0), TG = list(H = -8.4, S = -22.4), TT = list(H = -7.6, S = -21.3)
)

# Function to calculate the nearest-neighbor melting temperature (Tm) and Gibbs free energy (ΔG)
calculate_tm_and_dg <- function(primer, target_seq) {
  alignments <- pairwiseAlignment(primer, target_seq, type = "local")
  alignment_str <- as.character(alignments)
  
  # Define Na+ concentration (1 M)
  Na_concentration <- 1
  
  # Initialize ΔH (enthalpy) and ΔS (entropy)
  delta_H <- 0
  delta_S <- 0
  
  # Calculate the nearest-neighbor values
  for (i in 1:(nchar(alignment_str) - 1)) {
    pair <- substr(alignment_str, i, i + 1)
    if (pair %in% names(nn_values)) {
      delta_H <- delta_H + nn_values[[pair]]$H
      delta_S <- delta_S + nn_values[[pair]]$S
    }
  }
  
  # Calculate melting temperature (Tm) using the nearest-neighbor model
  tm <- (delta_H / delta_S) - 273.15 + 16.6 * log10(Na_concentration)
  
  # Calculate Gibbs free energy (ΔG)
  delta_G <- delta_H - tm * delta_S
  
  return(list(
    Tm = tm,
    Delta_G = delta_G
  ))
}

# Compute Tm and ΔG for each sequence
binding_seq1 <- calculate_tm_and_dg(primer_seq, seq1_CLK1)
binding_seq2 <- calculate_tm_and_dg(primer_seq, seq2_CLK4)

# Print results for Tm and ΔG
print(binding_seq1)
print(binding_seq2)
