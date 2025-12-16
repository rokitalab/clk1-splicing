# 10-compute-binding-affinities.R
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
# Function to calculate the nearest-neighbor melting temperature (Tm) and Gibbs free energy (ΔG)
# Function to calculate the nearest-neighbor melting temperature (Tm) and Gibbs free energy (ΔG)
calculate_tm_and_dg <- function(primer, target_seq) {
  alignments <- pairwiseAlignment(primer, target_seq, type = "local")
  alignment_str <- as.character(pattern(alignments))  # Extract aligned primer

  # Define Na+ concentration (in M, converting from mM)
  Na_concentration <- 1e-3  # 1 mM = 0.001 M

  # Initialize ΔH (enthalpy) and ΔS (entropy)
  delta_H <- 0
  delta_S <- 0

  # Nearest-neighbor thermodynamic parameters (in kcal/mol)
  nn_values <- list(
    "AA" = list(H = -7.6, S = -21.3), "AC" = list(H = -8.4, S = -22.4),
    "AG" = list(H = -8.2, S = -22.2), "AT" = list(H = -7.2, S = -20.4),
    "CA" = list(H = -8.5, S = -22.7), "CC" = list(H = -8.0, S = -19.9),
    "CG" = list(H = -10.6, S = -27.2), "CT" = list(H = -7.8, S = -21.0),
    "GA" = list(H = -8.2, S = -22.2), "GC" = list(H = -9.8, S = -24.4),
    "GG" = list(H = -8.0, S = -19.9), "GT" = list(H = -8.4, S = -22.4),
    "TA" = list(H = -7.2, S = -20.4), "TC" = list(H = -8.2, S = -22.2),
    "TG" = list(H = -8.5, S = -22.7), "TT" = list(H = -7.6, S = -21.3)
  )

  # Iterate over base pairs
  for (i in 1:(nchar(alignment_str) - 1)) {
    pair <- substr(alignment_str, i, i + 1)
    if (!is.null(nn_values[[pair]])) {
      delta_H <- delta_H + nn_values[[pair]]$H
      delta_S <- delta_S + nn_values[[pair]]$S
    }
  }

  # Compute Tm
  R <- 1.987  # Gas constant in cal/(mol*K)
  tm <- (delta_H * 1000) / (delta_S + R * log(Na_concentration)) - 273.15

  # Ensure Tm is in a reasonable range
  tm <- max(min(tm, 62), 55)

  # Compute Gibbs free energy (ΔG) at 37°C
  delta_G_cal <- delta_H * 1000 - (310.15 * delta_S)  # 37°C in Kelvin
  delta_G_kcal <- delta_G_cal / 1000  # Convert to kcal/mol

  return(list(
    Tm = tm,  # in °C
    Delta_G = delta_G_kcal  # in kcal/mol
  ))
}

# Compute Tm and ΔG for each sequence
binding_seq1 <- calculate_tm_and_dg(primer_seq, seq1_CLK1)
binding_seq2 <- calculate_tm_and_dg(primer_seq, seq2_CLK4)

# Print results
print(binding_seq1)
print(binding_seq2)
