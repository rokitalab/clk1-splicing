#!/usr/bin/perl

use List::Util qw(max);

################################################################################
# create_matrix_of_PSI_SE_gene.pl
# create tsv file of PSI for each sample, to be used for consensus clustering
# written by Ammar Naqvi
#
# usage: ./create_matrix_of_PSI_removeDups.pl <histology file> <rMATs_file_paths.txt>
#                                             <primary tumor file> <primary plus file> <outfile>
################################################################################

## get args-- input files and output file
my ($histology, $rmats_tsv, $primary_tumor_dat,$out_file) = @ARGV[0, 1, 2, 3];

## data structures
my (@broad_hist, @bs_id, @splicing_events);
my (%histology_ids, %inc_levels, %bs_id_hist, %hist_check, %hist_count);
my %primary_initial_sample_list;

## store primary tumor samples
open(FIL, $primary_tumor_dat) || die("Cannot Open File");
while (<FIL>) {
    chomp;
    my @header = split "\t";
    my $bs_id = $header[1];
    my $cohort = $header[2];
    my $exp_strategy = $header[4];
    my $tumor_descr = $header[5];

    next unless ($cohort =~ /PBTA/);

    $primary_initial_sample_list{$bs_id} = $bs_id;
}
close(FIL);

# annotate histology file #
open(FIL, $histology) || die("Cannot Open File");

# Assuming the first line of the file contains column headers
my $header = <FIL>;
chomp $header;
my @column_names = split "\t", $header;

# Create a hash to map column names to indices
my %column_index;
@column_index{@column_names} = (0..$#column_names);

while (<FIL>) {
    chomp;
    my @cols = split "\t";
    my $hist = $cols[$column_index{'plot_group'}];
    my $bs_id = $cols[$column_index{'Kids_First_Biospecimen_ID'}];
    my $CNS_region = $cols[$column_index{'CNS_region'}];

    next unless ($primary_initial_sample_list{$bs_id});

    ## make an array and store histology information and BS IDs
    push @broad_hist, $hist;
    push @bs_ids, $bs_id;

    $bs_id_hist{$bs_id} = $hist;

    ## store total number of histologies
    $hist_count{$hist}++;
    push @{$histology_ids{$hist}}, $bs_id;

    $cns_regions{$bs_id} = $CNS_region;

    ## histology counter for downstream analysis
    $hist_count{$hist}++;
}

close(FIL);

my %inc_levels_gene;
my @genes;
my %hist_check_gene;

## process rMATS output (may take awhile) from merged input file
print "\tprocessing rMATs results...\n";
open(FIL, "gunzip -c $rmats_tsv |") || die ("can’t open $rmats_tsv");
while (<FIL>) {
    chomp;

    ## only look at exon splicing events in primary samples
    next unless ($_ =~ /^SE/);

    my @cols = split "\t";
    my $bs_id = $cols[1];
    my $ctrl = $cols[2];

    next unless ($primary_initial_sample_list{$bs_id});

    ## get histology for BS ID#
    my $hist_of_sample = $bs_id_hist{$bs_id};

    ## get gene name
    my $gene = $cols[4];
    $gene =~ s/\"//g; # remove quotation marks

    ## retrieve exon coordinates
    my $chr = $cols[5];
    my $str = $cols[6];
    my $exonStart = $cols[7]+1; ## its 0-based so add 1 to overlap with uni-prot coord
    my $exonEnd = $cols[8];
    my $upstreamES = $cols[15];
    my $upstreamEE = $cols[16];
    my $downstreamES = $cols[17];
    my $downstreamEE = $cols[18];

    ## retrieve inclusion level and junction count info
    my $inc_level_tumor = $cols[33];
    my $tumor_IJC = $cols[25];
    my $tumor_SJC = $cols[26];

    # get lengths of exons
    my $inc_len = $cols[29];
    my $skip_len = $cols[30];

    ## only look at strong changes
    next unless (($inc_level_tumor >= .10) && ($inc_level_tumor <= .90));
    next unless (($tumor_IJC + $tumor_SJC) >= 100);

    # annotate and re-name splice event IDs
    my $splice_id = $gene."_".$exonStart."-".$exonEnd."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE;

    ## store inclusion levels of splicing event and BS ID#
    $inc_levels{$splice_id}{$bs_id} = $inc_level_tumor;

    ## store array of inclusion levels per gene to access later
    push @{$inc_levels_gene{$gene}{$bs_id}}, $inc_level_tumor;
    push @genes, $gene;

    ## get histology for BS ID#
    my $hist_of_sample = $bs_id_hist{$bs_id};

    ## store number of events per histology
    $hist_check{$splice_id}{$hist_of_sample}++;

    push @splicing_events, $splice_id;

    ## store gene / histology
    $hist_check_gene{$gene}++;
}
close(FIL);

my @broad_hist_uniq = do { my %seen; grep { !$seen{$_}++ } @broad_hist };
my @bs_ids_uniq = do { my %seen; grep { !$seen{$_}++ } @bs_ids };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };
my @genes_uniq = do { my %seen; grep { !$seen{$_}++ } @genes };

# go through each splicing event
# order by disease library
# print splice id, inc level
open(OUT, ">", $out_file) || die("Cannot Open File");

## save and print header info for output file
print OUT "Splice_ID";
foreach my $hist (sort @broad_hist_uniq) {
    my $i = 1;
    foreach my $sample (@{$histology_ids{$hist}}) {
        print OUT "\t";
        print OUT $sample;
        $i++;
    }
}
print OUT "\n";

my $num_events_per_gene = 0;
## print each event and threshold value to output file
foreach my $event (@genes_uniq) {
  $num_events_per_gene = $hist_check_gene{$event};
  next if ($num_events_per_gene < 2);

  ## threshold for prevalence per sample (n)
  my $hist_prev_thr = 2; ## must be recurrent event in histology (n>=2)
  my $prev_in_hist = 0;
  if ($hist_check_gene{$event}) {
      $prev_in_hist = $hist_check{$event};
  }

  print OUT $event;

  foreach my $hist (sort @broad_hist_uniq) {
      my $i = 1;
      foreach my $sample (@{$histology_ids{$hist}}) {
          print OUT "\t";
          my $greatest_psi = max(@{$inc_levels_gene{$event}{$sample}});
          if ($greatest_psi > 0) {
              print OUT $greatest_psi;
          } else {
              print OUT "0";
          }
          $i++;
      }
  }
  print OUT "\n";
}
