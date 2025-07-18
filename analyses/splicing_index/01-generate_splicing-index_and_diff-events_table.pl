#!/usr/bin/perl

use Statistics::Lite qw(:all);
use IO::Compress::Gzip qw(gzip $GzipError);
############################################################################################################
# 03-generate_splicing-index_and_diff-events_table.SE.pl
#
# Compute splicing index for each sample and generate splicing burden index tables
############################################################################################################
my ($histology,$rmats_tsv,$primary_tumor_dat, $splice_case) = ($ARGV[0], $ARGV[1], $ARGV[2],$ARGV[3]);
my (@broad_hist, @bs_id, @splicing_events);
my (%histology_ids, %inc_levels, %bs_id_hist, %hist_check, %hist_count);
my @splicing_events;
my %splicing_psi;
my (%splice_totals_per_sample, %splice_totals);

$splice_case=~s/\s+//;
print $splice_case,"\n";

## check to see if splicing case is correct
unless ($splice_case=~/SE$|A3SS$|A5SS$|RI$/)
{
  die("Splicing case does not exist, please use either 'SE', 'A5SS', 'A3SS', or 'RI'");
}
my %primary_initial_sample_list;

## store primary tumor samples
open(FIL,$primary_tumor_dat) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my @header = split "\t";
  my $bs_id = $header[1];
  my $cohort = $header[2];
  my $exp_strategy = $header[4];
  my $tumor_descr = $header[5];

  next unless ($cohort=~/PBTA/);

  $primary_initial_sample_list{$bs_id} = $bs_id;

}
close(FIL);

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
  my @cols       = split "\t";
  my $hist      = $cols[$column_index{'plot_group'}];
  my $bs_id      = $cols[$column_index{'Kids_First_Biospecimen_ID'}];
  my $CNS_region = $cols[$column_index{'CNS_region'}];
  my $RNA_libr   = $cols[$column_index{'RNA_library'}];

  next unless ($primary_initial_sample_list{$bs_id});
  next unless ($RNA_libr eq "stranded");

  ## make an array and store histology information and BS IDs
  push @broad_hist, $hist;
  push @bs_ids, $bs_id;

  $bs_id_hist{$bs_id} = $hist;

  ## store total number of histologies
  $hist_count{$hist}++;
  push @{$histology_ids{$hist}}, $bs_id;

  $cns_regions{$bs_id} = $CNS_region;
}

close(FIL);

# build a lookup hash of the bs_ids saved
my %bs_ids_keep = map { $_ => 1 } @bs_ids;

## process rMATS output (may take awhile) from merged input file
print "processing rMATs results... ".(localtime)."\n";
open(FIL, "gunzip -c $rmats_tsv |") || die ("can’t open $rmats_tsv");
while(<FIL>)
{
  chomp;

  ##only look at exon splicing events
  next unless($_=~/^$splice_case/);

  my @cols  = split "\t";
  my $bs_id = $cols[1];
  
  next unless exists $bs_ids_keep{$bs_id};
  
  my $ctrl  = $cols[2];

  ## get gene name
  my $gene         = $cols[4];
  $gene=~s/\"//g; # remove quotation marks

  ## retrieve exon coordinates
  my $chr          = $cols[5];
  my $str          = $cols[6];

  my $Start    = "";
  my $End      = "";
  my $prevES   = "";
  my $prevEE   = "";
  my $nextES = "";
  my $nextEE = "";

  if($splice_case=~/SE/)
  {
    $Start    = $cols[7]+1; ## its 0-based so add 1
    $End      = $cols[8];
    $prevES   = $cols[15];
    $prevEE   = $cols[16];
    $nextES = $cols[17];
    $nextEE = $cols[18];
  }
  elsif($splice_case=~/RI/)
  {
    $Start    = $cols[13]+1; ## its 0-based so add 1
    $End      = $cols[14];
    $prevES   = $cols[15];
    $prevEE   = $cols[16];
    $nextES = $cols[17];
    $nextEE = $cols[18];
  }
  elsif($splice_case=~/A5SS/)
  {
    $Start    = $cols[19]+1; ## its 0-based so add 1
    $End      = $cols[20];
    $prevES   = $cols[21];
    $prevEE   = $cols[22];
    $nextES = $cols[23];
    $nextEE = $cols[24];

  }
  else{
    $Start    = $cols[19]+1; ## its 0-based so add 1
    $End      = $cols[20];
    $prevES   = $cols[21];
    $prevEE   = $cols[22];
    $nextES = $cols[23];
    $nextEE = $cols[24];
  }

  ## retrieve inclusion level and junction count info
  my $inc_level  = $cols[33];
  my $IJC        = $cols[25];
  my $SJC        = $cols[26];

  #get lengths of exons
  my $inc_len  = $cols[29];
  my $skip_len = $cols[30];
  my $thr_diff = $cols[-1];

  ## create unique ID for splicing change

  my $splice_id = $gene.":".$Start."-".$End."_".$prevES."-".$prevEE."_".$nextES."-".$nextEE;
  $splice_id=~s/\.0//g;

  $inc_levels{$splice_id}{$bs_id} = $inc_level;

  #print $splice_id,"\t",$inc_level,"\n";
  push @splicing_events, $splice_id;
  push @{$splicing_psi{$splice_id}}, $inc_level;

  $splice_totals_per_sample{$bs_id}++;
  $splice_totals{$splice_id}++;
}
close(FIL);

print "## make sure all arrays are unique...\n";
## make  all arrays are unique
my @broad_hist_uniq      = do { my %seen; grep { !$seen{$_}++ } @broad_hist };
my @bs_ids_uniq          = do { my %seen; grep { !$seen{$_}++ } @bs_ids };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };

print "## calculate mean and sd for each splicing event...\n";
## calculate mean and sd for each splicing event
my (%std_dev_psi, %count_psi, %mean_psi);
foreach my $splice_event(@splicing_events_uniq)
{
  my $mean     = mean  (@{$splicing_psi{$splice_event}});
  my $count    = count (@{$splicing_psi{$splice_event}});
  my $std_dev  = stddev(@{$splicing_psi{$splice_event}});

  $std_dev_psi{$splice_event} = $std_dev;
  $count_psi{$splice_event}   = $count;
  $mean_psi{$splice_event}    = $mean;
}

# Find abberrant splicing/differential events
my %absplice_totals_per_sample;
my %absplice_totals_per_sample_pos;
my %absplice_totals_per_sample_neg;

open(EVENTS,">results/splice_events.diff.".$splice_case.".txt");
print EVENTS "Splice ID\tCase\tType\tTumor_PSI\tMean_PSI\tSample\tHistology\n";
foreach my $sample(@bs_ids_uniq)
{
  foreach my $splice_event(@splicing_events_uniq)
  {
    next unless $inc_levels{$splice_event}{$sample};
    my $psi_tumor = $inc_levels{$splice_event}{$sample};
    my $std_psi   = $std_dev_psi{$splice_event};
    my $mean_psi   = $mean_psi{$splice_event};

    #> +2 z-scores
    if($psi_tumor > ($mean_psi + ($std_psi + $std_psi)) )
    {
      $absplice_totals_per_sample_pos{$sample}++;
      print EVENTS $splice_event,"\t".$splice_case."\tSkipping\t";
      print EVENTS $psi_tumor,"\t",$mean_psi,"\t",$sample,"\t";
      print EVENTS $bs_id_hist{$sample},"\n";

    }
    # < -2 z-scores
    if($psi_tumor < ($mean_psi - ($std_psi + $std_psi)) )
    {
      $absplice_totals_per_sample_neg{$sample}++;
      print EVENTS $splice_event,"\t".$splice_case."\tInclusion\t";
      print EVENTS $psi_tumor,"\t",$mean_psi,"\t",$sample,"\t";
      print EVENTS $bs_id_hist{$sample},"\n";

    }
  }
}
close(EVENTS);

#make table for plotting of splice_index
if (!-d "results")
{
  mkdir "results";
}

print "make table for plotting of splice_index ".$splice_case."...\n";
open(TAB,">results/splicing_index.".$splice_case.".txt");
print TAB "Sample\tTotal\tAS_neg\tAS_pos\tAS_total\tSI\tHistology\n";
foreach my $sample(@bs_ids_uniq)
{
  next unless $splice_totals_per_sample{$sample};
  my $total_absplice_totals_per_sample = $absplice_totals_per_sample_neg{$sample} + $absplice_totals_per_sample_pos{$sample};
  my $splice_index = $total_absplice_totals_per_sample/$splice_totals_per_sample{$sample};

  print TAB $sample,"\t";
  print TAB $splice_totals_per_sample{$sample},"\t";
  print TAB $absplice_totals_per_sample_neg{$sample},"\t",$absplice_totals_per_sample_pos{$sample},"\t";
  print TAB $total_absplice_totals_per_sample,"\t";
  print TAB $splice_index,"\t";
  print TAB $bs_id_hist{$sample},"\n";
}
close(TAB);
