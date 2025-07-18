#!/usr/bin/perl

use Statistics::Lite qw(:all);
#use warnings;
############################################################################################################
# 04-generate_hist_spec_events_tab.pl
#
# ./04-generate_hist_spec_events_tab.pl <hist file> <rmats file> <indep_samples-plus>
############################################################################################################
my ($histology,$rmats_tsv,$primary_tumor_dat) = ($ARGV[0], $ARGV[1],$ARGV[2]);
my (@broad_hist, @bs_id, @splicing_events);
my (%histology_ids, %inc_levels, %bs_id_hist, %hist_check, %hist_count);
my %splicing_psi;

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

print "annotate and store histology information... ".localtime(time)."\n";
  # annotate histology file #
    # hash with BS and disease
    # make arrays of each histology of BS IDs

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


print "process rMATs information... ".localtime(time)."\n";
my (%splice_totals_per_sample, %splice_totals);
## process rMATS output (may take awhile) from merged input file
#print "processing rMATs results...\n";
open(FIL, "gunzip -c $rmats_tsv |") || die ("can’t open $rmats_tsv");
while(<FIL>)
{
  chomp;

  ##only look at exon splicing events
  next unless($_=~/^SE/);

  my @cols  = split "\t";
  my $bs_id = $cols[1];
  my $ctrl  = $cols[2];

  next unless  $bs_id_hist{$bs_id};
  ## get gene name
  my $gene         = $cols[4];
  $gene=~s/\"//g; # remove quotation marks

  ## retrieve exon coordinates
  my $chr          = $cols[5];
  my $str          = $cols[6];
  my $exonStart    = $cols[7]+1; ## its 0-based so add 1
  my $exonEnd      = $cols[8];
  my $upstreamES   = $cols[15];
  my $upstreamEE   = $cols[16];
  my $downstreamES = $cols[17];
  my $downstreamEE = $cols[18];

  ## remove decimals if needed
  $exonStart =~s/\.0//;
  $exonEnd =~s/\.0//;
  $upstreamES =~s/\.0//;
  $upstreamEE =~s/\.0//;
  $downstreamES =~s/\.0//;
  $downstreamEE =~s/\.0//;

  ## retrieve inclusion level and junction count info
  my $inc_level = $cols[33];
  my $IJC        = $cols[25];
  my $SJC        = $cols[26];

  #get lengths of exons
  my $inc_len  = $cols[29];
  my $skip_len = $cols[30];
  my $thr_diff = $cols[-1];

  ## create unique ID for splicing change
  my $splice_id = $gene.":".$exonStart."-".$exonEnd."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE;
  $inc_levels{$splice_id}{$bs_id} = $inc_level;

  #print $splice_id,"\t",$inc_level,"\n";
  push @splicing_events, $splice_id;
  push @{$splicing_psi{$splice_id}}, $inc_level;

  $splice_totals_per_sample{$bs_id}++;
  $splice_totals{$splice_id}++;
}
close(FIL);

#print "## make sure all arrays are unique...\n";
## make sure all arrays are unique
my @broad_hist_uniq      = do { my %seen; grep { !$seen{$_}++ } @broad_hist };
my @bs_ids_uniq          = do { my %seen; grep { !$seen{$_}++ } @bs_ids };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };

#print "calculate mean and sd for each splicing event...\n";
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

## identify aberrant splicing events
my %absplice_totals_per_sample;
my %absplice_totals_per_sample_pos;
my %absplice_totals_per_sample_neg;
my %splice_event_per_neg_hist_count;
my @ab_splicing_events_pos;
my @ab_splicing_events_neg;

my (%check_events_neg, %check_events_pos);
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
      my $histology = $bs_id_hist{$sample};
      #print $splice_event,"\t",$histology,"***\n";
      $splice_event_per_pos_hist_count{$splice_event}{$histology}++;
      push @ab_splicing_events_pos, $splice_event;

    }
    # < -2 z-scores
    if($psi_tumor < ($mean_psi - ($std_psi + $std_psi)) )
    {
      $absplice_totals_per_sample_neg{$sample}++;
      my $histology = $bs_id_hist{$sample};
      $splice_event_per_neg_hist_count{$splice_event}{$histology}++;
      push @ab_splicing_events_neg, $splice_event;

    }
  }
}

print "report results...".localtime(time)."\n";

my @ab_splicing_events_pos_uniq = do { my %seen; grep { !$seen{$_}++ } @ab_splicing_events_pos };
my @ab_splicing_events_neg_uniq = do { my %seen; grep { !$seen{$_}++ } @ab_splicing_events_neg };

my $output_file = "results/recurrent_splice_events_by_histology.tsv";
open(TAB,">".$output_file) || die("Cannot Open File");

print TAB "splicing_event\thistology\ttype\tfreq\n";

## save and report skipping events
foreach $hist (@broad_hist_uniq)
{
  foreach my $event (@ab_splicing_events_pos_uniq)
  {
        my $total_hist_count = $hist_count{$hist};

        if($splice_event_per_pos_hist_count{$event}{$hist}){
          my $event_count = $splice_event_per_pos_hist_count{$event}{$hist};
          #print $event,"\t",$hist,"\t",$total_hist_count,"\n";
          #if( ($event_count/$total_hist_count) >= .10 )
          if( ($event_count) >= 2 )


          {
            print TAB $event,"\t",$hist,"\tinclusion\t";
            print TAB $event_count,"\n";


          }
        }
  }
}

## save and report inclusion events
foreach $hist (@broad_hist_uniq)
{
  foreach my $event (@ab_splicing_events_neg_uniq)
  {
        my $total_hist_count = $hist_count{$hist};
        if($splice_event_per_neg_hist_count{$event}{$hist}){
          my $event_count = $splice_event_per_neg_hist_count{$event}{$hist};
          #print $event,"\t",$hist,"\t",$total_hist_count,"*\n";
          #if( ($event_count/$total_hist_count) >= .10 )
          if( ($event_count) >= 2 )

          {
            print TAB $event,"\t",$hist,"\tskipping\t";
            print TAB $event_count,"\n";

          }
        }
  }
}

close(TAB);
print "finish... ".localtime(time)."\n";
