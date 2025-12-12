#!/usr/bin/perl

use Statistics::Lite qw(:all);
############################################################################################################
# 01-extract_recurrent_splicing_events_cluster.pl
#
# Compute splicing index for each sample and generate splicing burden index tables
############################################################################################################
my ($cluster_file,$rmats_tsv, $splice_case) = ($ARGV[0], $ARGV[1], $ARGV[2]);
my (@broad_hist, @bs_id, @splicing_events);
my (%histology_ids, %inc_levels, %bs_id_hist, %hist_check, %hist_count);
my @splicing_events;
my %splicing_psi;
my %cns_regions;
my @clusters;
my (%bs_id_cluster,%cluster_count);

$splice_case=~s/\s+//;
print $splice_case,"\n";

## check to see if splicing case is correct
unless ($splice_case=~/SE$|A3SS$|A5SS$|RI$/)
{
  die("Splicing case does not exist, please use either 'SE', 'A5SS', 'A3SS', or 'RI'");
}
# annotate histology file #
# hash with BS and disease
# make arrays of each histology of BS IDs

my %primary_initial_sample_list;



## annotate and store histology information
open(FIL, $cluster_file) || die("Cannot Open File");

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
  my $cluster_assigned      = $cols[$column_index{'cluster'}];
  my $bs_id      = $cols[$column_index{'Kids_First_Biospecimen_ID'}];

  #next unless $cluster_assigned=~/6/;
  
  ## make an array and store histology information and BS IDs
  push @clusters, $cluster_assigned;
  push @bs_ids, $bs_id;

  $bs_id_cluster{$bs_id} = $cluster_assigned;

  ## store total number of histologies
  push @{$cluster_ids{$cluster_assigned}}, $bs_id;

  $cns_regions{$bs_id} = $molecular_subtype_display;

  ## histology counter for downstream analysis
  $cluster_count{$cluster_assigned}++;

}
close(FIL);


## process rMATs output
my (%splice_totals_per_sample, %splice_totals,%chr, %str);
my %splice_event_hist;

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
  my $ctrl  = $cols[2];

  ## filter for samples assigned to a cluster
  next unless $bs_id_cluster{$bs_id};
  my $cluster = $bs_id_cluster{$bs_id};

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
    #print $Start,"\t",$End,"\t",$prevES,"\t",$prevEE,"\t",$nextES,"\t",$nextEE,"\n";

  }

  ## retrieve inclusion level and junction count info
  my $inc_level = $cols[33];
  my $IJC        = $cols[25];
  my $SJC        = $cols[26];

  #get lengths of exons
  my $inc_len  = $cols[29];
  my $skip_len = $cols[30];

  ## create unique ID for splicing change
  my $splice_id = $gene.":".$Start."-".$End."_".$prevES."-".$prevEE."_".$nextES."-".$nextEE;

  ##remove .0 in coord
  $splice_id=~s/\.0//g;

  ## chr and strand
  $chr{$splice_id} = $chr;
  $str{$splice_id} = $str;

  ## store PSI for event and sample
  $inc_levels{$splice_id}{$bs_id} = $inc_level;

  push @splicing_events, $splice_id;
  push @{$splicing_psi{$splice_id}}, $inc_level;

  $splice_totals_per_sample{$bs_id}++;
  $splice_totals{$splice_id}++;

  $splice_event_hist{$splice_id}{$cluster}=$cluster;
  #print $splice_id,"\n";
}
close(FIL);

print "## make sure all arrays are unique...\n";
## make  all arrays are unique
my @cluster_uniq      = do { my %seen; grep { !$seen{$_}++ } @clusters };
my @bs_ids_uniq          = do { my %seen; grep { !$seen{$_}++ } @bs_ids };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };
my (%std_dev_psi,%count_psi,%mean_psi);

print "## calculate mean and sd for each splicing event...\n";

## calculate mean and sd for each splicing event
foreach my $splice_event(@splicing_events_uniq)
{
  my $mean     = mean  (@{$splicing_psi{$splice_event}});
  my $count    = count (@{$splicing_psi{$splice_event}});
  my $std_dev  = stddev(@{$splicing_psi{$splice_event}});

  $std_dev_psi{$splice_event} = $std_dev;
  $count_psi{$splice_event}   = $count;
  $mean_psi{$splice_event}    = $mean;
}

## assess each tumor samples to identify if it is aberrant (2 standard deviations from the mean)
if($subtype =~/\w+/){
  open(EVENTS,">results/splice_events.diff.".$splice_case.".".$subtype.".txt");
  print EVENTS "Splice ID\tCase\tSample\tCluster\tCNS\tType\tdPSI\n";
  open(BEDPOS, ">results/splicing_events.SE.total.".$subtype.".pos.bed");
  open(BEDNEG, ">results/splicing_events.SE.total.".$subtype.".neg.bed");
}
else{
  open(EVENTS,">results/splice_events.diff.".$splice_case.".txt");
  print EVENTS "Splice ID\tCase\tSample\tCluster\tCNS\tType\tdPSI\n";
  open(BEDPOS, ">results/splicing_events.".$splice_case.".total.pos.bed");
  open(BEDNEG, ">results/splicing_events.".$splice_case.".total.neg.bed");
}
foreach my $sample(@bs_ids_uniq)
{
  # Only reports on cluster 6 values
  # next unless $bs_id_cluster{$sample} == 6;
  foreach my $splice_event(@splicing_events_uniq)
  {
    ## look at splicing events that are present in sample and are recurrent
    next unless $inc_levels{$splice_event}{$sample};
    next if $count_psi{$splice_event} < 2;
    next unless $mean_psi{$splice_event};

    ## report only histology of interest if appropriate
    #next unless $bs_id_hist{$sample}=~/HGAT/;

    my $psi_tumor = $inc_levels{$splice_event}{$sample};
    my $std_psi   = $std_dev_psi{$splice_event};
    my $mean_psi  = $mean_psi{$splice_event};

    #> +2 z-scores
    if($psi_tumor > ($mean_psi + ($std_psi + $std_psi)) )
    {
      
      # dPSI filter
      my $dpsi = $psi_tumor-$mean_psi;
      next unless $dpsi > 0.1;
      
      print EVENTS $splice_event,"\t".$splice_case,"\t",$sample,"\t",$bs_id_cluster{$sample},"\t",$cns_regions{$sample},"\tInclusion\t",$dpsi,"\n";
      print BEDPOS $chr{$splice_event},"\t";

      ## get exon coords for bed to intersect w Uniprot later
      if($splice_event=~/(\w+)\:(\d+\-\d+)/)
      {
        $exon_coord = $2;
        my($start_exon_coord,$end_exon_coord) = split/\-/,$exon_coord;
        print BEDPOS $start_exon_coord,"\t",$end_exon_coord,"\t",$splice_event,"\t",$dpsi,"\t",$str{$splice_event},"\n";

      }
    }
    # < -2 z-scores, inclusion events
    if($psi_tumor < ($mean_psi - ($std_psi + $std_psi)) )
    {
      
      # dPSI filter
      my $dpsi = $mean_psi-$psi_tumor;
      next unless $dpsi > 0.1;
      
      print EVENTS $splice_event,"\t".$splice_case,"\t",$sample,"\t",$bs_id_cluster{$sample},"\t",$cns_regions{$sample},"\tSkipping\t",$dpsi,"\n";
      print BEDNEG $chr{$splice_event},"\t";

      ## get exon coords for bed to intersect w Uniprot later
      if($splice_event=~/(\w+)\:(\d+\-\d+)/)
      {
        $exon_coord = $2;
        my($start_exon_coord,$end_exon_coord) = split/\-/,$exon_coord;
        print BEDNEG $start_exon_coord,"\t",$end_exon_coord,"\t",$splice_event,"\t",$dpsi,"\t",$str{$splice_event},"\n";
      }
    }
  }
}
close(EVENTS);
close(BEDPOS);
close(BEDNEG);
