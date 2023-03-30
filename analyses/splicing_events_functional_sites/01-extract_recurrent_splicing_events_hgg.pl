#!/usr/bin/perl

######################################################################################
# extract_recurrent_splicing_events.pl
# descr: extract and append splicing ids from rMATS to be used to overlap with Uniprot
# usage: perl extract_recurrent_splicing_events.pl <pbta-histologies.tsv>
######################################################################################

## input argument for histology file from OpenPBTA project
my $histology = $ARGV[0]; ## ../../data/v19_plus_20210311_pnoc_rna.tsv

## merged rMATS file generated by previous script/hardcoded
my $rmats_tsv = "../../data/rMATS_merged.comparison.tsv.gz";

## data structures to store values
my (@splicing_events, @samples);
my %splicing_event_counts;
my %inc_levels;
my (%chr,%num_exons,%str);
my %hgg_midline_samples;

## output to terminal
print "processing histology...\n";
my %flip_event;

## open histology file and annotate DMG samples
open(FIL,$histology) || die("Cannot Open File");
while(<FIL>)
{
  chomp;

  ## only select RNA-seq and PBTA samples
  next unless ($_=~/RNA-Seq/);
  next unless ($_=~/PBTA/);

  my @cols       = split "\t";
  my $hist       = $cols[44];
  my $bs_id      = $cols[0];
  my $patient_id = $cols[3];
  my $CNS_region = $cols[32];

  ##keep only HGG midline samples
  next unless $hist=~/HGAT/;
  next unless($CNS_region=~/Midline/);

  $hgg_midline_samples{$bs_id} = $bs_id;
}
close(FIL);

## process rMATS output (may take awhile) from merged input file
print "processing rMATs results...\n";
open(FIL, "gunzip -c $rmats_tsv |")|| die("Cannot Open File");
while(<FIL>)
{
  chomp;

  ##only look at exon splicing events
  next unless($_=~/^SE/);

  my @cols = split "\t";
  my $sample = $cols[1];
  my $ctrl   = $cols[2];

  ## skip non-midline HGGs and non BS ctrl comparisons
  next unless $hgg_midline_samples{$sample};
  next unless $ctrl =~/control\-BS/;
  #print "$sample\n";
  #print "bs_id: ",$sample,"\n";

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

  ## retrieve inclusion level and junction count info
  my $inc_level_ctrl = $cols[33];
  my $inc_level_tumor = $cols[34];
  my $ctrl_IJC        = $cols[25];
  my $ctrl_SJC        = $cols[26];
  my $tumor_IJC       = $cols[27];
  my $tumor_SJC       = $cols[28];

  #get lengths of exons
  my $inc_len  = $cols[29];
  my $skip_len = $cols[30];

  ## get dPSI values
  my $thr_diff = $cols[-1];

  ## only look at strong changes, remove any dPSI < .20 and tumor junction reads > 10 reads
  next unless ( abs($thr_diff) >= .20);
  next unless ( ( ($tumor_IJC + $tumor_SJC) >=10) ) ;



  # annotate and re-name splice event IDs
  my $splice_id= $gene."_".$exonStart."-".$exonEnd."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE;
  my $splice_id_short= $gene."_".$exonStart."-".$exonEnd;
  #print $sample,"\t",$splice_id_short,"\t",$splice_id,"\n";

  ##keep one splicing event per exon (choose larger dPSI), sometimes with different upstream/downstream exons
  my $related_event = 0;
  if($inc_levels{$splice_id_short}{$sample})
  {
    #print "*",$sample,"\t",$splice_id_short,"\t",$splice_id,"\t",$thr_diff,"\n";
    $related_event = $inc_levels{$splice_id_short}{$sample};
    #print "*",$sample,"\t",$splice_id_short,"\t",$splice_id,"\t",$thr_diff,"\t",$related_event,"\n";
  }
  #print $sample,"\t",$splice_id_short,"\t",$splice_id,"\tthr:",$thr_diff,"\trelated:",$related_event,"\n";
  next unless(abs($thr_diff) > abs($related_event));
  #print "splice_id_short:".$splice_id_short,"\tdPSI:",$thr_diff,"E\tRE:",$related_event,"\n";

  $chr{$splice_id_short} = $chr;
  $str{$splice_id_short} = $str;
  $inc_levels{$splice_id_short}{$sample} = $thr_diff;
  #print "keep: ",$sample,"\t",$splice_id_short,"\t",$splice_id,"\t",$thr_diff,"\n";

  if( ($inc_level_ctrl >= .80) && ($inc_level_tumor <=.10) || ($inc_level_ctrl <= .10) && ($inc_level_tumor >=.80) )
  {
    $flip_event{$splice_id_short} = 1;
  }
  else{
    $flip_event{$splice_id_short} = 0;

  }


  ## store splicing events, inclusion levels, and splicing event counts/freqs
  push @splicing_events, $splice_id_short;
  push @samples,$sample;
}
close(FIL);


## make arrays unique
my @samples_uniq = do { my %seen; grep { !$seen{$_}++ } @samples };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };

## hashes to store and keep dPSIs of splicing events separating by pos and neg
my %splicing_event_deltapsis_pos;
my %splicing_event_deltapsis_neg;
my %splicing_event_counts_pos;
my %splicing_event_counts_neg;

## push all deltaPSI values into associated splicing event
foreach my $event(@splicing_events_uniq)
{
  foreach my $sample(@samples_uniq)
  {
    next unless $inc_levels{$event}{$sample};
    my $dPSI = $inc_levels{$event}{$sample};

    ## divide into postive and negative dPSIs
    if($dPSI>0){
      $splicing_event_counts_pos{$event}++;
      push @{$splicing_event_deltapsis_pos{$event}}, $dPSI;
    }
    else{
      $splicing_event_counts_neg{$event}++;
      push @{$splicing_event_deltapsis_neg{$event}}, $dPSI;
    }
  }
}

## create results folder if it does not exists
if (!-d "results") { mkdir "results"; }

## write to files for summary tables (with a bed version for downstream analyses)
## print to one tsv table and a bed table
open(TABPOS, ">results/splicing_events.total.pos.tsv");
open(BEDPOS, ">results/splicing_events.total.pos.bed");

open(TABNEG, ">results/splicing_events.total.neg.tsv");
open(BEDNEG, ">results/splicing_events.total.neg.bed");

print "writing output...\n";
print TABPOS "gene\tsplice_event\tavg_dpsi\tstdev_dpsi\tfreq\tcoord\tstr\tflip\n";
print TABNEG "gene\tsplice_event\tavg_dpsi\tstdev_dpsi\tfreq\tcoord\tstr\tflip\n";

foreach my $event (@splicing_events_uniq)
{

  ## compute avg and stdevs of each lsv dPSI for positive dPSIs
  if(@{$splicing_event_deltapsis_pos{$event}})
  {
    ## only keep recurring events
    next unless $splicing_event_counts_pos{$event} >= 2;

    ## compute avg and std dPSI
    my $avg_dpsi_pos = &average(\@{$splicing_event_deltapsis_pos{$event}});
  	my $std_dpsi_pos = &stdev(\@{$splicing_event_deltapsis_pos{$event}});

    ## get the gene name and exon coords
    my ($gene,$exon_coord) = split/\_/,$event;

    ## print results
    print TABPOS $gene,"\t",$event,"\t",$avg_dpsi_pos,"\t",$std_dpsi_pos,"\t";
    print TABPOS $#{$splicing_event_deltapsis_pos{$event}}+1,"\t";
    print TABPOS "\t",$chr{$event},":",$exon_coord,"\t",$str{$event},"\t";
    print TABPOS  $flip_event{$event},"\n";

    print BEDPOS $chr{$event},"\t";

    ## get exon coords for bed to intersect w Uniprot later
    if($event=~/(\w+)\_(\d+\-\d+)/)
    {
      $exon_coord = $2;
      my($start_exon_coord,$end_exon_coord) = split/\-/,$exon_coord;
      print BEDPOS $start_exon_coord,"\t",$end_exon_coord,"\t",$event,"\t",$avg_dpsi_pos,"\t",$str{$event},"\t";
      print BEDPOS  $flip_event{$event},"\n";

    }
  }
  if(@{$splicing_event_deltapsis_neg{$event}})
  {
    ## only keep recurring events
    next unless $splicing_event_counts_neg{$event} >= 2;

    ## compute avg and std dPSI
    my $avg_dpsi_neg = &average(\@{$splicing_event_deltapsis_neg{$event}});
  	my $std_dpsi_neg = &stdev(\@{$splicing_event_deltapsis_neg{$event}});

    ## get the gene name and exon coords
    my ($gene,$exon_coord) = split/\_/,$event;

    ## print results
    print TABNEG $gene,"\t",$event,"\t",$avg_dpsi_neg,"\t",$std_dpsi_neg,"\t";
    print TABNEG $#{$splicing_event_deltapsis_neg{$event}}+1,"\t";
    print TABNEG "\t",$chr{$event},":",$exon_coord,"\t",$str{$event},"\t";
    print TABNEG  $flip_event{$event},"\n";


    print BEDNEG $chr{$event},"\t";

    ## get exon coords for bed to intersect w Uniprot later
    if($event=~/(\w+)\_(\d+\-\d+)/)
    {
      $exon_coord = $2;
      my($start_exon_coord,$end_exon_coord) = split/\-/,$exon_coord;
      print BEDNEG $start_exon_coord,"\t",$end_exon_coord,"\t",$event,"\t",$avg_dpsi_neg,"\t",$str{$event},"\t";
      print BEDNEG  $flip_event{$event},"\n";

    }
  }
}
close(TABPOS);
close(BEDPOS);
close(TABNEG);
close(BEDNEG);

## function to calculate stdev given an array of dPSI values
sub stdev{
        my($data) = @_;
        if(@$data == 1){

                return "1";
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

## function to calculate avg given an array of dPSI values
sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty arrayn");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}