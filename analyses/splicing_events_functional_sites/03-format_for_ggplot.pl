#!/usr/bin/perl

use strict;
use warnings;

open(my $out_pos, '>', 'results/splicing_events.total.pos.intersectunip.ggplot.txt') or die "Cannot open file: $!";
print $out_pos "splicing_case\tSpliceID\tdPSI\tUniprot\n";
  
open(my $out_neg, '>', 'results/splicing_events.total.neg.intersectunip.ggplot.txt') or die "Cannot open file: $!";
print $out_neg "splicing_case\tSpliceID\tdPSI\tUniprot\n";
  
# Generate for ggplot for all events corresponding to functional sites
sub process_file {
    my ($prefix, $suffix, $category, $output, $case) = @_;
  
    open(my $in, '<', "results/splicing_events.$case.total.$prefix.intersect_UP000005640_9606_$suffix.wo.txt") or die "Cannot open file: $!";
    my %seen;
    while (my $line = <$in>) {
        chomp $line;
        my ($splice_id, $dpsi) = (split /\s+/, $line)[3, 4];
        next if $seen{$splice_id}++; # Skip duplicates
        print $output "$case\t$splice_id\t$dpsi\t$category\n";
    }
    close($in);
}

# Loop through all splicing event cases
my @cases = ("SE", "RI", "A3SS", "A5SS");
foreach my $case (@cases) {
  # Process positive files
  process_file('pos', 'mod_res', 'Modifications', $out_pos, $case);
  process_file('pos', 'domain', 'Domains', $out_pos, $case);
  process_file('pos', 'disulfid', 'DisulfBond', $out_pos, $case);
  process_file('pos', 'signal', 'LocSignal', $out_pos, $case);
  
  # Process negative files
  process_file('neg', 'mod_res', 'Modifications', $out_neg, $case);
  process_file('neg', 'domain', 'Domains', $out_neg, $case);
  process_file('neg', 'disulfid', 'DisulfBond', $out_neg, $case);
  process_file('neg', 'signal', 'LocSignal', $out_neg, $case);
}

close($out_pos);
close($out_neg);