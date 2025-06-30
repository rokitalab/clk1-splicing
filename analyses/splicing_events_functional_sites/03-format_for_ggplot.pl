#!/usr/bin/perl

use strict;
use warnings;

# Generate for ggplot for all events corresponding to functional sites
open(my $out_pos, '>', 'results/splicing_events.SE.total.pos.intersectunip.ggplot.txt') or die "Cannot open file: $!";
print $out_pos "SpliceID\tdPSI\tUniprot\n";

open(my $out_neg, '>', 'results/splicing_events.SE.total.neg.intersectunip.ggplot.txt') or die "Cannot open file: $!";
print $out_neg "SpliceID\tdPSI\tUniprot\n";

sub process_file {
    my ($prefix, $suffix, $category, $output) = @_;

    open(my $in, '<', "results/splicing_events.SE.total.$prefix.intersect_UP000005640_9606_$suffix.wo.txt") or die "Cannot open file: $!";
    my %seen;
    while (my $line = <$in>) {
        chomp $line;
        my ($splice_id, $dpsi) = (split /\s+/, $line)[3, 4];
        next if $seen{$splice_id}++; # Skip duplicates
        print $output "$splice_id\t$dpsi\t$category\n";
    }
    close($in);
}

# Process positive files
process_file('pos', 'mod_res', 'Modifications', $out_pos);
process_file('pos', 'domain', 'Domains', $out_pos);
process_file('pos', 'disulfid', 'DisulfBond', $out_pos);
process_file('pos', 'signal', 'LocSignal', $out_pos);

# Process negative files
process_file('neg', 'mod_res', 'Modifications', $out_neg);
process_file('neg', 'domain', 'Domains', $out_neg);
process_file('neg', 'disulfid', 'DisulfBond', $out_neg);
process_file('neg', 'signal', 'LocSignal', $out_neg);

close($out_pos);
close($out_neg);
