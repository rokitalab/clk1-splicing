#!/usr/bin/perl

my @psi_expr = <scr/*_expr.tab.txt>;

foreach my $file (@psi_expr)
{
  print $file,"\n";
  my $entries = `wc -l $file`;
  next if $entries < 21;
  system("awk -f stdev.awk -f pearson.awk $file $file")
}