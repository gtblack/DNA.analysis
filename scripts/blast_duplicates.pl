#!/usr/bin/perl

# While most of the protein-ORF pairs are unique, there are some pairs that
# share the PI17 proteins. We need to address those pairs manually, but
# before that we need to figure out which pairs share the proteins.

use strict;
use warnings;

my %pairs;

open BLAST, "BLASTp_best_matches.txt";
while (<BLAST>) {
  chomp;
  my @block = split /\t/, $_;
  if ($block[0] =~ /^pint/) {
    if (defined $pairs{$block[1]}) {
      $pairs{$block[1]} = [@{$pairs{$block[1]}}, $block[0]];
    } else {
      $pairs{$block[1]} = [$block[0]];
    }
  }
}
close BLAST;

open OUTPUT, ">class1.duplicates.txt";
foreach my $key (sort keys %pairs) {
  my $nums = scalar (@{$pairs{$key}});
  if ($nums > 1) {
    print OUTPUT "Duplicate ORFs with best match protein $key :\n";
    foreach my $pkey (1..$nums) {
      print OUTPUT (${$pairs{$key}}[($pkey - 1)] . "\t");
    }
    print OUTPUT "\n\n";
  }
}
close OUTPUT;
