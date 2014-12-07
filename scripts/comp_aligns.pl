#!usr/bin/perl

use strict;
use warnings;

my %anchor_pairs;
my %align_pairs;

open ANCH, "mauve_matches.txt";
while (<ANCH>) {
  chomp;
  my @block = split /\t/, $_;
  $anchor_pairs{$block[0]} = $block[1];
}
close ANCH;

open ALIGN, "BLASTp_aligned_matches.txt";
while (<ALIGN>) {
  chomp;
  my @block = split /\t/, $_;
  if ($block[0] ne "PI25611") {
    $align_pairs{$block[0]} = $block[1];
  }
}
close ALIGN;

my $count_anch = 0;
my $count_align = 0;
my $count_both = 0;

foreach my $key (keys %anchor_pairs) {
  if (defined $align_pairs{$key}) {
    $count_both++;
  } else {
    $count_anch++;
  }
}

$count_align = (keys %align_pairs) - $count_both;

print "$count_both pairs are in both files.\n";
print "$count_anch pairs are only in the first file.\n";
print "$count_align pairs are only in the second file.\n";
