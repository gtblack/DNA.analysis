#!usr/bin/perl

use strict;
use warnings;

my %anchor_pairs;
my %align_pairs;
my %anchor_line;

open ANCH, "BLASTp_mauve_best.txt";
while (<ANCH>) {
  chomp;
  my $line = $_;
  my @block = split /\t/, $line;
  if ($block[0] ne "PI25611") {
    $anchor_pairs{$block[0]} = $block[1];
    $anchor_line{$block[0]} = $line;
  }
}
close ANCH;

open ALIGN, "BLASTp_best_matches.txt";
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
my %anchor_only;

foreach my $key (keys %anchor_pairs) {
  if (defined $align_pairs{$key}) {
    $count_both++;
  } else {
    $count_anch++;
    $anchor_only{$key} = $anchor_pairs{$key};
  }
}

$count_align = (keys %align_pairs) - $count_both;

print "$count_both pairs are in both files.\n";
print "$count_anch pairs are only in the first file.\n";
print "$count_align pairs are only in the second file.\n";

open OUTPUT, ">MAUVE_unique.txt";
print OUTPUT "PI25611\t\tPI17\t\tScore\tE Value\n";
foreach my $key (sort {$anchor_only{$a} <=> $anchor_only{$b}} (keys %anchor_only)) {
  print OUTPUT "$anchor_line{$key}\n";
}
close OUTPUT;
