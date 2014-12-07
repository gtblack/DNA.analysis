#!/usr/bin/perl

# The exclusive file contains class III protein-ORF pairs, and we need to
# separate class IV and V pairs from them.

use strict;
use warnings;

my %pairs;
my %cpairs;

open EXC, "BLASTp_exclusive.txt";
while (<EXC>) {
  chomp;
  my $line = $_;
  my @block = split /\s/, $line;
  if ($block[0] ne 'Query') {
    $pairs{$block[0]} = $line;
  }
}
close EXC;

open CONS, "BLASTp_consecutive.txt";
while (<CONS>) {
  chomp;
  my $line = $_;
  my $result = ($line =~/^(?<pi25611>[^\s]+)\s+/);
  if ($result) {
    $cpairs{$+{pi25611}} = $line;
  }
}
close CONS;

open OUTPUT, ">BLASTp_nonmatching.txt";
print OUTPUT "Query\t\t\tDatabase\tScore\tE Value\n";
foreach my $key (sort keys %pairs) {
  if (!(defined $cpairs{$key})) {
    print OUTPUT "$pairs{$key}\n";
  }
}
close OUTPUT;
