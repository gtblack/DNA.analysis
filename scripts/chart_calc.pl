#!usr/bin/perl

# The result of the alignment scores are so scattered that we cannot predict
# considerably accurate points of alignment for the sequences.
# Therefore we are now looking for the most possible match by distributing
# parts of score to the nearby areas. This would be .95*(score)*(distance).
# Alignment scores with negative positions will be discarded.

use strict;
use warnings;

my %cscores;

open SCORE, "alignment_scores.txt";

my @ref_95;

foreach my $power (1..150) {
  $ref_95[$power] = .95 ** $power;
}

while (<SCORE>) {
  chomp;
  my @line = split /\t/, $_;
  if ($line[2] > 0 && $line[0] ne "") {
    $cscores{"$line[0]_$line[1]_$line[2]"} += $line[4];
    foreach my $dist (1..150) {
      my $mspos = $line[2] - $dist;
      my $pspos = $line[2] + $dist;
      $cscores{"$line[0]_$line[1]_$mspos"} += $ref_95[$dist] * $line[4];
      $cscores{"$line[0]_$line[1]_$pspos"} += $ref_95[$dist] * $line[4];
    }
  }
}

close SCORE;

my $file_indicator = 0;

open OUTPUT, ">chartscores/scores_0.txt";

foreach my $key (sort keys %cscores) {
  my @elem = split /_/, $key;
  if ($elem[0] eq $file_indicator) {
    close OUTPUT;
    open OUTPUT, ">chartscores/scores_$elem[0].txt";
  }
  print OUTPUT "$elem[1]_$elem[2]\t$cscores{$key}\n";
}

close OUTPUT;
