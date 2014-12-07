#!usr/bin/perl

use strict;
use warnings;

my %pairs;
my %scores;
my %evalues;
my $query_tag;
my $flag = 0;

open BLAST, "PI17.blastp.result.txt";
while (<BLAST>) {
  chomp;
  my $line = $_;
  if ($line =~ m/^Query= (?<qtag>\w+)/) {
    if ($flag > 0) {
      $flag = 0;
    }
    $flag++;
    $query_tag = $+{qtag};
  } elsif (($line =~ m/ gi\|(?<rtag>[^\|]+)[^\s]+ .{33} +(?<scr>[^\s]+) +(?<ev>[^\s]+)/) && ($flag == 1)) {
    $flag = 0;
    $pairs{$query_tag} = $+{rtag};
    $scores{$query_tag} = $+{scr};
    $evalues{$query_tag} = $+{ev};
  }
}
close BLAST;

my %ev_decimal;

foreach my $key (keys %evalues) {
  $ev_decimal{$key} = sprintf("%.10g", $evalues{$key});
}

open OUTPUT, ">BLASTp_reverse.txt";
print OUTPUT "PI25611\tPI17\tBit Score\tE value\n";
foreach my $key (sort {$ev_decimal{$b} <=> $ev_decimal{$a}} (keys %evalues)) {
  print OUTPUT "$key\t$pairs{$key}\t$scores{$key}\t$evalues{$key}\n";
}
close OUTPUT;
