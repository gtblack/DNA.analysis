#!usr/bin/perl

# Each contig of PI25611 has several anchors, some of which crosses
# the boundaries between contigs. Some are not so well aligned that
# the distances between the anchors in the reference and the shotgun
# might be more than 200k nucleotides long. To account this difference
# this script is going to readujust the contig locations defined by
# MAUVE by reading the start/end locations and then calculating
# the average distances between the reference and shotgun anchors in
# each contig so that the difference can be added to establish
# a final referential alignment text.

use strict;
use warnings;

open LOC, "PI25611.contigs.tab";

my %locs;
my $flag = 0;

# get the starting/ending locations of each contig relative to the reference
while (<LOC>) {
  chomp;
  my $line = $_;
  my @info = split /\t/, $line;
  if (/Ordered Contigs/) {
    $flag = 1;
  } elsif (($flag == 1) && ($info[0] eq "contig")) {
    $locs{$info[1]} = "$info[4]_$info[5]";
  } elsif ($line =~ m/Contigs with conflicting/) {
    $flag = 0;
  }
}

close LOC;

open ALIGN, "alignment3";

my %anchors;
my $temp_key;

# get each anchor's position in the reference and shotgun sequences
while (<ALIGN>) {
  chomp;
  my $line = $_;
  # set the position in the reference sequences as the key
  if (($line =~ /> 1:(?<key>[^\s]+)/) && ($flag == 0)) {
    $temp_key = $+{key};
    $flag = 1;
  # set the position in the shotgun sequences as the value
  } elsif (($line =~ /> 2:(?<key>[^\s]+)/) && ($flag == 1)) {
    $anchors{$temp_key} = $+{key};
    $flag = 0;
  # non-matched anchors are discarded
  } elsif (($line =~ /> 1:/) && ($flag == 1)) {
    $flag = 2;
  }
}

close ALIGN;

my @newlocs;

foreach my $chrom (1..24) {
  my @diffs;

  # get the relative location of each shotgun sequence
  $locs{"pint2554_c_$chrom"} =~ m/(?<st>\d+)_(?<end>\d+)/;
  my $cont_s = $+{st};
  my $cont_e = $+{end};

  foreach my $key (keys %anchors) {
    $anchors{$key} =~ m/(?<st>\d+)-(?<end>\d+)/;
    my $anch_s = $+{st};
    my $anch_e = $+{end};

    # the anchor should be in the range of the contig and the differences
    # between the shotgun should be not too far away.
    if (($cont_s <= $anch_e) && ($cont_e >= $anch_s)) {
      $key =~ m/(?<ref_st>\d+)-/;
      my $anch_diff = $+{ref_st} - $anch_s;
      if ((abs $anch_diff) < 300000) {
	push @diffs, $anch_diff;
      }
    }
  }

  # calculate the average of the difference to adjust the locations
  my $average = 0;
  foreach my $value (@diffs) {
    $average += ($value / @diffs);
  }

  my $newstart = int($cont_s + $average);
  my $newend = int($cont_e + $average);
  push @newlocs, "${newstart}_${newend}";
}

open OUTPUT, ">adjusted_locations.txt";

foreach my $chrom (1..24) {
  $newlocs[$chrom - 1] =~ m/(?<st>\d+)_(?<end>\d+)/;
  print OUTPUT "$chrom\t$+{st}\t$+{end}\n";
}

close OUTPUT;
