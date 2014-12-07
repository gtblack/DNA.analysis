#!/usr/bin/perl

# Many tools accept only one FASTA file as input and then aligns the genomes
# included in the file. We need to parse the original files to put only
# the nucleotidal sequences interested in being parsed into one file.
# Before that, the original file should be parsed into separate files.

use strict;

open INPUT, $ARGV[0];

my $add = $ARGV[1];
my $count = 1;

while (<INPUT>) {
  my $line = $_;
  if (/^>/) {
    close OUTPUT;
    open OUTPUT, ">${add}seg.${count}.fna";
    $count++;
    print OUTPUT $line;
  } else {
    print OUTPUT $line;
  }
}

close OUTPUT;
close INPUT;
