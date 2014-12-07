#!/usr/bin/perl

#

use strict;
use warnings;

my $file_contigs = 'PI25611.contigs.tab';
my $file_alignment = 'alignment3';
my $file_refptt = 'PI17.ptt';
my $file_queryfasta = 'PI_ATCC_25611_DSM_20706.aa.fasta';
my $file_blast = 'PI17.blastp.result.txt';

# NEEDS TO BE FIXED: CONSTANTS USED FOR PI17
my $chrom2_add = 579647;
my $chrom2_first_pid = 387132090;

open CONTIGS, $file_contigs;

my %contig_locs;
my %new_contig;
my $flag = 0;

# get the starting/ending locations of each contig relative to the reference
while (<CONTIGS>) {
  chomp;
  my $line = $_;
  my @block = split /\t/, $line;
  if (/Ordered Contigs/) {
    $flag = 1;
  } elsif (($flag == 1) && ($line eq "")) {
    $flag = 0;
  } elsif (($flag == 1) && ($block[0] eq "contig")) {
    $contig_locs{$block[1]} = [$block[4], $block[5]];
  }
}

$flag = 0;
close CONTIGS;


open ALIGN, $file_alignment;

my %anchor_pairs;
my $temp_key;

# get each anchor's position in the reference and shotgun sequences
while (<ALIGN>) {
  chomp;
  my $line = $_;
  # set the position in the reference sequence as the key
  if (($line =~ /> 1:(?<key>[^\s]+)/) && ($flag == 0)) {
    $temp_key = $+{key};
    $flag = 1;
  # set the position in the query sequences as the value
  } elsif (($line =~ /> 2:(?<st>\d+)-(?<end>\d+)/) && ($flag == 1)) {
    $anchor_pairs{$temp_key} = [$+{st}, $+{end}];
    $flag = 0;
  # non-matched anchors are discarded
  } elsif (($line =~ /> 1:/) && ($flag == 1)) {
    $flag = 2;
  }
}

$flag = 0;
close ALIGN;

foreach my $chrom (keys %contig_locs) {
  my @differences;

  # get the relative location of each shotgun sequence
  my $cont_start = $contig_locs{$chrom}[0];
  my $cont_end = $contig_locs{$chrom}[1];

  foreach my $key (keys %anchor_pairs) {
    my $anch_start = $anchor_pairs{$key}[0];
    my $anch_end = $anchor_pairs{$key}[1];

    # the anchor should be in the range of the contig and the differences
    # between the shotgun should be not too far away
    if (($cont_start <= $anch_end) && ($cont_end >= $anch_start)) {
      $key =~ m/(?<ref_st>\d+)-/;
      my $anch_diff = $+{ref_st} - $anch_start;
      if ((abs $anch_diff) < 300000) {
	push @differences, $anch_diff;
      }
    }
  }

  # calculate the average of the difference to adjust the locations
  my $average = 0;
  foreach my $value (@differences) {
    $average += ($value / @differences);
  }

  my $newstart = int($cont_start + $average);
  my $newend = int ($cont_end + $average);
  $new_contig{$chrom} = [$newstart, $newend];
}

my %reference_locs;

# read the location of each annotated proteins in the reference genome
open PTT, $file_refptt;
while (<PTT>) {
  chomp;
  # This part is based on the assumption that this chromosome 
  # NEEDS TO BE FIXED: DEPENDENT ON SPECIFIC PTT FORMAT
  if (/(?<sloc>\d+)\.\.(?<eloc>\d+)\t.\t[^\s]+\t(?<reftag>\d+)/) {
    if ($+{reftag} < $chrom2_first_pid) {
      $reference_locs{$+{reftag}} = [$+{sloc}, $+{eloc}];
    } else {
      my $ridsloc = $+{sloc} + $chrom2_add;
      my $rideloc = $+{eloc} + $chrom2_add;
      $reference_locs{$+{reftag}} = [$ridsloc, $rideloc];
    }
  }
}
close PTT;

my %query_locs;

# read location of each ORF on query strain
# NEEDS TO BE FIXED: DEPENDANT ON FASTA HEADER FORMAT
open QUERYFASTA, $file_queryfasta;
while (<QUERYFASTA>) {
  chomp;
  my $line = $_;
  if (/>(?<sid>[^\s]+)[^\[]+\[(?<st>\d+) - (?<end>\d+)/) {
    if ($+{st} < $+{end}) {
      $query_locs{$+{sid}} = [$+{st}, $+{end}];
    } else {
      $query_locs{$+{sid}} = [$+{end}, $+{st}];
    }
  }
}
close QUERYFASTA;

my $temp_match = "";
my $temp_score = 0;
my $temp_eval = 0;
my $query_tag;
my $chrom;

my %class1_pairs;
my %class1_scores;
my %class1_evalues;
my %class2_pairs;
my %class2_scores;
my %class2_evalues;
my %nonmatching_pairs;
my %nonmatching_scores;
my %nonmatching_evalues;

open BLAST, $file_blast;
while (<BLAST>) {
  chomp;
  my $line = $_;

  # first read the query protein segment and store the tag
  # NEEDS TO BE FIXED: FINDING CONIG DEPENDANT ON PROTEIN/ORF NAMING CONVENTION
  if ($line =~ m/^Query= (?<qtag>\w+)/) {
    $query_tag = $+{qtag};
    $flag = 1;
    $query_tag =~ m/(?<chr>[^\s]+)_\d+/;
    $chrom = $+{chr};
  } elsif (($flag > 0) && ($line =~ m/ gi\|(?<rtag>\d+)[^\s]+ .{33} +(?<scr>[^\s]+) +(?<ev>[^\s]+)/)) {
    my $reference_tag = $+{rtag};
    my $blast_score = $+{scr};
    my $blast_evalue = $+{ev};
    my $contig_start = $new_contig{$chrom}[0];
    my $contig_end = $new_contig{$chrom}[1];
    my $protein_start = $reference_locs{$reference_tag}[0];
    my $protein_end = $reference_locs{$reference_tag}[1];

    # if the protein is in a proper range and the E value is low enough
    # the value should be stored.

    if ($contig_start < $protein_end && $protein_start < $contig_end && $blast_evalue < 0.1) {
      $class1_pairs{$query_tag} = $reference_tag;
      $class1_scores{$query_tag} = $blast_score;
      $class1_evalues{$query_tag} = $blast_evalue;
      $flag = 0;
    } else {
      # each anchor pair should be checked to see whether a certain anchor
      # pair matches with the protein-ORF pair.
      # NEEDS TO BE FIXED: FINDING PAIRS DEPENDENT ON CHROMOSOME NUMBERS
      foreach my $key (keys %anchor_pairs) {
	my $anchor_query_start = $anchor_pairs{$key}[0];
	my $anchor_query_end = $anchor_pairs{$key}[1];
	$key =~ /(?<rst>\d+)-(?<rend>\d+)/;
	my $anchor_reference_start = $+{rst};
	my $anchor_reference_end = $+{rend};
	my $cond_ss = $anchor_query_start <
	  ($contig_locs{$chrom}[0] + $query_locs{$query_tag}[1]);
	my $cond_se = $anchor_query_end >
	  ($contig_locs{$chrom}[0] + $query_locs{$query_tag}[0]);
	my $cond_rs = $anchor_reference_start < $protein_end;
	my $cond_re = $anchor_reference_end > $protein_start;
	if ($cond_ss && $cond_se && $cond_rs && $cond_re && ($blast_evalue < 0.1)) {
	  $class2_pairs{$query_tag} = $reference_tag;
	  $class2_scores{$query_tag} = $blast_score;
	  $class2_evalues{$query_tag} = $blast_evalue;
	  $flag = 0;
	  last;
	}
      }
      if ($flag > 0) {
	$nonmatching_pairs{$query_tag} = $reference_tag;
	$nonmatching_scores{$query_tag} = $blast_score;
	$nonmatching_evalues{$query_tag} = $blast_evalue;
	$flag = 0;
      }
    }
  }
}

$flag = 0;
close BLAST;

open OUTPUT, ">class1.parsed.result.txt";
print OUTPUT "Query\t\tReference\tBit Score\tE value\n";
foreach my $key (sort {$class1_pairs{$a} cmp $class1_pairs{$b}} (keys %class1_pairs)) {
  print OUTPUT "$key\t$class1_pairs{$key}\t$class1_scores{$key}\t$class1_evalues{$key}\n";
}
close OUTPUT;

open OUTPUT, ">class2.parsed.result.txt";
print OUTPUT "Query\t\tReference\tBit Score\tE value\n";
foreach my $key (sort {$class2_pairs{$a} cmp $class2_pairs{$b}} (keys %class2_pairs)) {
  print OUTPUT "$key\t$class2_pairs{$key}\t$class2_scores{$key}\t$class2_evalues{$key}\n";
}
close OUTPUT;

open OUTPUT, ">nonmatching.parsed.result.txt";
print OUTPUT "Query\t\tReference\tBit Score\tE value\n";
foreach my $key (sort {$nonmatching_pairs{$a} cmp $nonmatching_pairs{$b}} (keys %nonmatching_pairs)) {
  print OUTPUT "$key\t$nonmatching_pairs{$key}\t$nonmatching_scores{$key}\t$nonmatching_evalues{$key}\n";
}
close OUTPUT;
