#!/usr/bin/perl

use strict;
use warnings;

my $file_contigs = 'PI25611.contigs.tab';
my $file_alignment = 'alignment3';
my $file_refptt = 'PI17.ptt';
my $file_queryfasta = 'PI_ATCC_25611_DSM_20706.aa.fasta';
my $file_blast = 'PI17.blastp.result.txt';
my $file_refnucseq = 'PI17.ffn';
my $file_refamnseq = 'PI17.faa';
my $file_querynucseq = 'Prevotella_Intermedia_ATCC25611_DSM20706.na.fasta';
my $file_queryamnseq = 'PI_ATCC_25611_DSM_20706.aa.fasta';

# NEEDS TO BE FIXED: CONSTANTS USED FOR PI17 may be different in other strains
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
  # this array will store the anchor location differences that will later
  # be summed up for contig readjustment.
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
      # the constant 300000 is the maximum distance allowed between a pair
      # of anchors; otherwise they will be considered to be on different
      # contigs and be disregarded.
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

  # the new locations of the contigs should be stored separately
  # in order for original contig locations to be utilized in the case of
  # anchor location adjustment in class II classification process.
  my $newstart = int($cont_start + $average);
  my $newend = int ($cont_end + $average);
  $new_contig{$chrom} = [$newstart, $newend];
}

# reference_locs stores the location of each protein as value, while
# reference_tags stores each protein id as value and its location as key.
# The first hash is needed for getting locations to judge whether a certain
# pair fits in the alignment, while the second is needed for guessing the
# protein id based on locations (used for tagging fasta sequences).
my %reference_locs;
my %reference_tags;

# read the location of each annotated proteins in the reference genome
open PTT, $file_refptt;
while (<PTT>) {
  chomp;
  # NEEDS TO BE FIXED: DEPENDENT ON SPECIFIC PTT FORMAT
  if (/(?<sloc>\d+)\.\.(?<eloc>\d+)\t(?<order>.)\t[^\s]+\t(?<reftag>\d+)/) {
    if ($+{reftag} < $chrom2_first_pid) {
      $reference_locs{$+{reftag}} = [$+{sloc}, $+{eloc}];
      if ($+{order} eq '+') {
	$reference_tags{"1_$+{sloc}-$+{eloc}"} = $+{reftag};
      } else {
	$reference_tags{"1_c$+{eloc}-$+{sloc}"} = $+{reftag};
      }
    } else {
      my $ridsloc = $+{sloc} + $chrom2_add;
      my $rideloc = $+{eloc} + $chrom2_add;
      $reference_locs{$+{reftag}} = [$ridsloc, $rideloc];
      if ($+{order} eq '+') {
	$reference_tags{"2_$+{sloc}-$+{eloc}"} = $+{reftag};
      } else {
	$reference_tags{"2_c$+{eloc}-$+{sloc}"} = $+{reftag};
      }
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
  # NEEDS TO BE FIXED: FINDING CONTIG DEPENDANT ON PROTEIN/ORF NAMING CONVENTION
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
    # the value should be stored as class I.

    if ($contig_start < $protein_end && $protein_start < $contig_end && $blast_evalue < 0.1) {
      $class1_pairs{$query_tag} = $reference_tag;
      $class1_scores{$query_tag} = $blast_score;
      $class1_evalues{$query_tag} = $blast_evalue;
      $flag = 0;
    } else {
      # each anchor pair should be checked to see whether a certain anchor
      # pair matches with the protein-ORF pair.
      # Given an anchor pair and a query-reference protein/ORF pair,
      # the reference protein should be within range of the reference anchor
      # and the query protein within that of the query anchor.
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

# With the pairs divided into classes, we need nucleotidal/amino sequences
# from fasta files for global alignment. The sequences need to be in the
# correct order so that the Needleman-Wunsch algorithm can see which
# sequence to compare against.

my %sequence_querynucl;
my %sequence_queryamino;
my %sequence_refnucl;
my %sequence_refamino;

my $temp_tag;
my $temp_line;

open QUERYNUCL, $file_querynucseq;
while (<QUERYNUCL>) {
  chomp;
  my $line = $_;
  if (/^>(?<stag>[^\s]+) /) {
    $temp_tag = $+{stag};
    $temp_line = $line . "\n";
  } else {
    $temp_line .= ($line . "\n");
    $sequence_querynucl{$temp_tag} = $temp_line;
  }
}
close QUERYNUCL;

$temp_tag = '';
$temp_line = '';

open QUERYAMINO, $file_queryamnseq;
while (<QUERYAMINO>) {
  chomp;
  my $line = $_;
  if (/^>(?<stag>[^\s]+) /) {
    $temp_tag = $+{stag};
    $temp_line = $line . "\n";
  } else {
    chop $line;
    $temp_line .= ($line . "\n");
    $sequence_queryamino{$temp_tag} = $temp_line;
  }
}
close QUERYAMINO;

$temp_tag = '';
$temp_line = '';

# the fasta file of the reference nucleotidal sequence did not contain
# proper protein ids, so the id had to be manually searched based on
# nucleotidal locations.
open REFNUCL, $file_refnucseq;
while (<REFNUCL>) {
  chomp;
  my $line = $_;
  if (/>g.+\|NC_01786(?<chrom>.)\.1\|:(?<loc>[^\s]+) /) {
    $sequence_refnucl{$temp_tag} = ($temp_line . "\n");
    my $chrom = $+{chrom} + 1;
    $temp_tag = $reference_tags{"${chrom}_$+{loc}"};
    $temp_line = ">gi\|${temp_tag}|:$+{loc}\n";
  } else {
    $temp_line .= $line;
  }
}
$sequence_refnucl{$temp_tag} = $temp_line;
close REFNUCL;

$temp_tag = '';
$temp_line = '';

open REFAMINO, $file_refamnseq;
while (<REFAMINO>) {
  chomp;
  my $line = $_;
  if (/>gi\|(?<rtag>\d+)\|/) {
    $sequence_refamino{$temp_tag} = ($temp_line . "\n");
    $temp_tag = $+{rtag};
    $temp_line = $line . "\n";
  } else {
    $temp_line .= $line;
  }
}
$sequence_refamino{$temp_tag} = $temp_line;
close REFAMINO;

open OUTREFAMINO, ">class1.reference.ordered.aa.fasta";
open OUTREFNUCL, ">class1.reference.ordered.na.fasta";
open OUTQUERYAMINO, ">class1.query.ordered.aa.fasta";
open OUTQUERYNUCL, ">class1.query.ordered.na.fasta";
foreach my $tag (sort keys %class1_pairs) {
  print OUTREFAMINO $sequence_refamino{$class1_pairs{$tag}};
  print OUTREFNUCL $sequence_refnucl{$class1_pairs{$tag}};
  print OUTQUERYAMINO $sequence_queryamino{$tag};
  print OUTQUERYNUCL $sequence_querynucl{$tag};
}
close OUTREFAMINO;
close OUTREFNUCL;
close OUTQUERYAMINO;
close OUTQUERYNUCL;

open OUTREFAMINO, ">class2.reference.ordered.aa.fasta";
open OUTREFNUCL, ">class2.reference.ordered.na.fasta";
open OUTQUERYAMINO, ">class2.query.ordered.aa.fasta";
open OUTQUERYNUCL, ">class2.query.ordered.na.fasta";
foreach my $tag (sort keys %class2_pairs) {
  print OUTREFAMINO $sequence_refamino{$class2_pairs{$tag}};
  print OUTREFNUCL $sequence_refnucl{$class2_pairs{$tag}};
  print OUTQUERYAMINO $sequence_queryamino{$tag};
  print OUTQUERYNUCL $sequence_querynucl{$tag};
}
close OUTREFAMINO;
close OUTREFNUCL;
close OUTQUERYAMINO;
close OUTQUERYNUCL;

open OUTREFAMINO, ">nonmatching.reference.ordered.aa.fasta";
open OUTREFNUCL, ">nonmatching.reference.ordered.na.fasta";
open OUTQUERYAMINO, ">nonmatching.query.ordered.aa.fasta";
open OUTQUERYNUCL, ">nonmatching.query.ordered.na.fasta";
foreach my $tag (sort keys %nonmatching_pairs) {
  print OUTREFAMINO $sequence_refamino{$nonmatching_pairs{$tag}};
  print OUTREFNUCL $sequence_refnucl{$nonmatching_pairs{$tag}};
  print OUTQUERYAMINO $sequence_queryamino{$tag};
  print OUTQUERYNUCL $sequence_querynucl{$tag};
}
close OUTREFAMINO;
close OUTREFNUCL;
close OUTQUERYAMINO;
close OUTQUERYNUCL;
