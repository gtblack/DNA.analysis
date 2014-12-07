#!usr/bin/perl

# Two files need to be generated in order to analyze the two strains of
# Prevotella Intermedia: one is a list of scores between protein-ORF pairs;
# the other is a list of indels and mismatches between the two strains.
# The scores list would contain the ids, total length, aligned length,
# unaligned length, identities, mismatches, and gaps.
# The indels list would contain the ids, location of indels,
# and the nucleic acid on each strain.
# This script would only generate the scores list.

use strict;
use warnings;

# finding all indices of a substring of string
sub find_index {
  my ($str, $ptrn) = @_;
  my $offset = 0;
  my $pos;
  my @ret;
  while (1) {
    $pos = index($str, $ptrn, $offset);
    last if ($pos < 0);
    push @ret, $pos;
    $offset = $pos + 1;
  }
  return @ret;
}
  
# The ids of ORF-protein pairs should be stored for later reference.
# Locations of each sequence should be stored to pinpoint the location
# of a certain indel.
# Length will be calculated using the raw data.
# The aligned portion and the unaligned portion will be separated so that
# the length of each portion can be calculated.
my %pairs;
my $shot_id;
my $ref_id;
my %shot_lens;
my %ref_lens;
my $flag = 0;
my $shot_orf;
my $align_result;
my %alens;
my %shot_unlens;
my %ref_unlens;
my %iden_counts;
my %mis_counts;
my %gap_counts_shot;
my %gap_counts_ref;
my @mismatches;
my @gaps;

open PROT, "amino.global.aligned.txt";
while (<PROT>) {
  chomp;
  my $line = $_;
  if ($flag == -1) {
    $flag = 0;
  # the first line marks PI25611 and its id and location
  } elsif (/>(?<sid>pint\d+_c_\d+_\d+) [^\[]+\[(?<sloc>\d+) - (?<eloc>\d+)/) {
    $shot_id = $+{sid};
  # the second line marks PI17 and its id
  } elsif (/>gi\|(?<refid>\d+)/) {
    $ref_id = $+{refid};
    $pairs{$shot_id} = $ref_id;
  # the third line contains the contig sequences,
  # which should be stored for later use
  } elsif ($flag == 0) {
    $shot_orf = $line;
    $shot_lens{$shot_id} = length $line;
    $flag = 1;
  # the fourth line contains the result of Needleman-Wunsch
  } elsif ($flag == 1) {
    # protein length cannot be measured by nucleotidal locations
    # mismatch and gap indices will later be used for listing indels
    @mismatches = &find_index ($line, '*');
    @gaps = &find_index ($line, ' ');
    $flag = 2;
    $align_result = $line;
    my $align = $line;
    $align =~ s/^ +| +$//g;
    $alens{$shot_id} = length $align;
    $shot_unlens{$shot_id} = $shot_lens{$shot_id} - length($align);
    $iden_counts{$shot_id} = ($align =~ s/\|/\|/g) || 0;
    $mis_counts{$shot_id} = ($align =~ s/\*/\*/g) || 0;
  # the fifth line contains the reference sequences, which is needed to
  # complete the indel lists.
  } else {
    $ref_lens{$ref_id} = length $line;
    $ref_unlens{$ref_id} = $ref_lens{$ref_id} - $alens{$shot_id};
    $align_result =~ m/(^ +)/;
    my $front_sp = defined $1 ? length $1 : 0;
    $align_result =~ m/( +$)/;
    my $tail_sp = defined $1 ? length $1 : 0;
    $gap_counts_shot{$shot_id} = 0;
    $gap_counts_ref{$shot_id} = 0;
    foreach my $gap (@gaps) {
      if ($gap >= $front_sp && ($gap < ((length $align_result) - $tail_sp))) {
	if (substr ($shot_orf, $gap, 1) eq '-') {
	  $gap_counts_shot{$shot_id}++;
	} else {
	  $gap_counts_ref{$shot_id}++;
	}
      }
    }
    $flag = -1;
  }
}
close PROT;

open PROUT, ">PI.amino.scores.txt";
print PROUT "PI17\tPI25611\tPI17 Length\tPI25611 Length\tAligned\t";
print PROUT "PI17 Nonaligned\tPI25611 Nonaligned\tIdentities\t";
print PROUT "Mismatches\tPI17 Gaps\tPI25611 Gaps\n";
foreach my $sid (sort keys %pairs) {
  my $rkey = $pairs{$sid};
  print PROUT "$sid\t$rkey\t$shot_lens{$sid}\t$ref_lens{$rkey}\t";
  print PROUT "$alens{$sid}\t$ref_unlens{$rkey}\t$shot_unlens{$sid}\t";
  my $iden_per = ($iden_counts{$sid} / $alens{$sid} * 100);
  my $mis_per = ($mis_counts{$sid} / $alens{$sid} * 100);
  my $gap_shotper = ($gap_counts_shot{$sid} / $alens{$sid} * 100);
  my $gap_refper = ($gap_counts_ref{$sid} / $alens{$sid} * 100);
  my $format = "%.2f\t%.2f\t%.2f\t%.2f\n";
  printf PROUT $format, $iden_per, $mis_per, $gap_refper, $gap_shotper;
}
close PROUT;
