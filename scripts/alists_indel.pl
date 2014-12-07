#!usr/bin/perl

# Two files need to be generated in order to analyze the two strains of
# Prevotella Intermedia: one is a list of scores between protein-ORF pairs;
# the other is a list of indels and mismatches between the two strains.
# The scores list would contain the ids, total length, aligned length,
# unaligned lengths, identities, mismatches, and gaps.
# The indels list would contain the ids, locations of indels,
# and the nucleic acid/protein on each strain.
# This script would only generate the indels list.

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

# The location of PI 17 proteins should be read so that we can identify
# the location of each protein-ORF pair and their indels.

my %ref_st;
my %ref_end;
my %ref_order;
open REFL, "PI17.ptt";
while (<REFL>) {
  chomp;
  my $line = $_;
  if (/(?<st>\d+)\.\.(?<end>\d+)\t(?<order>.)\t\d+\t(?<rid>\d+)/) {
    $ref_st{$+{rid}} = $+{st};
    $ref_end{$+{rid}} = $+{end};
    $ref_order{$+{rid}} = $+{order};
  }
}
close REFL;

# The list would include locations, gene/ORF ids, respective protein/ORF
# amino acid. To access the indel information pair by pair, we need to
# specify how many indels a protein-ORF pair has.
my %pairs;
my %num_indels;
my %adj_loc;
my %ref_seq;
my %shot_seq;
my %ref_index;
my %shot_index;
my %indices;
my %locations;
my %indel_types;
my @mismatches;
my @gaps;
my %blocks;
my $shot_id;
my $ref_id;
my $shot_full;
my $result;
my $al_res;
my $flag = 0;

# We need to read the Needleman-Wunsch results, analyzing the relative locations
# of indels and changes.
open PROT, "amino.global.aligned.txt";
while (<PROT>) {
  chomp;
  my $line = $_;
  if ($flag == -1) {
    $flag = 0;
  } elsif (/>(?<sid>pint2554_c_\d+_\d+)/) {
    $shot_id = $+{sid};
  } elsif (/>gi\|(?<rid>\d+)/) {
    $ref_id = $+{rid};
    $pairs{$shot_id} = $ref_id;
  } elsif ($flag == 0) {
    $shot_full = $line;
    $flag = 1;
  } elsif ($flag == 1) {
    $result = $line;
    $al_res = $line;
    $al_res =~ s/(^ +)|( +$)//g;
    $adj_loc{$shot_id} = index($result, $al_res, 0);
    $flag = 2;
  } else {
    @mismatches = &find_index($al_res, '*');
    @gaps = &find_index($al_res, ' ');
    # gaps should be connected to blocks
    my $start = 0;
    my $len = 0;
    my $changed = 0;
    foreach my $key (@gaps) {
      if ($key > $start + $len) {
	if ($changed) {
	  $blocks{$start} = $len;
	}
	$start = $key;
	$len = 1;
	$changed = 1;
      } else {
	$len++;
	$changed = 1;
      }
    }
    if ($changed) {
      $blocks{$start} = $len;
    }
    $num_indels{$ref_id} = @mismatches + (keys %blocks);
    my %index_sorter;
    foreach my $ind (@mismatches) {
      $index_sorter{$ind} = 'm';
    }
    foreach my $ind (keys %blocks) {
      $index_sorter{$ind} = 'i';
    }
    my @sorted_indices = sort {$a <=> $b} (keys %index_sorter);
    foreach my $key (1..$num_indels{$ref_id}) {
      my $ptag = "${ref_id}_$key";
      $indices{$ptag} = $sorted_indices[$key - 1];
      $indel_types{$ptag} = $index_sorter{$sorted_indices[$key - 1]};
      if ($indel_types{$ptag} eq 'm') {
	$ref_seq{$ptag} =
	  substr ($line, ($indices{$ptag} + $adj_loc{$shot_id}), 1);
	$shot_seq{$ptag} = 
	  substr ($shot_full, ($indices{$ptag} + $adj_loc{$shot_id}), 1);
      } else {
	my $ref_indic = 
	  substr ($line, ($indices{$ptag} + $adj_loc{$shot_id}), 1);
	my $blen = $blocks{$sorted_indices[$key - 1]};
	if ($ref_indic eq '-') {
	  $ref_seq{$ptag} = $ref_indic;
	  $shot_seq{$ptag} = 
	    substr ($shot_full, ($indices{$ptag} + $adj_loc{$shot_id}), $blen);
	} else {
	  $ref_seq{$ptag} =
	    substr ($line, ($indices{$ptag} + $adj_loc{$shot_id}), $blen);
	  $shot_seq{$ptag} = '-';
	}
      }
      # finding the location of a particular indel requires knowing the 
      # blanks inserted by Needleman-Wunsch and translating the indices
      # to proper locations.
      my $front_part = substr $line, 0, ($indices{$ptag} + $adj_loc{$shot_id});
      my $blanks = scalar (&find_index($front_part, '-'));
      my $adjust = $adj_loc{$shot_id} - $blanks + $indices{$ptag};
      if ($ref_order{$ref_id} eq '+') {
	$locations{$ptag} = $ref_st{$ref_id} + (3 * $adjust);
      } else {
	$locations{$ptag} = $ref_end{$ref_id} - (3 * $adjust);
      }
    }
    %blocks = ();
    $flag = -1;
  }
}
close PROT;

open OUTPUT, ">PI.amino.indels.txt";
foreach my $sid (sort {$pairs{$a} cmp $pairs{$b}} (keys %pairs)) {
  my $rid = $pairs{$sid};
  print OUTPUT "$rid :: $sid\n";
  foreach my $num (1..$num_indels{$rid}) {
    my $pid = "${rid}_$num";
    print OUTPUT "$locations{$pid}\t$ref_seq{$pid}\t$shot_seq{$pid}\n";
  }
  print OUTPUT "//\n\n";
}
