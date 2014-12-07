#!usr/bin/perl

# The result of the BLASTP analysis needs to be reorganized so that
# the matching protein segments between the two amino acid database
# can be used to align the shotgun sequences from PL25611.
# We should first find out whether the matching sequence is trustworthy
# enough to be used for alignment. That should be figured out using the
# E-value and Identities percentage.
# Then, once the matching segment is good enough, we find the segment
# in the ffn file of PI17 (reference data) and find out the exact
# location using the tag given from the BLASTP output.
# The pair from PI25611, the shotgun sequences, will be looked up
# from the FASTA file from HOMD. Matching the positions between PI17
# and PI25611 will be a simple subtraction from there.
# We should sum up all the scores from the protein sequences based on
# credibility so that we can judge where the 24 shotgun sequences can
# be placed.

use strict;

my %pairs;
my %scores;
my %bscores;
my %evalues;
my $query_tag;
my $score;
my $line_flag = 0;

open BLAST, "PI17.blastp.result.txt";

while (<BLAST>) {
  chomp;
  my $line = $_;

  # first read the query protein segment and store the tag
  if ($line =~ m/^Query= (?<qtag>\w+)/) {
    
    # there can be cases where there is only one segment that can be matched.
    # in this case we need to store the score and reset the flag.

    if ($line_flag eq 3) {
      $scores{$query_tag} = $score;
      $line_flag = 0;
    }
    $line_flag++;
    $query_tag = $+{qtag};

  # then find the highest ranked reference protein segment and store its tag
  } elsif (($line =~ m/ gi\|(?<rtag>[^\|]+)[^\s]+ .{33} +(?<scr>[^\s]+) +(?<ev>[^\s]+)/) && ($line_flag eq 1)) {
    $line_flag++;
    $pairs{$query_tag} = $+{rtag};
    $bscores{$query_tag} = $+{scr};
    $evalues{$query_tag} = $+{ev};

    # find the Identities score of the highest ranked segment
  } elsif (($line_flag eq 2) && ($line =~ m/Identities/)) {
    my @pline = split / /, $line;
    $score = substr $pline[4], 1, 2;
    $line_flag++;

  # find the Identities score of the second highest ranked segment
  } elsif (($line_flag eq 3) && ($line =~ m/Identities/)) {
    my @pline = split / /, $line;
    my $div = substr $pline[4], 1, 2;
    $score = ($score ** 2) / $div;
    $scores{$query_tag} = $score;
    $line_flag = 0;
  } else {
    next;
  }
}

close BLAST;

my %refloc;

open PTT, "PI17.ptt";

while (<PTT>) {
  if (/(?<loc>[^\s]+)\t(?<order>.)\t[^\s]+\t(?<reftag>\d+)/) {
    if ($+{order} eq '+') {
      $refloc{$+{reftag}} = $+{loc};
    } else {
      my @loc_temp = split /\.\./, $+{loc};
      $refloc{$+{reftag}} = "$loc_temp[1]..$loc_temp[0]";
    }
  }
}

close PTT;

my %querysloc;
my %queryeloc;

open FAA, "PI_ATCC_25611_DSM_20706.aa.fasta";

while (<FAA>) {
  if (/>(?<tag>[^\s]+) [^\[]+\[(?<fdig>\d+) - (?<sdig>\d+)/) {
    $querysloc{$+{tag}} = $+{fdig};
    $queryeloc{$+{tag}} = $+{sdig};
  }
}

close FAA;

# Checking in the middle to see whether we got the right data

=begin testcode
open OUTPUT, ">test_output.txt";
print OUTPUT "Query\tDatabase\tScore\tQuery Location\tData Start\tData End\n";

foreach my $key (keys %pairs) {
  print OUTPUT "${key}\t$pairs{$key}\t$scores{$key}\t$refloc{$pairs{$key}}\t";
  print OUTPUT "$querysloc{$key}\t$queryeloc{$key}\n";
}

close OUTPUT;
=end testcode
=cut

=begin parseout
open OUTPUT, ">BLASTp_parsed.txt";
print OUTPUT "Query\t\t\tDatabase\tScore\tE Value\n";

foreach my $key (keys %pairs) {
  if (length $key > 15) {
    print OUTPUT "${key}\t$pairs{$key}\t$bscores{$key}\t$evalues{$key}\n";
  } else {
    print OUTPUT "${key}\t\t$pairs{$key}\t$bscores{$key}\t$evalues{$key}\n";
  }
}

close OUTPUT;
=end parseout
=cut

# We have to calculate the relative positions between the shotgun sequences
# and the reference sequences. There will be two hashes that store the
# information about the positions; one will store possible positions for each
# shotgun sequence, with the key as 1, 2, ..., 23, 24 and the value as
# 3-1-14015, 7-2-498, 21-1-48012, and so on; the other will store the score
# of each position, with the key from the value of the former hash
# 3-1-14015, ..., and the value as the score.

my %pscores;

foreach my $key (keys %pairs) {
  $key =~ m/[^_]+_c_(?<chr>\d+)/;
  my $shotchr = $+{chr};
  my $refchr = ($pairs{$key} < 387132090) ? 1 : 2;
  $refloc{$pairs{$key}} =~ m/(?<sloc>\d+)\.\.(?<eloc>\d+)/;
  my $refmin = ($+{sloc} > $+{eloc}) ? $+{eloc} : $+{sloc};
  my $hyppos = ($querysloc{$key} > $queryeloc{$key}) ?
    $refmin - $queryeloc{$key} + 1 : $refmin - $querysloc{$key} + 1;
  my $postag = $shotchr . "_" . $refchr . "_" . $hyppos;
  $pscores{$postag} = $pscores{$postag} + $scores{$key};
}

open OUTPUT, ">alignment_scores.txt";
print OUTPUT "PI25611\tPI17\tPosition\tScore\n";

foreach my $key (1..24) {
  foreach my $poskey (sort {$pscores{$b} <=> $pscores{$a}} (keys %pscores)) {
    $poskey =~ m/(?<shot>\d+)_(?<ref>\d+)_(?<pos>.+)/;
    if ($key eq $+{shot}) {
      print OUTPUT "$key\t$+{ref}\t$+{pos}\t\t$pscores{$poskey}\n";
    }
  }
  print OUTPUT "\n";
}

close OUTPUT;
