#!usr/bin/perl

# Now that the supposed locations of each PI25611 contig is known, we need to
# find out only the best match protein-ORF pairs that fit in the locations.
# Since the ORFs are based on the nucleotide data of PI25611, we need
# to sort out the suggested proteins that are not in the reference range.

use strict;
use warnings;

my %pairs;
my %best_pairs;
my %un_pairs;
my %scores;
my %evalues;
my %align;
my %pidlocs;
my $query_tag;
my $score;
my $line_flag = 0;
my $chrom;
my $chrom2_add = 579647;
my $chrom2_first_pid = 387132090;

# read the alignment locations of PI25611 shotgun sequences
open ALIGN, "adjusted_locations.txt";
while (<ALIGN>) {
  chomp;
  my @line = split /\t/, $_;
  $align{$line[0]} = "$line[1]_$line[2]";
}
close ALIGN;

# read the locations of each annotated proteins in PI17
open PTT, "PI17.ptt";
while (<PTT>) {
  chomp;
  if (/(?<sloc>\d+)\.\.(?<eloc>\d+)\t.\t[^\s]+\t(?<reftag>\d+)/) {
    if ($+{reftag} < $chrom2_first_pid) {
      $pidlocs{$+{reftag}} = "$+{sloc}_$+{eloc}";
    } else {
      my $pidsloc = $+{sloc} + $chrom2_add;
      my $pideloc = $+{eloc} + $chrom2_add;
      $pidlocs{$+{reftag}} = "${pidsloc}_${pideloc}";
    }
  }
}

close PTT;

my $temp_match = "";
my $temp_score = 0;
my $temp_eval = 0;

open BLAST, "PI17.blastp.result.txt";
while (<BLAST>) {
  chomp;
  my $line = $_;

  # first read the query protein segment and store the tag
  if ($line =~ m/^Query= (?<qtag>\w+)/) {
    $query_tag = $+{qtag};
    $line_flag = 1;
    $query_tag =~ m/_c_(?<chr>\d+)/;
    $chrom = $+{chr};
  } elsif ($line_flag > 0 && $line =~ m/ gi\|(?<rtag>\d+)[^\s]+ .{33} +(?<scr>[^\s]+) +(?<ev>[^\s]+)/) {
    my $rtag = $+{rtag};
    my $scr = $+{scr};
    my $eval = $+{ev};
    $align{$chrom} =~ m/(?<conts>\d+)_(?<conte>\d+)/;
    my $cont_s = $+{conts};
    my $cont_e = $+{conte};
    $pidlocs{$rtag} =~ m/(?<prots>\d+)_(?<prote>\d+)/;
    my $prot_s = $+{prots};
    my $prot_e = $+{prote};

    # if the protein is in a proper range and the evalue is low enough
    # the value should be stored
    if ($cont_s < $prot_e && $prot_s < $cont_e && $eval < 0.1) {
      $pairs{$query_tag} = $rtag;
      $scores{$query_tag} = $scr;
      $evalues{$query_tag} = $eval;
      $best_pairs{$query_tag} = $rtag if $line_flag == 1;
      $line_flag = 0;
    # if the protein is not in proper range the next value should be checked
    # in case we do not find a good match we should store the best match
    } elsif ($eval < 0.1 && $line_flag == 1) {
      $temp_match = $rtag;
      $temp_score = $scr;
      $temp_eval = $eval;
      $line_flag = 2;
    # if the evalue is not low enough there is no longer a good match left.
    # save the temporary best match but tag it
    } elsif ($eval >= 0.1) {
      $un_pairs{$query_tag} = "$temp_match";
      $scores{$query_tag} = $temp_score;
      $evalues{$query_tag} = $temp_eval;
      $line_flag = 0;
    }
  }
}
close BLAST;

open OUTPUT, ">BLASTp_aligned_matches.txt";
print OUTPUT "PI25611\t\tPI17\t\tBit Score\tE value\n";
foreach my $key (sort keys %pairs) {
  print OUTPUT "${key}\t$pairs{$key}\t$scores{$key}\t$evalues{$key}\n";
}
close OUTPUT;

open OUTPUT, ">BLASTp_unaligned_matches.txt";
print OUTPUT "PI25611\t\tPI17\t\tBit Score\tE value\n";
foreach my $key (sort keys %un_pairs) {
  print OUTPUT "${key}\t$un_pairs{$key}\t$scores{$key}\t$evalues{$key}\n";
}
close OUTPUT;

open OUTPUT, ">BLASTp_best_matches.txt";
print OUTPUT "PI25611\tPI17\t\tBit Score\tE value\n";
foreach my $key (sort keys %best_pairs) {
  print OUTPUT "${key}\t$best_pairs{$key}\t$scores{$key}\t$evalues{$key}\n";
}
close OUTPUT;
