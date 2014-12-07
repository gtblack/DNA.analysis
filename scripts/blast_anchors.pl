#!/usr/bin/perl

use strict;
use warnings;

my $chrom2_add = 579647;
my $chrom2_first_pid = 387132090;

my %align_start;
my $flag = 0;

sub ptrue {
  if ($_[0]) {
    return "1";
  } else {
    return "0";
  }
}

open ALOC, "PI25611.contigs.tab";
while (<ALOC>) {
  chomp;
  my $line = $_;
  my @info = split /\t/, $line;
  if (/Ordered Contigs/) {
    $flag = 1;
  } elsif (($flag == 1) && ($line eq "")) {
    $flag = 0;
  } elsif (($flag == 1) && ($info[0] eq "contig")) {
    $align_start{$info[1]} = $info[4];
  }
}
close ALOC;

my %anchors;
my $temp_key;

open ALIGN, "alignment3";
while (<ALIGN>) {
  chomp;
  my $line = $_;
  if (($line =~ /> 1:(?<key>[^\s]+)/) && ($flag == 0)) {
    $temp_key = $+{key};
    $flag = 1;
  } elsif (($line =~ /> 2:(?<key>[^\s]+)/) && ($flag == 1)) {
    $anchors{$temp_key} = $+{key};
    $flag = 0;
  } elsif (($line =~ /> 1:/) && ($flag == 1)) {
    $flag = -1;
  }
}
close ALIGN;

my %ref_st;
my %ref_end;

open REFL, "PI17.ptt";
while (<REFL>) {
  chomp;
  my @tokens = split /\t/, $_;
  if (($tokens[0] =~ /(?<st>\d+)\.\.(?<end>\d+)/) && (defined $tokens[1])) {
    $ref_st{$tokens[3]} = $+{st};
    $ref_end{$tokens[3]} = $+{end};
  }
}
close REFL;

my %shot_st;
my %shot_end;

open SHOTL, "PI_ATCC_25611_DSM_20706.aa.fasta";
while (<SHOTL>) {
  chomp;
  my $line = $_;
  if (/>(?<sid>[^\s]+)[^\[]+\[(?<st>\d+) - (?<end>\d+)/) {
    if ($+{st} < $+{end}) {
      $shot_st{$+{sid}} = $+{st};
      $shot_end{$+{sid}} = $+{end};
    } else {
      $shot_st{$+{sid}} = $+{end};
      $shot_end{$+{sid}} = $+{st};
    }
  }
}
close SHOTL;

my %pairs;
my %scores;
my %evalues;
my %bpairs;
my $qr_tag;
my $qr_chrom;

open BLAST, "PI17.blastp.result.txt";
while (<BLAST>) {
  chomp;
  my $line = $_;
  if ($line =~ /^Query= (?<qtag>\w+)/) {
    $qr_tag = $+{qtag};
    $flag = 1;
    $qr_tag =~ /(?<chr>pint2554_c_\d+)/;
    $qr_chrom = $+{chr};
  } elsif ($flag > 0 && $line =~ / gi\|(?<rtag>\d+)[^\s]+ .{33} +(?<scr>[^\s]+) +(?<ev>[^\s]+)/) {
    my $rtag = $+{rtag};
    my $scr = $+{scr};
    my $eval = $+{ev};
    foreach my $key (keys %anchors) {
      $anchors{$key} =~ /(?<shst>\d+)-(?<shend>\d+)/;
      my $anc_sst = $+{shst};
      my $anc_send = $+{shend};
      $key =~ /(?<rst>\d+)-(?<rend>\d+)/;
      my $anc_rst = $+{rst};
      my $anc_rend = $+{rend};
      my $cond_chr = $rtag > $chrom2_first_pid;
      my $cond_ss = $anc_sst < $align_start{$qr_chrom} + $shot_end{$qr_tag};
      my $cond_se = $anc_send > $align_start{$qr_chrom} + $shot_st{$qr_tag};
      my $cond_r1 = ($anc_rst < $ref_end{$rtag}) &&
	($anc_rend > $ref_st{$rtag});
      my $cond_r2 = ($anc_rst < ($ref_end{$rtag} + $chrom2_add)) &&
	($anc_rend > ($ref_st{$rtag} + $chrom2_add));
      if ($cond_ss && $cond_se && (($cond_r1 && !$cond_chr) || ($cond_r2 && $cond_chr))) {
	$pairs{$qr_tag} = $rtag;
	$scores{$qr_tag} = $scr;
	$evalues{$qr_tag} = $eval;
	if ($flag == 1) {
	  $bpairs{$qr_tag} = $rtag;
	}
	$flag = 0;
	last;
      }
    }
    $flag++;
  }
}
close BLAST;

open OUTPUT, ">BLASTp_mauve_matches.txt";
print OUTPUT "PI25611\t\tPI17\t\tBit Score\tE value\n";
foreach my $key (sort {$pairs{$a} cmp $pairs{$b}} (keys %pairs)) {
  print OUTPUT "$key\t$pairs{$key}\t$scores{$key}\t$evalues{$key}\n";
}
close OUTPUT;

open OUTPUT, ">BLASTp_mauve_best.txt";
print OUTPUT "PI25611\t\tPI17\t\tBit Score\tE value\n";
foreach my $key (sort {$bpairs{$a} cmp $bpairs{$b}} (keys %bpairs)) {
  if ($evalues{$key} < 0.1 ) {
    print OUTPUT "$key\t$bpairs{$key}\t$scores{$key}\t$evalues{$key}\n";
  }
}
close OUTPUT;
