#!/usr/bin/perl

use strict;
use warnings;

open NUCL, "PI17.ffn";

my %pi_nucl;
my %pi_st;
my %pi_end;
my $temp_pid;
my $temp_st;
my $temp_end;
my $temp_line = "";
my $flag = 0;

while (<NUCL>) {
  chomp;
  my $line = $_;
  if (/>gi\|387131618[^:]+:(?<loc>[^ ]+)/) {
    my $loc = $+{loc};
    $loc =~ /(?<flc>\d+)-(?<slc>\d+)/;
    my $floc = $+{flc};
    my $sloc = $+{slc};
    if ($temp_line ne '') {
      $pi_nucl{$temp_pid} = $temp_line;
      $pi_st{$temp_pid} = $temp_st;
      $pi_end{$temp_pid} = $temp_end;
    }
    $temp_line = '';
    $flag = 0;
    if (($floc > 8647) && ($sloc < 18300)) {
      $temp_pid = $loc;
      if ($loc =~ /c/) {
	$temp_st = $sloc;
	$temp_end = $floc;
      } else {
	$temp_st = $floc;
	$temp_end = $sloc;
      }
      $flag = 1;
    }
  } elsif ($flag == 1) {
    $temp_line .= $line;
  }
}
close NUCL;

open VCF, "../raw/20141103/ftp.bioftp.org/20141103_WGRS/04_variants_vcf_file/TN1410D3298.dedup.AddRG.realn.vcf";

my %indels;
my %indel_pairs;
while (<VCF>) {
  chomp;
  my @line = split /\t/, $_;
  if ($line[0] eq 'I') {
    foreach my $key (keys %pi_nucl) {
      if (($line[1] > $pi_st{$key}) && ($line[1] < $pi_end{$key})) {
	$indels{$line[1]} = [$line[3], $line[4]];
	$indel_pairs{$line[1]} = $key;
      }
    }
  }
}
close VCF;

my %pi_result;
$temp_pid = '';
$temp_line = '';
my $adj_total = 0;

foreach my $key (sort {$a <=> $b} (keys %indels)) {
  my $adj_index = (length $indels{$key}[1]) - (length $indels{$key}[0]);
  if ($temp_pid eq '') {
    $temp_pid = $indel_pairs{$key};
    $temp_line = $pi_nucl{$temp_pid};
    my $ind = $key - $pi_st{$temp_pid};
    substr ($temp_line, $ind, (length $indels{$key}[0]), $indels{$key}[1]);
    $adj_total += $adj_index;
  } elsif ($temp_pid eq $indel_pairs{$key}) {
    my $ind = $key - $pi_st{$temp_pid} + $adj_total;
    substr ($temp_line, $ind, (length $indels{$key}[0]), $indels{$key}[1]);
    $adj_total += $adj_index;
  } else {
    $pi_result{$temp_pid} = $temp_line;
    $temp_pid = $indel_pairs{$key};
    $temp_line = $pi_nucl{$temp_pid};
    my $ind = $key - $pi_st{$temp_pid};
    substr ($temp_line, $ind, (length $indels{$key}[0]), $indels{$key}[1]);
    $adj_total = $adj_index;
  }
}

open OUTPUT, ">PI.primers.txt";
foreach my $key (sort {$pi_st{$a} <=> $pi_st{$b}} (keys %pi_result)) {
  print OUTPUT "$key\n$pi_result{$key}\n";
}
close OUTPUT;
