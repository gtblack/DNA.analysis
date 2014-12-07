#!usr/bin/perl

# The fasta file should be divided into pairs which the needleman-wunsch
# algorithm can accept as input.

open ALIGN, "BLASTp_consecutive.txt";

my @shot_tags;
my @ref_tags;

while (<ALIGN>) {
  chomp;
  my $result = /^(?<pi25611>[^\s]+)\s+(?<pi17>\d+)/;
  if ($result) {
    push @shot_tags, $+{pi25611};
    push @ref_tags, $+{pi17};
  }
}
close ALIGN;

my %shot_amino;
my $temp_tag;
my $temp_input;

open SHOTA, "PI_ATCC_25611_DSM_20706.aa.fasta";
while (<SHOTA>) {
  chomp;
  my $line = $_;
  if (/^>(?<stag>[^\s]+) /) {
    $temp_tag = $+{stag};
    $temp_input = $line . "\n";
  } else {
    chop $line;
    $temp_input .= ($line . "\n");
    $shot_amino{$temp_tag} = $temp_input;
  }
}
close SHOTA;

my %ref_amino;
$temp_tag = "";
$temp_input = "";

open REFA, "PI17.faa";
while (<REFA>) {
  chomp;
  my $line = $_;
  if (/>gi\|(?<rtag>\d+)\|/) {
    $ref_amino{$temp_tag} = ($temp_input . "\n");
    $temp_tag = $+{rtag};
    $temp_input = $line . "\n";
  } else {
    $temp_input .= $line;
  }
}
$ref_amino{$temp_tag} = $temp_input;
close REFA;

my %shot_nucl;

open SHOTN, "Prevotella_Intermedia_ATCC25611_DSM20706.na.fasta";
while (<SHOTN>) {
  chomp;
  my $line = $_;
  if (/^>(?<stag>[^\s]+) /) {
    $temp_tag = $+{stag};
    $temp_input = $line . "\n";
  } else {
    $temp_input .= ($line . "\n");
    $shot_nucl{$temp_tag} = $temp_input;
  }
}
close SHOTN;

my %ref_locs;

open LOC, "PI17.ptt";
while (<LOC>) {
  chomp;
  my $line = $_;
  if (/(?<sloc>\d+)\.\.(?<eloc>\d+)\t(?<order>.)\t\d+\t(?<pid>\d+)/) {
    if ($+{order} eq "+") {
      if ($+{pid} < 387132090) {
	$ref_locs{"1_$+{sloc}-$+{eloc}"} = $+{pid};
      } else {
	$ref_locs{"2_$+{sloc}-$+{eloc}"} = $+{pid};
      }
    } else {
      if ($+{pid} < 387132090) {
	$ref_locs{"1_c$+{eloc}-$+{sloc}"} = $+{pid};
      } else {
	$ref_locs{"2_c$+{eloc}-$+{sloc}"} = $+{pid};
      }
    }
  }
}
close LOC;

my %ref_nucl;
$temp_tag = "";
$temp_input = "";

open REFN, "PI17.ffn";
while (<REFN>) {
  chomp;
  my $line = $_;
  if (/>g.+\|NC_01786(?<chrom>.)\.1\|:(?<loc>[^\s]+) /) {
    $ref_nucl{$temp_tag} = ($temp_input . "\n");
    my $chrom = $+{chrom} + 1;
    $temp_tag = $ref_locs{"${chrom}_$+{loc}"};
    $temp_input = ">gi\|${temp_tag}|:$+{loc}\n";
  } else {
    $temp_input .= $line;
  }
}
$ref_nucl{$temp_tag} = $temp_input;
close REFN;

open OUTPUT, ">PI25611.cons.aa.fasta";
foreach my $stag (@shot_tags) {
  print OUTPUT $shot_amino{$stag};
}
close OUTPUT;

open OUTPUT, ">PI17.cons.faa";
foreach my $rtag (@ref_tags) {
  print OUTPUT $ref_amino{$rtag};
}
close OUTPUT;


open OUTPUT, ">PI25611.cons.na.fasta";
foreach my $stag (@shot_tags) {
  print OUTPUT $shot_nucl{$stag};
}
close OUTPUT;

open OUTPUT, ">PI17.cons.fna";
foreach my $rtag (@ref_tags) {
  print OUTPUT $ref_nucl{$rtag};
}
close OUTPUT;
