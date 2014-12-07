#!usr/bin/perl

open TEST, "test_output.txt";

my @right_count;
my @wrong_count;

$right_count[24] = 0;
$wrong_count[24] = 0;

while (<TEST>) {
  @line = split /\t/, $_;
  @loc_temp = split /\.\./, $line[3];
  $line[0] =~ m/.+_c_(?<shg>[^_]+)_/;
  $seq = $+{shg};
  $query_order = ($loc_temp[0] < $loc_temp[1]);
  $data_order = ($line[4] < $line[5]);
  if ($query_order eq $data_order) {
    $right_count[$seq]++;
  } else {
    $wrong_count[$seq]++;
  }
}

close TEST;

foreach (1..24) {
  print "The matching in sequence $_ contained $right_count[$_] protein";
  print "sequences in the right order and $wrong_count[$_] protein sequences";
  print " in the wrong order.\n";
}
