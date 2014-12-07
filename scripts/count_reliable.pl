#!/usr/bin/perl

use strict;
use warnings;

my $i_count;
my $t_count;
my @idts;

open SCR, $ARGV[0];
while (<SCR>) {
  chomp;
  my @block = split /\t/, $_;
  if ($block[0] ne 'PI17') {
    if ($block[7] > 90) {
      $i_count++;
    }
    $t_count++;
    push @idts, $block[7];
  }
}
close SCR;

print "Out of $t_count pairs $i_count pairs are found out to be over 90% reliable.\n";

my $i_mean = 0;
my $i_stdv = &standard_deviation(@idts);

foreach my $num (@idts) {
  $i_mean += $num;
}
$i_mean /= (scalar @idts);

print "The mean of the identity scores is $i_mean and its standard deviation ${i_stdv}.\n";


sub standard_deviation {
my(@numbers) = @_;
#Prevent division by 0 error in case you get junk data
return undef unless(scalar(@numbers));

# Step 1, find the mean of the numbers
my $total1 = 0;
foreach my $num (@numbers) {
$total1 += $num;
}
my $mean1 = $total1 / (scalar @numbers);

# Step 2, find the mean of the squares of the differences
# between each number and the mean
my $total2 = 0;
foreach my $num (@numbers) {
$total2 += ($mean1-$num)**2;
}
my $mean2 = $total2 / (scalar @numbers);

# Step 3, standard deviation is the square root of the
# above mean
my $std_dev = sqrt($mean2);
return $std_dev;
} 
