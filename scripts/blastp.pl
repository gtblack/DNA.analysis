#!/usr/bin/perl

$input = $ARGV[0];
$db = $ARGV[1];
$out = $ARGV[2];
$cmd = "/BiO/BiOTools/blast-2.2.28+/bin/blastp -db $db -query $input -out $out";
system $cmd;

