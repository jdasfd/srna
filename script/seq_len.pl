#!/usr/bin/perl

use warnings;
use strict;

my %length;

while(<>){
    chomp;
    my $seq_len = length($_);
    $length{$seq_len}++;
}

for my $keys (sort {$length{$a} <=> $length{$b}} keys %length){
    print "$keys\t$length{$keys}\n";
}
