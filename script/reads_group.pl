#!/usr/bin/perl

use warnings;
use strict;

my %group;
my ($name, $catgry, $all);

while(<>){
    chomp;
    my @a = split/\t/,$_;
    $name = $a[0];
    $catgry = $a[3];
    if($a[1] eq "all"){
        $all = $a[2];
        }
    else{
        $group{$a[1]} = $a[2];
        }
}

for my $key (keys %group){
    $group{$key} = $group{$key}*100/$all;
    printf "%s\t%s\t%.4f\t%s\n", $name, $key, $group{$key}, $catgry;
}
