#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $all_length;
my %rna;
open my $length_in, '<', "/mnt/e/project/srna/output/cover/length.tsv";
while(<$length_in>){
    chomp;
    my @a = split/\t/,$_;
    $rna{$a[1]} = $a[0];
}
close $length_in;

foreach my $len (keys %rna){
    $all_length += $rna{$len};
}

open my $tsv_in, '<', "/mnt/e/project/srna/output/cover/all_count.tsv";
while(<$tsv_in>){
    my @b = split/\t/,$_;
    if($b[0] eq "file"){
        next;
    }
    else{
        my $scalling = ($b[3]*1000000)/$all_length;
        my ($trpk,$rrpk,$mrpk);
        foreach my $len (keys %rna){
            if($len eq "trna"){
                $trpk = $b[2]*1000/$rna{$len};
            }
            elsif($len eq "rrna"){
                $rrpk = $b[5]*1000/$rna{$len};
            }
            else{
                $mrpk = $b[6]*1000/$rna{$len};
            }
        }
        my $ttpm = $trpk*1000/$scalling;
        my $rtpm = $rrpk*1000/$scalling;
        my $mtpm = $mrpk*1000/$scalling;
        print "$b[0]\t$b[1]\t$b[4]\t$ttpm\t$rtpm\t$mtpm\n";
    }
}