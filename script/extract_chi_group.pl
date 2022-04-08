#!/usr/bin/perl

use strict;
use warnings;

my $file = 0;
my (@p1, @p2, @p3, @p4);

while(<>){
    chomp;
    my @inarr = split/\t/, $_;
    if($inarr[0] eq $file){
        if($inarr[2] == 1){
            push (@p1, $inarr[1], $inarr[3], $inarr[4]);
        }
        elsif($inarr[2] == 2){
            push (@p2, $inarr[1], $inarr[3], $inarr[4]);
        }
        elsif($inarr[2] == 3){
            push (@p3, $inarr[1], $inarr[3], $inarr[4]);
        }
        elsif($inarr[2] == 4){
            push (@p4, $inarr[1], $inarr[3], $inarr[4]);
        }
    }
    else{
        print join ("\t", @p1);
        print "\n";
        print join ("\t", @p2);
        print "\n";
        print join ("\t", @p3);
        print "\n";
        print join ("\t", @p4);
        print "\n";
        @p1 = @p2 = @p3 = @p4 = ();
        $file = $inarr[0];
        push(@p1, $inarr[0], "1");
        push(@p2, $inarr[0], "2");
        push(@p3, $inarr[0], "3");
        push(@p4, $inarr[0], "4");
        push(@p1, $inarr[1], $inarr[3], $inarr[4]);
    }
}

END{
    print join ("\t", @p1);
    print "\n";
    print join ("\t", @p2);
    print "\n";
    print join ("\t", @p3);
    print "\n";
    print join ("\t", @p4);
    print "\n";;
}
