#!/usr/bin/perl

use Data::Dumper;
use strict;
use warnings;
use Getopt::Long;

=head1 NAME
ratio2count.pl - Count different ratios with bins.

=head 1 SYNOPSIS

perl ratio2count.pl -b <bin_num> -r <ratio.tsv> -o <output_file>
Options:
    -b,--bin    number of the bin range
    -r,--ratio  ratio of the plant reads 
    -o,--out    output files, default is stdout
    -h,--help   help information
=cut

GetOptions(
    "b|bin=n"       => \(my $bin),
    "r|ratio=s"     => \(my $ratio),
    "o|out=s"       => \(my $output = 'stdout'),
    "h|help"        => \(my $help),
);

sub usage{
    my $help_str = <<"EOF";
    perl ratio2count.pl -b <bin_num> -r <ratio.tsv> -o <output_file>
Options:
    -b,--bin    number of the bin range
    -r,--ratio  ratio of the plant reads 
    -o,--out    output files, default is stdout
    -h,--help   help information
EOF
    return $help_str;
}

die usage() if defined $help;
die ("Give a bin number") if not defined $bin;
die ("Give a reasonble bin") if $bin > 100 || $bin <= 0;
die ("Input a file with ratio") if not defined $ratio;

my @ratio;
my $i;
my %count;
my @for_print;

open my $ratio_in, "<", "$ratio";
while(<$ratio_in>){
    chomp;
    my @array = split/\t/,$_;
    if($array[0] =~ /name/){
        next;
    }
    else{
        push (@ratio, $array[1]);
    }
}

my $range = 100/$bin;

for($i = 0; $i < $range; $i++){
    my $n = $i*$bin;
    my $m = ($i + 1)*$bin;
    $count{$i} = 0;
    for(@ratio){
        if($_ <= $m && $_ > $n){
            $count{$i}++;
        }
    }
}

foreach my $a (keys %count){
    my $k = ($a + 1)*$bin;
    my $for_print = "$k\t$count{$a}\n";
    push (@for_print, $for_print);
}

my $out_tsv;
if (lc($output) eq "stdout"){
    $out_tsv = *STDOUT;
}
else{
    open $out_tsv, ">", $output;
}

for(@for_print){
    print {$out_tsv} "$_";
}
close $out_tsv;
