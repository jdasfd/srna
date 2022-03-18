#!/usr/bin/perl

use Data::Dumper;
use strict;
use warnings;
use Getopt::Long;

=head1 NAME
count.pl - Count every ratio of reads.

=head 1 SYNOPSIS

perl count.pl -a <gene.bed> -r <bedcov.tsv> -o <output_file>
Options:
    -a,--all    bed file giving region
    -r,--rest   per-base file contain depth from samtools
    -o,--out    output files, default is stdout
    -h,--help   help information
=cut

GetOptions(
    "a|all=s"   => \(my $all),
    "r|rest=s"  => \(my $rest),
    "o|out=s"   => \(my $output = 'stdout'),
    "h|help"    => \(my $help),
);

sub usage{
    my $help_str = <<"EOF";
    perl retrive_depth.pl -b <.bed> -t <depth.tsv> -o <output_file>
Options:
    -a,--all    bed file giving region
    -r,--rest   per-base file contain depth from samtools
    -o,--out    output files, default is stdout
    -h,--help   help information
EOF
    return $help_str;
}

die usage() if defined $help;
die ("Select all file.") if not defined $all;
die ("Select rest file.") if not defined $rest;

my %num;
my ($sum, $renum);
my @for_print;

open my $all_in, "<", "$all";
while(<$all_in>){
    chomp;
    my @allarray = split/\t/,$_;
    $num{$allarray[0]} = $allarray[1];
    $sum += $allarray[1];
}
close $all_in;

open my $rest_in, "<", "$rest";
while(<$rest_in>){
    chomp;
    $renum = $_;
}
close $rest_in;

$sum = $sum + $renum;
foreach my $keys (keys %num){
    $num{$keys} = $num{$keys}*100/$sum;
    my $for_print = "$keys\t$num{$keys}\n";
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
