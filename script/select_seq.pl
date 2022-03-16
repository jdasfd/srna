#!/usr/bin/perl

use Data::Dumper;
use strict;
use warnings;
use Getopt::Long;

=head1 NAME
select_seq.pl - Filter out the aligned sequence with different tier.

=head 1 SYNOPSIS

perl select_seq.pl -i <.tsv> -t <tier.tsv> -o <.tsv>
Options:
    -i,--input  input tsv file giving chrom and seq
    -t,--tier   tier file contain different sequence
    -o,--out    output files, default is stdout
    -h,--help   help information
=cut

GetOptions(
    "i|input=s" => \(my $input),
    "t|tier=s"  => \(my $tier),
    "o|out=s"   => \(my $output = 'stdout'),
    "h|help"    => \(my $help),
);

sub usage{
    my $help_str = <<"EOF";
    perl select_seq.pl -i <.tsv> -t <tier.tsv> -o <.tsv>
Options:
    -i,--input  input tsv file giving chrom and seq
    -t,--tier   tier file contain different sequence
    -o,--out    output files, default is stdout
    -h,--help   help information
EOF
    return $help_str;
}

die usage() if defined $help;
die ("Please input a tsv") if not defined $input;
die ("Please provide a tsv contained sequence") if not defined $tier;

my (@tier,@for_print);

open my $tier_in, "<", "$tier";
while(<$tier_in>){
    chomp;
    my @array = split/\t/,$_;
    push(@tier, $array[0]);
}
close $tier_in;

open my $input_in, "<", "$input";
while(<$input_in>){
    chomp;
    my @inarray = split/\t/,$_;
    if(grep {$inarray[1] eq $_} @tier){
        my $for_print = "$inarray[0]\t$inarray[1]\n";
        push (@for_print, $for_print)
    }
    else{
        next;
    }
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