#!/usr/bin/perl
use Data::Dumper;
use strict;
use warnings;
use Getopt::Long;
use List::Util;

GetOptions(
    "f|file=s" => \(my $file),
    "o|out=s" => \(my $output = 'stdout'),
    "h|help"  => \(my $help),
);

sub usage{
    my $help_str = <<"EOF";
    perl occurrence.pl -f <file> -o <output_file>
Options:
    -f,--file   file contain all species
    -o,--out    output files, default is stdout
    -h,--help   help information
EOF
    return $help_str;
}

die usage() if defined $help;
die ("Please provide at least one file cotained coverage.") if not defined $file;

open my $file_in, "<", "$file" or die $!;

my %count;
my @input;
my @species;

while(<$file_in>){
    chomp;
    @input = split/\t/,$_;
    if(grep {$input[0] eq $_} @species){
        if($input[1] != 0){
            $count{$input[0]}++;
        }
        else{
            next;
        }
    }
    else{
        push (@species, $input[0]);
    }
}
close $file_in;

my $out_tsv;
if (lc($output) eq "stdout"){
    $out_tsv = *STDOUT;
}
else{
    open $out_tsv, ">", $output;
}

for my $keys (sort {$count{$a} <=> $count{$b}} keys %count){
    print {$out_tsv} "$keys\t","$count{$keys}\n";
}
