#!/usr/bin/perl

use Data::Dumper;
use strict;
use warnings;
use Getopt::Long;

=head1 NAME
extract_genelist.pl - Extract genelist from gff (in tsv format) using bam (in tsv format).

=head 1 SYNOPSIS

perl extract_genelist.pl -g <gff.tsv> -i <bam.tsv> -o <output_file>
Options:
    -b,--bed    bed file giving region
    -f,--file   per-base file contain depth from samtools
    -o,--out    output files, default is stdout
    -h,--help   help information
=cut

GetOptions(
    "g|gene=s"  => \(my $glist),
    "i|input=s" => \(my $tsv),
    "o|out=s"   => \(my $output = 'stdout'),
    "h|help"    => \(my $help),
);

sub usage{
    my $help_str = <<"EOF";
    perl retrive_depth.pl -b <.bed> -t <depth.tsv> -o <output_file>
Options:
    -g,--gene   genelist from a gff (in tsv format)
    -i,--input  input file extract from bam files (in tsv format)
    -o,--out    output files, default is stdout
    -h,--help   help information
EOF
    return $help_str;
}

die usage() if defined $help;
die ("Please provide a genelist file.") if not defined $glist;
die ("Please provide a bam.tsv file.") if not defined $tsv;

my %genelist;

open my $gene_input, "<", "$glist" or die $!;

while(<$gene_input>){
    chomp;
    my @g_array = split/\t/,$_;
    my @g_name = split/;/,$g_array[4];
    if($g_name[0] =~ s/^ID=gene-(.+)/$1/){
        $genelist{$g_name[0]} = [$g_array[0],$g_array[1],$g_array[2],$g_array[3]];
    }
    else{
        next;
    }
}
close $gene_input;

=pod
END{
    print Dumper \%genelist;
}
=cut

open my $bam_input, "<", "$tsv" or die $!;

while(<$tsv>){
    chomp;
    my @tsvarray = split/\t/,$_;
}