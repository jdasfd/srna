#!/usr/bin/perl

use Data::Dumper;
use strict;
use warnings;
use Getopt::Long;

=head1 NAME
bed_coverage.pl - Give out absolute and mean coverage from the bed file.

=head 1 SYNOPSIS

perl bed_coverage.pl -b <gene.bed> -i <bedcov.tsv> -o <output_file>
Options:
    -b,--bed    bed file giving region
    -t,--tsv    tsv file as result of samtools bedov result
    -o,--out    output files as tsv
    -h,--help   help information
=cut

GetOptions(
    "b|bed=s"   => \(my $bed),
    "i|input=s" => \(my $tsv = 'stdin'),
    "o|out=s"   => \(my $output = 'stdout'),
    "h|help"    => \(my $help),
);

sub usage{
    my $help_str = <<"EOF";
    perl retrive_depth.pl -b <.bed> -t <depth.tsv> -o <output_file>
Options:
    -b,--bed    bed file giving region
    -t,--tsv    tsv file contain depth from samtools
    -o,--out    output files, default is stdout
    -h,--help   help information
EOF
    return $help_str;
}

die usage() if defined $help;
die ("Please provide a bed file.") if not defined $bed;

my %length;
my (@chr,@for_print);
my ($trna,$mean);

open my $bed_input, "<", "$bed" or die $!;

while(<$bed_input>){
    chomp;
    my @locus = ();
    my @bedarray = split/\t/,$_;
    @locus = split/;/,$bedarray[9];
    foreach(@locus){
        if($_ =~ /locus_tag=(.+)/){
            $trna = $1;
        }
        else{
            next;
        }
    }
    $#bedarray = 9 || die ("bed should contain 10 cols which extracted from convert2bed");
    if (grep {$bedarray[0] eq $_} @chr){
        my $trna_len = $bedarray[2] - $bedarray[1];
        $length{$bedarray[0]}->{$trna} = $trna_len;
    }
    else{
        my $trna_len = $bedarray[2] - $bedarray[1];
        $length{$bedarray[0]}->{$trna} = $trna_len;
        push (@chr, $bedarray[0]);
    }
}
close $bed_input;

=pod
END{
    print Dumper \%length;
}
=cut

my $input_tsv;
if (lc($tsv) eq "stdin"){
    $input_tsv = *STDIN;
}
else{
    open $input_tsv, "<", "$tsv";
}

while(<$input_tsv>){
    chomp;
    my @tsvarray = split/\t/,$_;
    my @test = ();
    @test = split/;/,$tsvarray[9];
    foreach(@test){
        if($_ =~ /locus_tag=(.+)/){
            $trna = $1;
        }
        else{
            next;
        }
    }
    my $a = $length{$tsvarray[0]}->{$trna};
    my $mean = $tsvarray[10]/$a;
    my $for_print = "$tsvarray[0]\t$trna\t$tsvarray[10]\t$mean\n";
    push (@for_print, $for_print);
}
close $input_tsv;

my $out_tsv;
if (lc($output) eq "stdout"){
    $out_tsv = *STDOUT;
}
else{
    open $out_tsv, ">", $output;
}

for(@for_print){
    chomp;
    my @b = split/\t/,$_;
    if($b[2] == 0 || $b[3] == 0){
        next;
    }
    else{
        print {$out_tsv} "$_\n";
    }
}
close $out_tsv;
