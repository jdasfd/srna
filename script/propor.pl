#!/usr/bin/perl

use Data::Dumper;
use strict;
use warnings;
use Getopt::Long;

=head1 NAME
bed_coverage.pl - Give out absolute and mean coverage from the bed file.

=head 1 SYNOPSIS

perl bed_coverage.pl -b <gene.bed> -i <tsv> -o <output_file>
Options:
    -b,--bed    bed file giving region
    -i,--in     match info extract from samtools with 5 cols
    -o,--out    output files, default is stdout
    -h,--help   help information
=cut

GetOptions(
    "b|bed=s"   => \(my $bed),
    "i|input=s" => \(my $file),
    "o|out=s"   => \(my $output = 'stdout'),
    "h|help"    => \(my $help),
);

sub usage{
    my $help_str = <<"EOF";
    perl retrive_depth.pl -b <.bed> -t <depth.tsv> -o <output_file>
Options:
    -b,--bed    bed file giving region
    -i,--in     per-base file contain depth from samtools
    -o,--out    output files, default is stdout
    -h,--help   help information
EOF
    return $help_str;
}

die usage() if defined $help;
die ("Please provide a bed file.") if not defined $bed;
#die ("Please provide a tsv file.") if not defined $file;

my $trna,$trf5,$trf3,$trh,$other;
my @chr;
my %pos;

open my $bed_input, "<", "$bed" or die $!;

while(<$bed_input>){
    chomp;
    my @locus = ();
    my @bedarray = split/\t/,$_;
    @locus = split/;/,$bedarray[9];
    foreach(@locus){
        if($_ =~ /^ID=gene-(.+)/){
            $trna = $1;
        }
        else{
            next;
        }
    }
    $#bedarray = 9 || die ("bed should contain 10 cols which originated from convert2bed");
    my $trna1 = $bedarray[1] + 25;
    my $trna2 = $bedarray[2] - 25;
    if (grep {$bedarray[0] eq $_} @chr){
        $pos{$bedarray[0]}->{$trna} = [$bedarray[1],$trna1,$trna2,$bedarray[2],$bedarray[5]];
    }
    else{
        $pos{$bedarray[0]}->{$trna} = [$bedarray[1],$trna1,$trna2,$bedarray[2],$bedarray[5]];
        push (@chr, $bedarray[0]);
    }
}
close $bed_input;

=pod
END{
    print Dumper \%pos;
}
=cut

open my $tsv_input, "<", "$file" or die $!;
while(<$tsv_input>){
    my @tsvarray = split/\t/,$_;
    my $ref = \%pos;
    foreach keys in %pos{
        if($tsvarray[2] eq $pos{$keys}){
            my $value = $ref->{$keys};
            foreach my $comp (keys %$value){
                my $seq_len = length($tsvarray[4]);
                if($seq_len >= 30 ){
                    $trh++;
                    print 
                }
                else{
                    my $end = $tsvarray[3] + $seq_len;
                    if($comp->[4] eq "+"){
                        if($tsvarray[3] >= $comp->[0] && $end <= $comp->[1]){
                            $trf5++;
                        }
                        elsif($tsvarray[3] >= $comp->[2] && $end <= [3]){
                            $trf3++;
                        }
                        else{
                            $other++;
                        }
                    }
                    else($comp->[4] eq "-"){
                        if($tsvarray[3] >= $comp->[0] && $end <= $comp->[1]){
                            $trf3++;
                        }
                        elsif($tsvarray[3] >= $comp->[2] && $end <= [3]){
                            $trf5++;
                        }
                        else{
                            $other++;
                        }
                    }
                }
            }
        }
    }
}

print "$\n";

my $out_tsv;
if (lc($output) eq "stdout"){
    $out_tsv = <STDOUT>;
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