#!/usr/bin/perl
use Data::Dumper;
use strict;
use warnings;
use Getopt::Long;
use List::Util;

=head1 NAME
retrive_depth.pl - Retrive relative depth information from bed files

=head 1 SYNOPSIS

perl retrive_depth.pl -b <.bed> -t <depth.tsv> -o <output_file>
Options:
    -b,--bed    bed file giving region
    -t,--tsv    tsv file contain depth from samtools
    -o,--out    output files
    -h,--help   help information
=cut

GetOptions(
    "b|bed=s" => \(my $bed),
    "t|tsv=s" => \(my $tsv),
    "o|out=s" => \(my $output = 'stdout'),
    "h|help"  => \(my $help),
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
die ("Please provide a tsv file.") if not defined $tsv;

open my $gff_in, "<", "$bed" or die $!;

my (@a,@b,@d);
my %data;
my $rec;

while (<$gff_in>){
    chomp;
    my @a = split/\t/,$_;
    $data{$a[0]} = $rec;
    if (grep {$_ eq $a[0]} @b){
        $rec->{$a[3]} = [$a[1],$a[2],$a[5]];
    }
    else{
        $rec = {};
        $rec->{$a[3]} = [$a[1],$a[2],$a[5]];
        push @b, "$a[0]";
    }
}
close $gff_in;

=pod
END{
    print Dumper \%data;
}
=cut

my $rpos;
my @for_print;

open my $tsv_in, "<", "$tsv" or die $!;
while(<$tsv_in>){
    chomp;
    my @d = split/\t/,$_;
    my $ref = \%data;
    foreach my $bac (keys %data){
        if($d[0] eq $bac){
            my $trna = $ref->{$bac};
            foreach my $posinfo (keys %$trna){
                my $start = $trna->{$posinfo}->[0];
                my $end = $trna->{$posinfo}->[1];
                my $strand = $trna->{$posinfo}->[2];
                my $unit = ($end-$start+1)/10;
                # print "$posinfo,$start,$end,$strand\n";
                if($d[1] >= $start && $d[1] <= $end){
                    if($strand eq "+"){
                        my $pos = $d[1] - $start + 1;
                        $rpos = int($pos/$unit);
                    }
                    else{
                        my $pos = $end - $d[1] + 1;
                        $rpos = int($pos/$unit);
                    }
                }
                else{
                    next;
                }
            }
        }
        else{
            next;
        }
    }
    my $a = "$d[0]\t$rpos\t$d[2]\n";
    push (@for_print,$a);
}
close $tsv_in;

my $out_tsv;
if (lc($output) eq "stdout"){
    $out_tsv = *STDOUT;
}
else{
    open $out_tsv, ">", $output;
}

for(@for_print){
    print {$out_tsv} $_;
}
close $out_tsv;
