#/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use YAML::Syck;

GetOptions(
    "c|chrom=s"     => \(my $chrom),
    "t|trna=s"      => \(my $trna),
    "i|intersect=s" => \(my $inter),
    "o|output=s"    => \(my $output = "result.csv"),
    "h|help"        => \(my $help),
);

sub usage {
    my $usage =<<EOF;

    This Script is used for combine different csv files into an one-line fourfold table

    OPTION
    ======

        -c|--chrom      the csv file contained chromosome results
        -t|--trna       the csv file contained trna length
        -i|--intersect  the csv file contained trna covered length
        -o|--output     the output file
        -h|--help       help info

EOF
}

die usage() if $help;

sub file_format {
    my $name = shift;
    if ( $name =~ m/.csv$/i ){
        return "csv";
    }else{
        return "txt";
    }
}

sub get_chrom {
    open my $f_h, "<", $
}

sub get_trna {

}
my ($data1,$data2,$data3);

my $file_type1 = file_format($chrom);
if ($file_type1 eq "csv"){
    $data1 = YAML::Syck::LoadFile($chrom);
}else{
    last;
}

my $file_type2 = file_format($trna);
if ($file_type2 eq "yml"){
    $data2 = YAML::Syck::LoadFile($trna);
}else{
    $data2 = LoadTxt($trna);
}

my $file_type3 = file_format($inter);
if ($file_type3 eq "yml"){
    $data3 = YAML::Syck::LoadFile($inter);
}else{
    $data3 = LoadTxt($inter);
}

print "data1\n";
