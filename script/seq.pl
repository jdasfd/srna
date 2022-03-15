#!/usr/bin/perl

my %sequence;
my @c;

while(<>){
    chomp;
    @a = split/\t/,$_;
    my $seq = $a[9];
    $sequence{$seq}++;
    }

for my $seq (sort {$sequence{$a} <=> $sequence{$b}} keys %sequence){
    print "$seq\t","$sequence{$seq}\n";
}