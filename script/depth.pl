#!/usr/bin/perl
use Data::Dumper;

open GFF, "<", "/mnt/e/project/srna/annotation/bacteria/bac_trna.bed";
while (<GFF>){
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
=pod
END{
    print Dumper \%data;
}
=cut

open TSV, "<", "/mnt/e/project/srna/output/depth/1mis.trna.tsv";
while(<TSV>){
    chomp;
    my @d = split/\t/,$_;
    my $ref = \%data;
    foreach $bac (keys %data){
        if($d[0] eq $bac){
            $trna = $ref->{$bac};
            foreach $posinfo (keys %$trna){
                $start = $trna->{$posinfo}->[0];
                $end = $trna->{$posinfo}->[1];
                $strand = $trna->{$posinfo}->[2];
                # print "$posinfo,$start,$end,$strand\n";
                if($d[1] >= $start && $d[1] <= $end){
                    if($strand eq "+"){
                        $pos = $d[1] - $start + 1;
                    }
                    else{
                        $pos = $end - $d[1] + 1;
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
    print "$d[0]\t$pos\t$d[2]\n";
}
close TSV;
