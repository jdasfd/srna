#!/usr/bin/perl -w

$index=shift;
open IN,"$index";
while(<IN>){
chomp;
$_=~s/\>//g;
($id,$genus,$spe)=(split/\s/,$_)[0,1,2];
$genus.="\_";
$genus.=$spe;
print"$id\t$genus\n";
}
close IN;
