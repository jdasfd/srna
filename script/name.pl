#!/usr/bin/perl -w

$index=shift;
open IN,"$index";
while(<IN>){
chomp;
$_=~s/\>//g;
($id,$genus,$spe,$cate)=(split/\s/,$_)[0,1,2,3];
$genus.="\_";
$genus.=$spe;
print"$id\t$genus\t$cate\n";
}
close IN;
