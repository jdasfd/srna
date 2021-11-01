#!/user/bin/perl -w

$/="\=\=\> ";

while(<>){
chomp;
if (!$_){
    next;
    }
else{
    $_=~s/^.?\n//g;
    $_=~s/\,//g;
    $_=~s/(NC\_[\d\,\.]+)/$1\#/g;
    $_=~s/(NZ\_[A|C]P[\d\,\.]+)/$1\#/g;
    $_=~s/(all)/$1#/g;
    $_=~s/(X\-squared)/\#$1/g;
    $_=~s/df\s=\s1\s?/\#/g;
    $_=~s/(\n.?\s)$/\#\n/g;
    ($name,$chi,$pval)=(split/\#/,$_)[0,2,3];
    print"$name\t$chi\t$pval\n";
    }
}
