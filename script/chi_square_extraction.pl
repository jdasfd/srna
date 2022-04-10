#!/user/bin/perl -w

$/="\=\=\> ";

while(<>){
chomp;
if (!$_){
    next;
    }
else{
    $_ =~ /(SRR.+)\./;
    $file = $1;
    $_ =~ /SRR.+\.(.+)\r?\n/;
    $group = $1;
    $_ =~ /X\-squared\s=\s(.+),\sdf/;
    $chi = $1;
    print"$file\t$group\t$chi\n";
    }
}