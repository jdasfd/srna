my %num;
open $name_in, "<", "/mnt/e/project/srna/rawname.tsv";
while(<$name_in>){
    chomp;
    @a = split/\t/,$_;
    $num{$a[1]} = 0;
    push (@c, $a[1]);
}
close $name_in;

while(<>){
    chomp;
    @b = split/\t/,$_;
    if(grep {$b[0] eq $_} @c){
        $num{$b[0]} = $b[1];
    }
}

END{
    foreach $name (keys %num){
        print "$name\t","$num{$name}\n";
    }
}