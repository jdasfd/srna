while(<>){
    chomp;
    my @a = split(/\t/,$_);
    if ($a[2] == 0 && $a[4] == 0){
        next;
    }
    else{
        print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\n";
    }
}