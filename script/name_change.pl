while(<>){
    chomp;
    my @a = split (/\t/,$_);
    if ($a[2] == 1){
        $a[2] = 3;
        print "$a[0]\t$a[1]\t$a[2]\n";
    }
    elsif ($a[2] == 2){
        $a[2] = 1;
        print "$a[0]\t$a[1]\t$a[2]\n";
    }
    elsif ($a[2] == 3){
        $a[2] = 2;
        print "$a[0]\t$a[1]\t$a[2]\n";
    }
    else {
        print "$a[0]\t$a[1]\t$a[2]\n";
    }
}