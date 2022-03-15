open $t2_in, "<", "/mnt/e/project/srna/output/sequence/tier2.tsv";
while(<$t2_in>){
    chomp;
    push(@tier2, $_);
}
close $t2_in;

while(<>){
    chomp;
    @a = split/\t/,$_;
    if(grep {$a[1] eq $_} @tier2){
        print "$a[0]\t$a[1]\n";
    }
    else{
        next;
    }
}
