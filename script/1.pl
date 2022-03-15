open $t1_in, "<", "/mnt/e/project/srna/output/sequence/tier1.tsv";
while(<$t1_in>){
    chomp;
    push(@tier1, $_);
}
close $t1_in;

while(<>){
    chomp;
    @a = split/\t/,$_;
    if(grep {$a[1] eq $_} @tier1){
        print "$a[0]\t$a[1]\n";
    }
    else{
        next;
    }
}
