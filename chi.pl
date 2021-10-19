#/usr/bin/perl

while (<>){
    chomp;
    @a=split(/\t/,$_);
    if ($a[0]=~/chrom/){
        next;
    }
    if ($a[0]=~/total.+?/){
        next;
    }
    if ($a[0]=~/[0-9]$/){
        print "$a[0]","\t","$a[1]","\t","$a[2]","\t";
    }
    elsif ($a[0]=~/_region/){
        print "$a[1]","\t","$a[2]\n";
    }
}
