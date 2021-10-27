#/usr/bin/perl

my @a;
my $b;
print "---\n";

while(<>){
    chomp;
    @a=split(/\t/,$_);
    if($a[3]!=0){
        if(defined($b)){
            if($a[0] eq $b){
                print ",","$a[1]-$a[2]";
            }
            else{
                print "\n","$a[0]: ","$a[1]-$a[2]";
            }
        }
        else{
            print "$a[0]: ","$a[1]-$a[2]";
        }
    $b=$a[0];
    }
    else{
        next;
    }
}
print "\n";