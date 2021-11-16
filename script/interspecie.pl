my $ctg1 = 0;
my $ctg2 = 0;
my $ctg3 = 0;
my $ctg4 = 0;
my $ctgtar1 = 0;
my $ctgtar2 = 0;
my $ctgtar3 = 0;
my $ctgtar4 = 0;

while(<>){
    chomp;
    my @a = split(/\t/,$_);
    if($a[6] == 1){
        if($a[4] == 0){
            $ctg1++;
        }
        else{
            $ctg1++;
            $ctgtar1++;
        }
    }
    elsif($a[6] == 2){
        if($a[4] == 0){
            $ctg2++;
        }
        else{
            $ctg2++;
            $ctgtar2++;
        }
    }
    elsif($a[6] == 3){
        if($a[4] == 0){
            $ctg3++;
        }
        else{
            $ctg3++;
            $ctgtar3++;
        }
    }
    else{
        if($a[4] == 0){
            $ctg4++;
        }
        else{
            $ctg4++;
            $ctgtar4++;
        }
    }
}

$mi1 = $ctg1 - $ctgtar1;
$mi2 = $ctg2 - $ctgtar2;
$mi3 = $ctg3 - $ctgtar3;
$mi4 = $ctg4 - $ctgtar4;
print "1vs2\t","$ctgtar1\t$mi1\t$ctgtar2\t$mi2\n";
print "1vs3\t","$ctgtar1\t$mi1\t$ctgtar3\t$mi3\n";
print "1vs4\t","$ctgtar1\t$mi1\t$ctgtar4\t$mi4\n";
print "2vs3\t","$ctgtar2\t$mi2\t$ctgtar3\t$mi3\n";
print "2vs4\t","$ctgtar2\t$mi2\t$ctgtar4\t$mi4\n";
print "3vs4\t","$ctgtar3\t$mi3\t$ctgtar4\t$mi4\n";