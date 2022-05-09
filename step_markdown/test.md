# sequence length distribution

cat tier1.tsv | perl -e 'my %len_count;while(<>){chomp;@a=split/\t/,$_;$len=length($a[0]);$len_count{$len}++;}END{for my $len(sort {$len_count{$a} <=> $len_count{$b}} keys %len_count){print"$len\t$len_count{$len}\n";}}' | sed '1ilength\tcount' > seq_length_tier1.tsv

Rscript -e '
library(readr)
library(ggplot2)
args <- commandArgs(T)
count <- read_tsv(args[1])

tplot <- ggplot(count, aes(x = length, y = count)) +
geom_line()
ggsave(tplot, file = "seq_length_tier1.pdf", width = 7, height = 4)
' seq_length_tier1.tsv

# TEST: Extract bam genelist with 1mis hit:

```bash
cat SRR1042171_plant1mis.tsv | \
perl -n -e 'chomp;if($_=~/XM:i:1/){print"$_\n";}else{next;}' | \
tsv-select -f 3,1,2,4,10 > ../SRR1042171_plant1mis.tsv
# 'XM:i:1' in sam/bam files represents 1mis alignment
```



```bash
parallel -j 6 " \
perl /mnt/e/project/srna/script/extract_gene.pl \
-g Atha_genelist.tsv -i SRR1042171_plant1mis.tsv \
-o {}.genelist.tsv \
" ::: $(ls )
```

## Only in tRF region

Remove the intersection of trf3_5 and other_trf.

```bash
mkdir -p /mnt/e/project/srna/output/bam/trf_tsv
cd /mnt/e/project/srna/output/bam/trf

parallel -j 6 " \
samtools view -@ 2 {}.sort.bam > ../trf_tsv/{}.tsv \
" ::: $(ls *.trf3_5.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
parallel -j 6 " \
samtools view -@ 2 {}.sort.bam > ../trf_tsv/{}.tsv \
" ::: $(ls *.other_trf.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
parallel -j 6 " \
cat {}.trf3_5.tsv | tsv-join --filter-file {}.other_trf.tsv --exclude \
> {}.only_trf.tsv \
" ::: $(ls *.trf3_5.tsv | perl -p -e 's/\.tr.+tsv$//')
```

```bash
mkdir -p /mnt/e/project/srna/output/trf/tsv
cd /mnt/e/project/srna/output/bam/trf_tsv

parallel -j 6 " \
cat {}.only_trf.tsv | \
tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-select -f 3,2 | perl ../../../script/filter.pl | \
tsv-join --filter-file ../../count/all/{}.all.tsv --key-fields 1 --append-fields 2 \
> ../../trf/tsv/{}.tsv
" ::: $(ls *.only_trf.tsv | perl -p -e 's/\.on.+tsv$//')
```

```bash
cd /mnt/e/project/srna/output/trf/tsv
bash ../../../script/group_rna_count.sh > ../name_count.only_trf.tsv
```

```bash
cd ..

cat name_count.only_trf.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[3];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.only_trf.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.only_trf.tsv -t "only_trf" -y "Bac-reads in tRNA" -o ../figure/only_trf_percent.pdf
```

Previous outdated method

```R
seq <- read_tsv("all_seq_10.count.tsv")
count <- read_tsv("all_seq_10.num.tsv")
ratio <- read_tsv("all_seq_10.ratio.tsv")

dens <- Mclust(ratio$ratio)

summary(dens)
summary(dens$BIC)
summary(dens, parameters = TRUE)


summary(densi, parameters = TRUE)

------------------------------------------------------- 
Density estimation via Gaussian finite mixture modeling 
------------------------------------------------------- 

Mclust V (univariate, unequal variance) model with 4 components: 

 log-likelihood     n df       BIC       ICL
      -124219.8 34307 11 -248554.5 -260858.1

Mixing probabilities:
        1         2         3         4 
0.1713062 0.3548931 0.3384245 0.1353762 

Means:
       1        2        3        4 
11.42694 14.78338 25.11425 67.11505 

Variances:
           1            2            3            4 
   0.2683464    4.6234189   55.1242285 1341.6427739 

x <- seq(11, 216, by = 1)
y1 <- dnorm(x, 11.42694, 0.2683464)
y2 <- dnorm(x, 14.78338, 4.6234189)
y3 <- dnorm(x, 25.11425, 55.1242285)
y4 <- dnorm(x, 67.11505, 1341.6427739)
data_fun <- data.frame(file = x, values1 = y1, values2 = y2, values3 = y3, values4 = y4, ratio = ratio$ratio, 
                       group = c("g1","g1", rep(c("g2"), each = 13), rep(c("g3"), each = 138), rep(c("g4"), each = 53)))

p <- ggplot(data_fun) +
    geom_bar(aes(x = file, y = ratio, fill = group), stat="identity", color = "black") +
    geom_line(aes(x = file, y = values1), size = 1, color = "hotpink") +
    geom_line(aes(x = file, y = values2), size = 1, color = "orchid3") +
    geom_line(aes(x = file, y = values3), size = 1, color = "royalblue3") +
    geom_line(aes(x = file, y = values4), size = 1, color = "tomato1") +
    coord_cartesian(ylim = c(0,0.02), xlim = c(11,216)) +
    theme(axis.ticks = element_blank()) +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black")) +
    scale_fill_manual(values = c("hotpink", "orchid3", "royalblue3", "tomato1"))



p <- ggplot(num, aes(x = num, y = count)) + 
geom_bar(stat="identity") + 
coord_cartesian(ylim = c(0, 1000)) # coord_cartesian will keep all the data rather than remove them

x <- seq(11, 216, by = 1)
y1 <- dnorm(x, 11.42704, 0.2684187)
y2 <- dnorm(x, 14.78558, 4.6310346)
y3 <- dnorm(x, 25.12289, 55.1790645)
y4 <- dnorm(x, 67.13416, 1341.9265511)
data_fun <- data.frame(file = x, values1 = y1, values2 = y2, values3 = y3, values4 = y4, count = count$count, 
                       group = c("g1","g1", rep(c("g2"), each = 13), rep(c("g3"), each = 138), rep(c("g4"), each = 53)))

ggsave(p, file = "Group_tier.pdf", width = 12, height = 4)

plot(densi, what = "density", data = seq$num, breaks = br)
cdens <- predict(densi, x, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*densi$parameters$pro))
matplot(x, cdens, type = "l", lwd = 1, add = TRUE, lty = 3)

table <- data.frame(file = x, cdens, ratio = ratio$ratio, group = c("g1","g1", rep(c("g2"), each = 6), rep(c("g3"), each = 26), rep(c("g4"), each = 172)))

p <- ggplot(table) +
    geom_bar(aes(x = file, y = ratio, fill = group), stat="identity", color = "black") +
    geom_line(aes(x = file, y = X1), size = 1, color = "hotpink") +
    geom_line(aes(x = file, y = X2), size = 1, color = "orchid3") +
    geom_line(aes(x = file, y = X3), size = 1, color = "royalblue3") +
    geom_line(aes(x = file, y = X4), size = 1, color = "tomato1") +
    scale_x_continuous(breaks = c(seq(11, 221, 10))) +
    coord_cartesian(ylim = c(0,0.03), xlim = c(11,216)) +
    theme(axis.ticks = element_blank()) +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black")) +
    scale_fill_manual(values = c("hotpink", "orchid3", "royalblue3", "tomato1"))

ggsave(p, file = "Group_tier.pdf", width = 12, height = 4)
ggsave(p, file = "Group_tier_all.pdf", width = 12, height = 4)
    
```


```bash
cd /mnt/e/project/srna/output/bam/plant

echo "name,all" > ../../count/plant_allreads.csv
for file in `ls SRR*_plantall.sort.bam | perl -p -e 's/_p.+bam$//'`
do
count=`samtools view -@ 10 --count ${file}_plantall.sort.bam`;
echo "${file},${count}" | tee ../../count/plant_allreads.csv;
done
```

```bash
for num in {mrna,rrna,ncrna,lncrna,snrna,snorna,trna,mirna}
do
cat plant_rna_ratio_50.csv | mlr --icsv --otsv cat | \
tsv-filter -H --str-eq group:${num} | \
tsv-select -f 1,2 | sed '1d' | sed "1iname\t${num}" > plant_${num}.tsv;
cat plant_allreads.csv | mlr --icsv --otsv cat | \
tsv-join -H --filter-file plant_${num}.tsv --key-fields name --append-fields ${num}

# sed "" can pass variable to sed command
```

We first want to summarize all sRNA-seq from different regions.

```bash
cd /mnt/e/project/srna/output/count

cat all_file_plant_rna_ratio_50.csv | mlr --icsv --otsv cat | tsv-summarize -H --group-by name --sum count > all.tmp.tsv
```

```bash
for num in {mrna,rrna,ncrna,lncrna,snrna,snorna,trna,mirna}
do
cat all_file_plant_rna_ratio_50.csv | mlr --icsv --otsv cat | \
tsv-filter -H --str-eq group:${num} | \
tsv-select -f 1,2 | sed '1d' | sed "1iname\t${num}" > plant_${num}.tmp.tsv;
done

cat all.tmp.tsv | tsv-join -H --filter-file plant_mirna.tmp.tsv --key-fields name --append-fields mirna | \
tsv-join -H --filter-file plant_lncrna.tmp.tsv --key-fields name --append-fields lncrna | \
tsv-join -H --filter-file plant_ncrna.tmp.tsv --key-fields name --append-fields ncrna | \
tsv-join -H --filter-file plant_mrna.tmp.tsv --key-fields name --append-fields mrna | \
tsv-join -H --filter-file plant_rrna.tmp.tsv --key-fields name --append-fields rrna | \
tsv-join -H --filter-file plant_trna.tmp.tsv --key-fields name --append-fields trna | \
tsv-join -H --filter-file plant_snorna.tmp.tsv --key-fields name --append-fields snorna | \
tsv-join -H --filter-file plant_snrna.tmp.tsv --key-fields name --append-fields snrna | \
tsv-join -H --filter-file plant_mirna.tmp.tsv --key-fields name --append-fields mirna | 
# sed "" can pass variable to sed command
```

```bash
mkdir -p /mnt/e/project/srna/output/count/bacteria
cd /mnt/e/project/srna/output/bam/bac_tsv

parallel -j 4 " \
cat {}.tsv | tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 \
> ../../count/bacteria/{}.bac.tsv \
" ::: $(ls *.tsv | perl -p -e 's/\.tsv$//')


parallel -j 6 " \
cat {}.bac.tsv | tsv-summarize --group-by 3 --sum 2 | \
sort -nk 1 > {}.group.tsv \
" ::: $(ls *.bac.tsv | perl -p -e 's/\.bac\.tsv$//')
```

```bash
cat SRR10049355_aliall.all.tsv | perl -n -e 'while(<>){@a=split/\t/,$_; $name = $a[0];$catgry = $a[3];
if($a[1] eq all){$all = $a[2];}else{$group{$a[1]} = $a[2];}}
for $key (keys %group){$group{$key} = $group{$key}*100/$all; print "$name\t$key\t$group{$key}\t$catgry";}
'
```

| group | catgry | ratio_mean     |
| ----- | ------ | -------------- |
| 1     | aliall | 3.48643137842  |
| 2     | aliall | 7.28741824153  |
| 3     | aliall | 2.51668786714  |
| 4     | aliall | 0.62532542953  |
| 1     | mis    | 7.17837896075  |
| 2     | mis    | 10.2846318487  |
| 3     | mis    | 10.2209894554  |
| 4     | mis    | 0.329244945806 |
| 1     | unali  | 11.7434182248  |
| 2     | unali  | 9.69985005573  |
| 3     | unali  | 6.82058079933  |
| 4     | unali  | 1.2904385973   |
