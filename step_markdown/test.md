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

## TEST: Extract bam genelist with 1mis hit

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

- Plant old count

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

- Old version

```bash
parallel -j 6 " \
samtools idxstats {}.trna.sort.bam | \
tsv-select -f 1,3 | grep -v '*' | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 \
> ../../count/trna/{}.trna.tsv \
" ::: $(ls *.trna.sort.bam | perl -p -e 's/\.trna.+bam$//')

parallel -j 6 " \
samtools idxstats {}.rrna.sort.bam | \
tsv-select -f 1,3 | grep -v '*' | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 \
> ../../count/rrna/{}.rrna.tsv \
" ::: $(ls *.rrna.sort.bam | perl -p -e 's/\.rrna.+bam$//')

parallel -j 6 " \
samtools idxstats {}.mrna.sort.bam | \
tsv-select -f 1,3 | grep -v '*' | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 \
> ../../count/mrna/{}.mrna.tsv \
" ::: $(ls *.mrna.sort.bam | perl -p -e 's/\.mrna.+bam$//')
```

All aligned reads to different bacteria of all regions.

```bash
cd /mnt/e/project/srna/output/count/trna

for file in `ls *.trna.tsv | perl -p -e 's/\.trna\.tsv//'`
do
cat ${file}.trna.tsv | \
tsv-join --filter-file ../all/${file}.all.tsv --key-fields 1 --append-fields 2 \
> ${file}.tsv
done

rm *.trna.tsv

bash ../../../script/group_rna_count.sh > ../name_count.trna.tsv
```

```bash
cd /mnt/e/project/srna/output/count/rrna

for file in `ls *.rrna.tsv | perl -p -e 's/\.rrna\.tsv//'`
do
cat ${file}.rrna.tsv | \
tsv-join --filter-file ../all/${file}.all.tsv --key-fields 1 --append-fields 2 \
> ${file}.tsv
done

rm *.rrna.tsv

bash ../../../script/group_rna_count.sh > ../name_count.rrna.tsv
```

```bash
cd /mnt/e/project/srna/output/count/mrna

for file in `ls *.mrna.tsv | perl -p -e 's/\.mrna\.tsv//'`
do
cat ${file}.mrna.tsv | \
tsv-join --filter-file ../all/${file}.all.tsv --key-fields 1 --append-fields 2 \
> ${file}.tsv
done

rm *.mrna.tsv

bash ../../../script/group_rna_count.sh > ../name_count.mrna.tsv
```

```bash
cd /mnt/e/project/srna/output/count

cat name_count.trna.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[3];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.trna.tsv

cat name_count.rrna.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[3];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.rrna.tsv

cat name_count.mrna.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[3];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.mrna.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.trna.tsv -t tRNA_region -y "Bac-reads in tRNA" -o ../figure/trna_reads.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.rrna.tsv -t rRNA_region -y "Bac-reads in rRNA" -o ../figure/rrna_reads.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.mrna.tsv -t mRNA_region -y "Bac-reads in mRNA" -o ../figure/mrna_reads.pdf
```

- After filtering

```bash
cat result.trna.tsv | tsv-join -H --filter-file plant_50.tsv --key-fields 1 > result_50.trna.tsv
cat result.rrna.tsv | tsv-join -H --filter-file plant_50.tsv --key-fields 1 > result_50.rrna.tsv
cat result.mrna.tsv | tsv-join -H --filter-file plant_50.tsv --key-fields 1 > result_50.mrna.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.trna.tsv -t tRNA_region -y "Bac-reads in tRNA" -o ../figure/trna_reads_50.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.rrna.tsv -t rRNA_region -y "Bac-reads in rRNA" -o ../figure/rrna_reads_50.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.mrna.tsv -t mRNA_region -y "Bac-reads in mRNA" -o ../figure/mrna_reads_50.pdf
```

## Bacteria occurred frequencies among different groups (waiting for update)

```bash
mkdir -p /mnt/e/project/srna/output/occurrence
cd /mnt/e/project/srna/output/count
```

```bash
for file in `ls`
do
if [ -d "$file" ]
then
if [[ "$file" = "all" ]]
then
continue
else
cd ${file};
cat *_1mis.tsv >> ../../occurrence/${file}/1mis.tsv;
cat *_unali.tsv >> ../../occurrence/${file}/unali.tsv;
cat *_aliall.tsv >> ../../occurrence/${file}/aliall.tsv;
cd ..;
fi
fi
done
```

```bash
cd /mnt/e/project/srna/output/occurrence/

for file in `ls`
do
if [[ -d "$file" ]]
then
cd $file;
parallel -j 3 " \
perl ../../../script/occurrence.pl -f {}.tsv -o {}.num.tsv \
" ::: $(ls *.tsv | perl -p -e 's/\.tsv//');
cd ..;
fi
done
```

```bash
for dir in `ls`
do
if [[ -d "$dir" ]]
then
cd $dir;
for file in `ls *.num.tsv | perl -p -e 's/\.num\.tsv//'`
do
cat ${file}.num.tsv | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
awk -v file=$file '{print $1"\t"$2"\t"$3"\t"file}'
done | \
sed '1ispecies\tnum\tgroup\tcatgry'> ../${dir}.num.tsv;
cd ..;
fi
done
```

```bash
parallel -j 3 " \
Rscript /mnt/e/project/srna/script/rna_plot.r -f {}.num.tsv \
-o ../figure/{}_freq.pdf -t {} -y frequencies \
" ::: $(ls *.tsv | perl -p -e 's/\.num\.tsv//')
```

- Old tier

```bash
parallel -j 4 " \
samtools view -@ 2 {}.trna.sort.bam | \
tsv-select -f 3,10 > ../../tier/trna/{}.trna.tsv \
" ::: $(ls *.trna.sort.bam | perl -p -e 's/\.trna.+bam$//')
```

```bash
cd ../among

parallel -j 6 " \
cat {}.tsv | tsv-summarize --group-by 1,2 --count \
> {}.all.tsv \
" ::: $(ls reads_*.tsv | perl -p -e 's/\.tsv$//')
```

```bash
mkdir -p /mnt/e/project/srna/output/tier/among
mkdir -p /mnt/e/project/srna/output/tier/file
cd /mnt/e/project/srna/output/tier/trna

parallel -j 6 " \
cat {}_aliall.trna.tsv {}_1mis.trna.tsv {}_unali.trna.tsv | \
tsv-summarize --group-by 2 --count > ../file/{}.trna.tsv \
" ::: $(ls *.trna.tsv | perl -p -e 's/_.+tsv$//' | uniq)

cd ../file

rm ../among/all_seq.tsv # ensure >> will not just output after the old file
cat *.tsv | tsv-select -f 1 >> ../among/all_seq.tsv
cd ../among
cat all_seq.tsv | tsv-summarize --group-by 1 --count > all_seq.count.tsv

cat all_seq.count.tsv | tsv-summarize --group-by 2 --count | sed '1inum\tcount'> seq_num.tsv
```

```bash
Rscript -e '
library(ggplot2)
library(readr)
library(ggforce)
args <- commandArgs(T)
count <- read_tsv(args[1])
s <- ggplot (data = count, aes(x = num, y = count)) +
geom_bar(stat = "identity") +
theme(axis.ticks.x = element_blank()) +
labs(x = "file num", y = "frequency") +
scale_x_continuous(limits = c(0,220)) +
facet_zoom(ylim = c(0, 2000)) +
geom_col()
ggsave(s, file = "../../figure/seq_distri.pdf", width = 10, height = 4)
' seq_num.tsv
```

```bash
cat all_seq.count.tsv | tsv-filter --ge 2:120 > tier1.tsv
cat all_seq.count.tsv | tsv-filter --ge 2:60 --lt 2:120 > tier2.tsv
cat all_seq.count.tsv | tsv-filter --lt 2:60 > tier3.tsv
```

```bash
mkdir tier1 tier2 tier3
cd /mnt/e/project/srna/output/tier/trna

parallel -j 10 " \
perl /mnt/e/project/srna/script/select_seq.pl -i {}.trna.tsv \
-t ../among/tier1.tsv -o ../tier1/{}.tier1.tsv \
" ::: $(ls *.trna.tsv | perl -p -e 's/\.trna\.tsv$//')

parallel -j 10 " \
perl /mnt/e/project/srna/script/select_seq.pl -i {}.trna.tsv \
-t ../among/tier2.tsv -o ../tier2/{}.tier2.tsv \
" ::: $(ls *.trna.tsv | perl -p -e 's/\.trna\.tsv$//')

parallel -j 10 " \
perl /mnt/e/project/srna/script/select_seq.pl -i {}.trna.tsv \
-t ../among/tier3.tsv -o ../tier3/{}.tier3.tsv \
" ::: $(ls *.trna.tsv | perl -p -e 's/\.trna\.tsv$//')
```

```bash
cd /mnt/e/project/srna/output/tier/tier1

parallel -j 6 " \
cat {}.tier1.tsv | \
tsv-summarize --group-by 1 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-select -f 3,2 | perl ../../../script/filter.pl | \
tsv-join --filter-file ../../count/all/{}.all.tsv --key-fields 1 --append-fields 2 \
> {}.tsv \
" ::: $(ls *.tier1.tsv | perl -p -e 's/\.tier.+tsv$//')

rm *.tier1.tsv

bash ../../../script/group_rna_count.sh > ../name_count.tier1.tsv
```

```bash
cd /mnt/e/project/srna/output/tier/tier2

parallel -j 6 " \
cat {}.tier2.tsv | \
tsv-summarize --group-by 1 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-select -f 3,2 | perl ../../../script/filter.pl | \
tsv-join --filter-file ../../count/all/{}.all.tsv --key-fields 1 --append-fields 2 \
> {}.tsv \
" ::: $(ls *.tier2.tsv | perl -p -e 's/\.tier.+tsv$//')

rm *.tier2.tsv

bash ../../../script/group_rna_count.sh > ../name_count.tier2.tsv
```

```bash
cd /mnt/e/project/srna/output/tier/tier3

parallel -j 6 " \
cat {}.tier3.tsv | \
tsv-summarize --group-by 1 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-select -f 3,2 | perl ../../../script/filter.pl | \
tsv-join --filter-file ../../count/all/{}.all.tsv --key-fields 1 --append-fields 2 \
> {}.tsv \
" ::: $(ls *.tier3.tsv | perl -p -e 's/\.tier.+tsv$//')

rm *.tier3.tsv

bash ../../../script/group_rna_count.sh > ../name_count.tier3.tsv
```

```bash
cd ..

cat name_count.tier1.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[3];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.tier1.tsv

cat name_count.tier2.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[3];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.tier2.tsv

cat name_count.tier3.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[3];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.tier3.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.tier1.tsv -t ">120_files" -y "Bac-reads in tRNA (T1)" -o ../figure/tier1_percent.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.tier2.tsv -t "60-120_files" -y "Bac-reads in tRNA (T2)" -o ../figure/tier2_percent.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.tier3.tsv -t "<60_files" -y "Bac-reads in tRNA (T3)" -o ../figure/tier3_percent.pdf
```

```bash
cat SRR10049355_aliall.trna.tsv | perl -n -e 'while(<>){
    chomp;
    m>     chomp;
>     my @a = split/\t/, $_;
>     my $len = length($a[9]);
>     my $start = $a[3];
>     my $end = $start + $len - 1;
    prin>     print "$a[2]:$start-$end|reads=$a[0]\n";
> }' | head -n 20
```

## tRF region

tRF3/5 regions and other_tRNA regions were extracted from the tRNA.bed file according to tRF characteristics.

- Get yml of tRF and other tRNA regions file

```bash
cat bac_trna.gff | spanr gff stdin -o bac_trna.yml
cat bac_trna.yml | spanr span stdin --op trim -n 22 -o bac_other_trna.yml
# trim 22 from each tRNA as other tRNA regions
spanr compare --op diff bac_trna.yml bac_other_trna.yml -o bac_trf.yml
```

```bash
mkdir -p /mnt/e/project/srna/output/seq/trna
cd /mnt/e/project/srna/output/bam/rna_tsv

parallel -j 6 " \
cat {}_aliall.trna.tsv | perl ../../../script/bam2yml.pl \
> ../../seq/trna/{}_aliall.trna.yml \
" ::: $(cat ../../count/plant_30.tsv | tsv-select -f 1 | sed '1d')

parallel -j 6 " \
cat {}_mis.trna.tsv | perl ../../../script/bam2yml.pl \
> ../../seq/trna/{}_mis.trna.yml \
" ::: $(cat ../../count/plant_30.tsv | tsv-select -f 1 | sed '1d')

parallel -j 6 " \
cat {}_unali.trna.tsv | perl ../../../script/bam2yml.pl \
> ../../seq/trna/{}_unali.trna.yml \
" ::: $(cat ../../count/plant_30.tsv | tsv-select -f 1 | sed '1d')
```

- Intersect from tRNA region and yml

Using spanr for getting intersect in `yml` format.

```bash
cd /mnt/e/project/srna/output/seq/trna

parallel -j 6 " \
spanr compare --op intersect ../../../annotation/bacteria/bac_trf.yml \
{}.trna.yml -o {}.trf.yml \
" ::: $(ls *.trna.yml | perl -p -e 's/\.trna\.yml//')

parallel -j 6 " \
spanr compare --op intersect ../../../annotation/bacteria/bac_other_trna.yml \
{}.trna.yml -o {}.other.yml \
" ::: $(ls *.trna.yml | perl -p -e 's/\.trna\.yml//')
```

- Transfer yml to tsv for join

```bash
mkdir -p ../tsv

for file in `ls *.trf.yml | perl -p -e 's/\.yml$//'`
do
cat ${file}.yml | spanr convert stdin | cut -d "-" -f 1 | \
perl -n -e 'chomp;@a = split/:/,$_;print "$a[0]\t$a[1]\n";' \
> ../tsv/${file}.tsv
done

for file in `ls *.other.yml | perl -p -e 's/\.yml$//'`
do
cat ${file}.yml | spanr convert stdin | cut -d "-" -f 1 | \
perl -n -e 'chomp;@a = split/:/,$_;print "$a[0]\t$a[1]\n";' \
> ../tsv/${file}.tsv
done
```

### Region count and plot

- Region ratio count

```bash
mkdir -p /mnt/e/project/srna/output/count/region
cd /mnt/e/project/srna/output/bam/rna_tsv

for file in `cat ../../count/plant_30.tsv | tsv-select -f 1 | sed '1d'`
do
for catgry in "aliall" "mis" "unali"
do
for region in "trf" "other"
do
cat ${file}_${catgry}.trna.tsv | tsv-select -f 3,4,10 | \
tsv-join --filter-file ../../seq/tsv/${file}_${catgry}.${region}.tsv --key-fields 1,2 | \
tsv-summarize --group-by 1 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
awk -v file=${file} -v catgry=$catgry '{print file"\t"$1"\t"$2"\t"catgry}' \
> ../../count/region/${file}_${catgry}.${region}.tsv
done
done
done
```

```bash
cd /mnt/e/project/srna/output/count/region

for file in `ls *.trf.tsv | perl -p -e 's/\.trf\.tsv$//'`
do
for region in "trf" "other"
do
cat ${file}.${region}.tsv | tsv-join --filter-file ../all/${file}.all.tsv \
--key-fields 2 --append-fields 3 >> ../bac_reads_group.${region}.tsv;
done
done
```

```bash
cd ..

for region in "trf" "other"
do
cat bac_reads_group.${region}.tsv | \
tsv-filter --ne 5:0 | \
perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[4];
printf"%s\t%s\t%.3f\t%s\n",$a[0],$a[1],$b,$a[3];
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > bac_ratio_group.${region}.tsv
done
```

- Plot

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-f bac_ratio_group.trf.tsv -t tRF -y "Bac-reads in tRF" -o ../figure/trf_reads.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f bac_ratio_group.other.tsv -t others -y "Bac-reads in other" -o ../figure/other_reads.pdf
```

- Summary

```bash
for region in "trf" "other"
do
cat bac_ratio_group.${region}.tsv | \
tsv-summarize -H --group-by group,catgry --mean ratio --median ratio | \
mlr --itsv --omd cat
done
```

### All seq files

```bash
cd /mnt/e/project/srna/annotation/bacteria
```

```bash
cat trna.bed | perl -e 'while(<>){
    chomp;
    @a = split/\t/,$_;
    $end1 = $a[1] + 25;
    $end2 = $a[2] - 25;
    print"$a[0]\t$a[1]\t$end1\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\n";
    print"$a[0]\t$end2\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\n";
}' > trf3_5.bed
```

```bash
cat trna.bed | perl -e 'while(<>){
    chomp;
    @a = split/\t/,$_;
    $start = $a[1] + 26;
    $end = $a[2] - 26;
    print"$a[0]\t$start\t$end\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\n";
}' > other_trf.bed
```

```bash
mkdir -p /mnt/e/project/srna/output/bam/trf
cd /mnt/e/project/srna/output/bam/rna

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/bacteria/trf3_5.bed \
{}.trna.sort.bam > ../trf/{}.trf3_5.bam \
" ::: $(ls *.trna.sort.bam | perl -p -e 's/\.trna.+bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/bacteria/other_trf.bed \
{}.trna.sort.bam > ../trf/{}.other_trf.bam \
" ::: $(ls *.trna.sort.bam | perl -p -e 's/\.trna.+bam$//')
```

``` bash
cd /mnt/e/project/srna/output/bam/trf

parallel -j 3 " 
samtools sort -@ 4 {1}.bam > {1}.sort.bam 
samtools index {1}.sort.bam
" ::: $(ls *.bam | perl -p -e 's/\.bam$//')

parallel -j 3 " \
rm {}.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
cd /mnt/e/project/srna/output/count
mkdir trf3_5 other_trf
cd /mnt/e/project/srna/output/bam/trf

parallel -j 6 " \
samtools idxstats {}.trf3_5.sort.bam | \
tsv-select -f 1,3 | grep -v '*' | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 \
> ../../count/trf3_5/{}.trf3_5.tsv \
" ::: $(ls *.trf3_5.sort.bam | perl -p -e 's/\.trf.+bam$//')

parallel -j 6 " \
samtools idxstats {}.other_trf.sort.bam | \
tsv-select -f 1,3 | grep -v '*' | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 \
> ../../count/other_trf/{}.other_trf.tsv \
" ::: $(ls *.other_trf.sort.bam | perl -p -e 's/\.other.+bam$//')
```

```bash
cd /mnt/e/project/srna/output/count/trf3_5

for file in `ls *.trf3_5.tsv | perl -p -e 's/\.trf.+tsv//'`
do
cat ${file}.trf3_5.tsv | \
tsv-join --filter-file ../all/${file}.all.tsv --key-fields 1 --append-fields 2 \
> ${file}.tsv
done

rm *.trf3_5.tsv

bash ../../../script/group_rna_count.sh > ../name_count.trf3_5.tsv
```

```bash
cd /mnt/e/project/srna/output/count/other_trf

for file in `ls *.other_trf.tsv | perl -p -e 's/\.other.+tsv//'`
do
cat ${file}.other_trf.tsv | \
tsv-join --filter-file ../all/${file}.all.tsv --key-fields 1 --append-fields 2 \
> ${file}.tsv
done

rm *.other_trf.tsv

bash ../../../script/group_rna_count.sh > ../name_count.other_trf.tsv
```

```bash
cd /mnt/e/project/srna/output/count

cat name_count.trf3_5.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[3];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.trf3_5.tsv

cat name_count.other_trf.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[3];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.other_trf.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.trf3_5.tsv -t tRF-3/5_region -y "tRF3/5 Bac-reads" -o ../figure/trf3_5_percent.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.other_trf.tsv -t other-tRNA_region -y "other tRNA Bac-reads" -o ../figure/other_trf_percent.pdf
```

#### After filtering

```bash
cat result.trf3_5.tsv | tsv-join -H --filter-file plant_50.tsv --key-fields 1 > result_50.trf3_5.tsv
cat result.other_trf.tsv | tsv-join -H --filter-file plant_50.tsv --key-fields 1 > result_50.other_trf.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.trf3_5.tsv -t tRF-3/5_region -y "tRF3/5 Bac-reads" -o ../figure/trf3_5_percent_50.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.other_trf.tsv -t other-tRNA_region -n 25 -y "other tRNA Bac-reads" -o ../figure/other_trf_percent_50.pdf
```

| group | catgry | ratio_mean      | ratio_median |
| ----- | ------ | --------------- | ------------ |
| 1     | aliall | 11.0068905473   | 6.242        |
| 2     | aliall | 5.79283251232   | 3.269        |
| 3     | mis    | 9.26852747253   | 8.0455       |
| 1     | mis    | 5.63471764706   | 4.6015       |
| 2     | mis    | 4.5817712766    | 3.314        |
| 4     | mis    | 1.61575324675   | 0.855        |
| 3     | unali  | 6.95255367232   | 6.564        |
| 1     | unali  | 10.0633611111   | 10.0515      |
| 2     | unali  | 7.05948333333   | 6.811        |
| 4     | unali  | 6.23801960784   | 3.846        |
| 3     | aliall | 4.70980693069   | 2.925        |
| 4     | aliall | 0.083           | 0.0435       |

| group | catgry | ratio_mean      | ratio_median |
| ---   | ---    | ---             | ---          |
| 2     | aliall | 0.884591133005  | 0.632        |
| 3     | aliall | 0.652705        | 0.44         |
| 4     | aliall | 0.0965391304348 | 0.045        |
| 3     | mis    | 3.56351412429   | 2.687        |
| 2     | mis    | 1.91545049505   | 1.546        |
| 1     | mis    | 1.74139181287   | 1.361        |
| 3     | unali  | 3.75737209302   | 3.4665       |
| 1     | unali  | 3.83183815029   | 3.629        |
| 2     | unali  | 3.42653         | 3.183        |
| 4     | unali  | 3.14217391304   | 2.565        |
| 1     | aliall | 0.300144385027  | 0.224        |
| 4     | mis    | 0.891146341463  | 0.61         |
