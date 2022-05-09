# sRNA reads characteristics

The sRNA sequences has some specific characteristics. We wanted to find if there were patterns of the sRNA distribution.

## Ratio of sRNA reads / bac-reads from different RNA regions of 4 groups

- Using bed of rna to extract mapping reads from different RNA regions

```bash
mkdir -p /mnt/e/project/srna/output/bam/rna
cd /mnt/e/project/srna/output/bam/bacteria

parallel -j 4 " \
samtools view -F 4 -bh -L ../../../annotation/bacteria/trna.bed \
{}.sort.bam > ../rna/{}.trna.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -F 4 -bh -L ../../../annotation/bacteria/rrna.bed \
{}.sort.bam > ../rna/{}.rrna.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -F 4 -bh -L ../../../annotation/bacteria/mrna.bed \
{}.sort.bam > ../rna/{}.mrna.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
cd /mnt/e/project/srna/output/bam/rna

parallel -j 3 " 
samtools sort -@ 4 {1}.bam > {1}.sort.bam 
samtools index {1}.sort.bam
" ::: $(ls *.bam | perl -p -e 's/\.bam$//')

parallel -j 3 " \
rm {}.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

- Use CIGAR to filter alll bac-reads from different region

Three different RNA regions `.bam` files from the previous step.

```bash
mkdir -p /mnt/e/project/srna/output/bam/rna_tsv
cd /mnt/e/project/srna/output/bam/rna

for file in `ls SRR*.bam | perl -p -e 's/\.sort\.bam$//'`
do
samtools view -@ 10 ${file}.sort.bam | \
tsv-filter --not-iregex 6:X > ../rna_tsv/${file}.tsv;
done
```

```bash
cd /mnt/e/project/srna/output/count/
mkdir all rna
cd /mnt/e/project/srna/output/bam/bac_tsv

for file in `ls *.tsv | perl -p -e 's/\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.tsv | \
tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
awk -v name=${name} -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
> ../../count/all/${file}.all.tsv;
done

cd /mnt/e/project/srna/output/bam/rna_tsv

for file in `ls *.trna.tsv | perl -p -e 's/\.trna\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.trna.tsv | \
tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
awk -v name=${name} -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
> ../../count/rna/${file}.trna.tsv;
done

for file in `ls *.rrna.tsv | perl -p -e 's/\.rrna\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.rrna.tsv | \
tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
awk -v name=${name} -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
> ../../count/rna/${file}.rrna.tsv;
done

for file in `ls *.mrna.tsv | perl -p -e 's/\.mrna\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.mrna.tsv | \
tsv-summarize --group-by 3 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 | \
awk -v name=${name} -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
> ../../count/rna/${file}.mrna.tsv;
done
```

- RNA reads ratio in all bac-reads

```bash
cd /mnt/e/project/srna/output/count/rna

for file in `ls *.trna.tsv | perl -p -e 's/\.trna\.tsv$//'`
do
cat ${file}.trna.tsv | tsv-join --filter-file ../all/${file}.all.tsv \
--key-fields 2 --append-fields 3 >> ../bac_reads_group.trna.tsv;
done

for file in `ls *.rrna.tsv | perl -p -e 's/\.rrna\.tsv$//'`
do
cat ${file}.rrna.tsv | tsv-join --filter-file ../all/${file}.all.tsv \
--key-fields 2 --append-fields 3 >> ../bac_reads_group.rrna.tsv;
done

for file in `ls *.mrna.tsv | perl -p -e 's/\.mrna\.tsv$//'`
do
cat ${file}.mrna.tsv | tsv-join --filter-file ../all/${file}.all.tsv \
--key-fields 2 --append-fields 3 >> ../bac_reads_group.mrna.tsv;
done
```

```bash
cd ..

for rna in "trna" "rrna" "mrna"
do
cat bac_reads_group.${rna}.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[4];
print"$a[0]\t$a[1]\t$b\t$a[3]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > bac_ratio_group.${rna}.tsv
done

for rna in "trna" "rrna" "mrna"
do
cat bac_ratio_group.${rna}.tsv | \
tsv-join -H --filter-file plant_30.tsv --key-fields name \
> bac_ratio_group_30.${rna}.tsv
done
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-f bac_ratio_group_30.trna.tsv -t tRNA_region -y "Bac-reads in tRNA" -o ../figure/trna_reads.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f bac_ratio_group_30.rrna.tsv -t rRNA_region -y "Bac-reads in rRNA" -o ../figure/rrna_reads.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f bac_ratio_group_30.mrna.tsv -t mRNA_region -y "Bac-reads in mRNA" -o ../figure/mrna_reads.pdf
```

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

## Sequence among all files

So the main reason I do this step is to select those frequently occurred among all sequence files, that is, the most likely sRNA appeared among *A. tha* sRNA-seq files.

#### All seq files

```bash
mkdir -p /mnt/e/project/srna/output/tier/trna
cd /mnt/e/project/srna/output/bam/rna

parallel -j 4 " \
samtools view -@ 2 {}.trna.sort.bam | \
tsv-select -f 3,10 > ../../tier/trna/{}.trna.tsv \
" ::: $(ls *.trna.sort.bam | perl -p -e 's/\.trna.+bam$//')
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

#### After filtering

```bash
cat result.tier1.tsv | tsv-join -H --filter-file ../count/plant_50.tsv --key-fields 1 > result_50.tier1.tsv
cat result.tier2.tsv | tsv-join -H --filter-file ../count/plant_50.tsv --key-fields 1 > result_50.tier2.tsv
cat result.tier3.tsv | tsv-join -H --filter-file ../count/plant_50.tsv --key-fields 1 > result_50.tier3.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.tier1.tsv -t ">120_files" -y "Bac-reads in tRNA (T1)" -o ../figure/tier1_percent_50.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.tier2.tsv -t "60-120_files" -y "Bac-reads in tRNA (T2)" -o ../figure/tier2_percent_50.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.tier3.tsv -t "<60_files" -y "Bac-reads in tRNA (T3)" -o ../figure/tier3_percent_50.pdf
```

## tRF region

tRF3/5 regions and other_tRNA regions were extracted from the tRNA.bed file according to tRF characteristics.

#### All seq files

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

## Extract tRF sequence from different bacteria

### Get the *A. thaliana* transcripts

```bash
mkdir -p /mnt/e/project/srna/genome/plant_CDS/Atha
cd /mnt/e/project/srna/genome/plant_CDS/Atha
gffread ../../../annotation/plant/Atha/Atha.gff -g ../../plant/Atha/Atha.fna \
-w Atha_transcript.fa -x Atha_CDS.fa
```

### Extract sequence from the tRF3/5 region and high frequency sequence

```bash
mkdir -p /mnt/e/project/srna/output/seq/trf
cd /mnt/e/project/srna/output/seq/trf
mkdir 1 2 3
cd /mnt/e/project/srna/output/bam/trf

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:1 \
> ../../seq/trf/1/{}.tsv \
" ::: $(ls *.trf3_5.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:2 \
> ../../seq/trf/2/{}.tsv \
" ::: $(ls *.trf3_5.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:3 \
> ../../seq/trf/3/{}.tsv \
" ::: $(ls *.trf3_5.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
mkdir -p /mnt/e/project/srna/output/seq/other
cd /mnt/e/project/srna/output/seq/other
mkdir 1 2 3
cd /mnt/e/project/srna/output/bam/trf

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:1 \
> ../../seq/other/1/{}.tsv \
" ::: $(ls *.other_trf.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:2 \
> ../../seq/other/2/{}.tsv \
" ::: $(ls *.other_trf.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:3 \
> ../../seq/other/3/{}.tsv \
" ::: $(ls *.other_trf.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
cd /mnt/e/project/srna/output/seq/trf/1

cat *_1mis.trf3_5.tsv | tsv-select -f 2,4,5 | \
perl -n -e 'chomp;@a=split/\t/,$_;$l=length($a[2]);print"$a[0]\t$a[1]\t$l\t$a[2]\n";' | \
tsv-summarize --group-by 1,2,3,4 --count | \
perl -n -e 'chomp;@a = split/\t/,$_;print">$a[0]_$a[1]_$a[2]\n";print"$a[3]\n"' \
> ../1mis.1.fasta
```

---

## Reads cover depth and position in different RNA region

### Reads depth and position distribution in different RNA

Divide tRNA regions to 10 separate domains. Every domains reads were counted as relative depth of RNA regions.

```bash
cd /mnt/e/project/srna/output/depth/trna

cat *_1mis.trna.txt >> ../1mis.trna.txt
cat *_unali.trna.txt >> ../unali.trna.txt
cat *_aliall.trna.txt >> ../aliall.trna.txt

cd ..

for file in `ls *.txt | perl -p -e 's/\.trna\.txt//'`
do
cat ${file}.trna.txt | tsv-summarize --group-by 1,2 --sum 3 \
> ${file}.trna.tsv
done

parallel -j 3 " \
perl ../../script/depth.pl -b ../../annotation/bacteria/bac_trna.bed \
-t {}.trna.tsv -o {}.trna.depth.tsv \
" ::: $(ls *.trna.tsv | perl -p -e 's/\.trna\.tsv//')

rm *.txt
```

```bash
parallel -j 3 " \
cat {}.depth.tsv | \
tsv-join --filter-file ../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 1,2,4 --sum 3 | tsv-select -f 3,2,4 | \
tsv-join --filter-file ../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 4,2 --sum 3 | sort -nk 1,1 -nk 2,2 | \
sed '1i\group\tpos\tdepth' > {}.sum.tsv \
" ::: $(ls *.depth.tsv | perl -p -e 's/\.dep.+\.tsv$//')
```

```bash
parallel -j 3 " \
Rscript ../../script/line_cover.r -d {}.sum.tsv \
-o {}.pdf -t {} \
" ::: $(ls *.sum.tsv | perl -p -e 's/\.sum\.tsv$//')
```

### Bacteria with tRNA coverage

```bash
cd /mnt/e/project/srna/output/count/trna

for file in `ls *_1mis.tsv | perl -p -e 's/\.tsv//'`
do
cat ${file}.tsv | \
tsv-filter --ne 3:0 | \
tsv-select -f 1,3 \
>> ../1mis_tRNA.name.tsv
done

cd ..

cat 1mis_tRNA.name.tsv | tsv-summarize --group-by 1 --sum 2 | \
sort -r -nk2 | tsv-join --filter-file ../../name.tsv --key-fields 1 --append-fields 2 \
> bac_count.trna.tsv

cat bac_count.trna.tsv | head -n 50 | sed '1i\name\tcount\tgroup'> top50.trna.tsv

cp top50.trna.tsv /mnt/c/Users/59717/Documents/
```

```R
library(ggplot2)
library(readr)
library(dplyr)
library(forcats)
library(ggforce)

t50t <- read_tsv("top50.trna.tsv")

t50t$group <- as.character(t50t$group)

tplot <- t50t %>%
mutate(name = fct_reorder(name, desc(count))) %>%
ggplot(aes(x = name, y = count, fill = group)) +
geom_bar(stat="identity") +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4)) +
geom_col() +
facet_zoom(ylim = c(0, 30000))
# as.character could change the numeric variables to character variables
# mutate in forcats could sort for the bar plot, desc could sort reversely
# facet_zoom could seperate the bar plot into 2 different y axis resolution
# facet_zoom function was from ggforce
```

### tRNA reads statistical information

Reads extracted from tRNA regions were considered as our target sRNAs. First, we used a script to check file

```bash
mkdir -p /mnt/e/project/srna/output/kegg
cd /mnt/e/project/srna/output/target/trna

for file in `ls *.fasta | perl -p -e 's/\.fasta//'`
do
cat ${file}.fasta | perl ../../../script/seq.pl | sort -nk 2 -r > ${file}_seq.tsv;
done

touch all_trna.tsv
cat *_seq.tsv | cut -f 1 >> all_trna.tsv

cat all_trna.tsv | perl -n -e 'while(<>){chomp; my $seq = $_; $hash{$seq}++}
END{
for my $seq (sort keys %hash){print "$seq\t","$hash{$seq}\n";}
}
' | tsv-filter --gt 2:20 | sort -nk 2 -r > ../../kegg/target_trna.tsv


cat *_rrna_seq.tsv | cut -f 1 >> all_rrna.tsv


cat all_rrna.tsv | perl -n -e 'while(<>){chomp; my $seq = $_; $hash{$seq}++}
END{
for my $seq (sort keys %hash){print "$seq\t","$hash{$seq}\n";}
}
' | tsv-filter --gt 2:20 | sort -nk 2 -r > ../kegg/target_rrna.tsv

rm all_trna.tsv all_rrna.tsv
```

```bash
cat target_trna.tsv | cut -f 1 | perl -n -e 'while(<>){chomp; my $seq = $_; $hash{$seq} = length($seq)}
END{
for my $seq (sort keys %hash){print "$seq\t","$hash{$seq}\n";}
}
' | cut -f 2 | perl -n -e 'while(<>){chomp; my $num = $_; $count{$num}++;}
END{
for my $num (sort keys %count){print "$num\t","$count{$num}\n";}
}
' | sed '1i length\tcount'> length.tsv
```

```R
library(ggplot2)
library(readr)

length <- read_tsv("length.tsv")
l <- ggplot(data = length, aes(x = length, y = count))+
geom_point()+
geom_line()+
scale_x_continuous(breaks=seq(23,35,1))

plot(l)
```

```bash
for file in `ls *.fasta | perl -p -e 's/\.fasta$//'`
do
cat ${file}.fasta | perl -n -e '
while(<>){chomp;
if($_ =~ /^>(.+?)$/)
{print"$1\n";}
else{next;}
}' > ../../sequence/${file}.tsv
done
```

```bash
cd /mnt/e/project/srna/output/bam/rna

parallel -j 10 " \
samtools view {}.trna.sort.bam | \
tsv-select -f 1,3 > ../../sequence/{}.bam.tsv \
" ::: $(ls *_1mis.trna.sort.bam | perl -p -e 's/\.trna.+?bam$//')

cd /mnt/e/project/srna/output/sequence

parallel -j 10 " \
tsv-join {}_1mis.tsv --filter-file {}_1mis.bam.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 2 --count | \
tsv-join --filter-file ../../rawname.tsv --key-fields 1 --append-fields 2 \
> {}.tsv \
" ::: $(ls *.tsv | perl -p -e 's/\_.+?tsv$//' | uniq)

rm *_1mis*.tsv
```

---

```R
bar <- ggplot(data = name, mapping = aes(x = 'category', y = count, fill = cate)) +
geom_bar(stat = 'identity', position = 'stack') +
coord_polar(theta = 'y') +
labs(x = '', y = '', title = '') +

```

```bash
cat plot.count | perl -n -e 'chomp;
my $p,$chr;
@a=split(/\t/,$_);
if($a[0]==$chr){
$p++;
print"$a[0]\t$p\t$a[2]\n";
}
else{
$chr=$a[0];
$p=1;
print"$a[0]\t$p\t$a[2]\n";
}
' > chr.count
```

```R
library(readr)
library(ggplot2)

trnad <- read_tsv("1mis.trna.depth")
names(trnad) = c("bac", "pos", "depth")

depth <- ggplot(data = trnad, aes(x = pos, y = depth, group = bac))+
geom_point()+
geom_line()+
theme(legend.position = 'none')

plot(depth)
```

```bash
Rscript -e '
library(readr)
library(ggplot2)
library(dplyr)
library(forcats)

count <- read_tsv("all_seq.count.tsv")
count$num <- as.character(count$num)

tplot <- count %>%
mutate(seq = fct_reorder(seq, desc(count))) %>%
ggplot(aes(x = seq, y = num)) +
geom_bar(stat="identity") +
theme(axis.text.x = NULL)
ggsave(tplot, file = "freq.pdf", width = 7, height = 4)
'
```

## Reads TPM among different groups (waiting for update)

The goal is to filter out those reads among one file.

```bash
cd /mnt/e/project/srna/annotation/bacteria

cat trna.bed | \
perl -e 'my $sum = 0;
while(<>){chomp;@a = split/\t/,$_;
$length = $a[2]-$a[1];
$sum = $sum+$length;
}
print "$sum";
' | awk '{print $1"\ttrna"}'> ../../output/cover/length.tsv

cat rrna.bed | \
perl -e 'my $sum = 0;
while(<>){chomp;@a = split/\t/,$_;
$length = $a[2]-$a[1];
$sum = $sum+$length;
}
print "$sum";
' | awk '{print $1"\trrna"}'>> ../../output/cover/length.tsv

cat mrna.bed | \
perl -e 'my $sum = 0;
while(<>){chomp;@a = split/\t/,$_;
$length = $a[2]-$a[1];
$sum = $sum+$length;
}
print "$sum";
' | awk '{print $1"\tmrna"}'>> ../../output/cover/length.tsv
```

```bash
cd /mnt/e/project/srna/output/cover

perl ../../script/tpm.pl | sed '1ifile\tgroup\tcatgry\tttpm\trtpm\tmtpm' > tpm.tsv

tsv-select -H -f file,group,catgry,ttpm tpm.tsv | sed '1d' | sed '1ifile\tgroup\tcatgry\tnum' > ttpm.tsv
tsv-select -H -f file,group,catgry,rtpm tpm.tsv | sed '1d' | sed '1ifile\tgroup\tcatgry\tnum' > rtpm.tsv
tsv-select -H -f file,group,catgry,mtpm tpm.tsv | sed '1d' | sed '1ifile\tgroup\tcatgry\tnum' > mtpm.tsv

Rscript /mnt/e/project/srna/script/rna_plot.r -f ttpm.tsv -o ../figure/ttpm.pdf -t trna -y TPM
Rscript /mnt/e/project/srna/script/rna_plot.r -f rtpm.tsv -o ../figure/rtpm.pdf -t rrna -y TPM
Rscript /mnt/e/project/srna/script/rna_plot.r -f mtpm.tsv -o ../figure/mtpm.pdf -t mrna -y TPM
```
