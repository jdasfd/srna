# Checking mistakes occured in the previous alignment

Alignment using `bowtie2` before the parameters were fully understood caused serious errors.

I decided to take following steps for troubleshooting.

## Replication previous results

- Aligning to plant

```bash
mkdir -p /mnt/e/project/srna/output/test/plant
mkdir /mnt/e/project/srna/output/test/fastq
cd /mnt/e/project/srna/trim

parallel -j 3 " \
bowtie2 -q {}_trimmed.fq.gz -N 0 --xeq \
-x ../genome/plant/Atha/Atha --threads 4 \
-S ../output/test/plant/{}_plantall.sam \
--al-gz ../output/test/fastq/{}_aliall.fq.gz \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_trimmed.+gz$//')

parallel -j 3 " \
bowtie2 -q {}_trimmed.fq.gz -N 1 --xeq \
-x ../genome/plant/Atha/Atha --threads 4 \
-S ../output/test/plant/{}_plantall1mis.sam \
--al-gz ../output/test/fastq/{}_plantall1mis.fq.gz \
--un-gz ../output/test/fastq/{}_unali.fq.gz \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_trimmed.+gz$//')
```

```bash
cd /mnt/e/project/srna/output/test/plant

parallel -j 3 " 
samtools sort -@ 4 {1}.sam > {1}.sort.bam 
samtools index {1}.sort.bam 
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')

rm *.sam
```

- Extract mis reads

```bash
cd /mnt/e/project/srna/output/test/fastq

parallel -j 3 " \
seqkit grep -j 2 --quiet -n -v \
-f <(seqkit seq -j 2 -n {}_aliall.fq.gz) \
{}_plantall1mis.fq.gz -o {}_1mis.fq.gz \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_.+\.gz$//' | uniq)
```

- Aligning to bacteria

```bash
mkdir /mnt/e/project/srna/output/test/bacteria
cd /mnt/e/project/srna/output/test/fastq

# unali
parallel -j 3 " \
bowtie2 -q {}_unali.fq.gz --xeq \
-x ../../genome/bacteria/bacteria --threads 4 -S ../bacteria/{}_unali.sam \
" ::: $(ls SRR*_unali.fq.gz | perl -p -e 's/_.+gz$//')

# mis
parallel -j 3 " \
bowtie2 -q {}_1mis.fq.gz --xeq \
-x ../../../genome/bacteria/bacteria --threads 4 -S ../bacteria/{}_1mis.sam \
" ::: $(ls SRR*_1mis.fq.gz | perl -p -e 's/_.+gz$//')

# aliall
parallel -j 3 " \
bowtie2 -q {}_aliall.fq.gz --xeq \
-x ../../genome/bacteria/bacteria --threads 4 -S ../bacteria/{}_aliall.sam \
" ::: $(ls SRR*_1mis.fq.gz | perl -p -e 's/_.+gz$//')
```

```bash
cd /mnt/e/project/srna/output/test/bacteria

parallel -j 3 " 
samtools sort -@ 4 {1}.sam > {1}.sort.bam 
samtools index {1}.sort.bam 
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')

rm *.sam
```

- Split bam according to RNA regions

```bash
mkdir /mnt/e/project/srna/output/test/rna
cd /mnt/e/project/srna/output/test/bacteria

parallel -j 4 " \
samtools view -bh -L ../../../annotation/bacteria/trna.bed \
{}.sort.bam > ../rna/{}.trna.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -bh -L ../../../annotation/bacteria/rrna.bed \
{}.sort.bam > ../rna/{}.rrna.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -bh -L ../../../annotation/bacteria/mrna.bed \
{}.sort.bam > ../rna/{}.mrna.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
cd /mnt/e/project/srna/output/test/rna

parallel -j 3 " 
samtools sort -@ 4 {1}.bam > {1}.sort.bam 
samtools index {1}.sort.bam 
" ::: $(ls *.bam | perl -p -e 's/\.bam$//')

rm *rna.bam
```

- Count reads

```bash
mkdir /mnt/e/project/srna/output/test/count
cd /mnt/e/project/srna/output/test/count
mkdir all trna rrna mrna
cd /mnt/e/project/srna/output/test/bacteria

# all bac_reads
parallel -j 6 " \
samtools idxstats {}.sort.bam | \
tsv-select -f 1,3 | grep -v '*' | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 \
> ../count/all/{}.all.tsv \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

cd /mnt/e/project/srna/output/test/rna

# trna
parallel -j 6 " \
samtools idxstats {}.trna.sort.bam | \
tsv-select -f 1,3 | grep -v '*' | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 \
> ../count/trna/{}.trna.tsv \
" ::: $(ls *.trna.sort.bam | perl -p -e 's/\.trna.+bam$//')

# rrna
parallel -j 6 " \
samtools idxstats {}.rrna.sort.bam | \
tsv-select -f 1,3 | grep -v '*' | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 \
> ../count/rrna/{}.rrna.tsv \
" ::: $(ls *.rrna.sort.bam | perl -p -e 's/\.rrna.+bam$//')

# mrna
parallel -j 6 " \
samtools idxstats {}.mrna.sort.bam | \
tsv-select -f 1,3 | grep -v '*' | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 \
> ../count/mrna/{}.mrna.tsv \
" ::: $(ls *.mrna.sort.bam | perl -p -e 's/\.mrna.+bam$//')
```

```bash
cd /mnt/e/project/srna/output/test/count/trna

for file in `ls *.trna.tsv | perl -p -e 's/\.trna\.tsv//'`
do
cat ../all/${file}.all.tsv | \
tsv-join --filter-file ${file}.trna.tsv --key-fields 1 --append-fields 2 \
> ${file}.tsv
done
rm *.trna.tsv

cd ../rrna
for file in `ls *.rrna.tsv | perl -p -e 's/\.rrna\.tsv//'`
do
cat ../all/${file}.all.tsv | \
tsv-join --filter-file ${file}.rrna.tsv --key-fields 1 --append-fields 2 \
> ${file}.tsv
done
rm *.rrna.tsv

cd ../mrna
for file in `ls *.mrna.tsv | perl -p -e 's/\.mrna\.tsv//'`
do
cat ../all/${file}.all.tsv | \
tsv-join --filter-file ${file}.mrna.tsv --key-fields 1 --append-fields 2 \
> ${file}.tsv
done
rm *.mrna.tsv
```

```bash
cd /mnt/e/project/srna/output/test/count

for rna in "trna" "rrna" "mrna"
do
cd ./${rna}
bash ../../count.sh > ../name_count.${rna}.tsv
cd ..
done

for rna in "trna" "rrna" "mrna"
do
cat name_count.${rna}.tsv | sed '1d' | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[3]*100/$a[2];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.${rna}.tsv
done
```

- Plot

```bash
cd /mnt/e/project/srna/output/test/count

for rna in "trna" "rrna" "mrna"
do
Rscript -e '
library(ggplot2)
library(readr)
library(gridExtra)

args <- commandArgs(T)
rna <- read_tsv(args[1], show_col_type = FALSE)

rna$group <- as.character(rna$group)

plot1 <- ggplot (data = rna, aes(x = group, y = ratio, group = group, fill = group)) +
geom_boxplot() + 
geom_jitter(color = "black", alpha = 0.1, show.legend = FALSE) +
facet_wrap(~catgry) +
theme(legend.position = "none") +
labs(x = "bacterial group", y = "RNA/all_RNA percent")

plot2 <- ggplot (data = rna, aes(x = catgry, y = ratio, fill = catgry)) +
geom_boxplot() +
theme(legend.position = "none") +
labs(x = "reads group", y = "RNA/all_RNA percent")

pdf(args[2], width = 12, height = 3)
grid.arrange(plot1, plot2, ncol = 2, layout_matrix = rbind(c(1,1,1,2)))
dev.off()
' result.${rna}.tsv ${rna}.pdf
done
```

After viewing the plot, it was pretty sure that the 1st alignment results were reproduced.

## Troubleshooting

Because of the 1st alignment had serious errors, I had to exclude all possible mistakes from anywhere.

- Alignment results count

```bash
cd /mnt/e/project/srna/output/test

bash check_trim.sh > ../check/check_trim.tmp.tsv
bash check_bam_1.sh > ../check/check_bam.1.tmp.tsv
bash check_fq_1.sh > ../check/check_fq.1.tmp.tsv
```

```bash
cd /mnt/e/project/srna/output/check

cat check_bam_1.tmp.sh | tsv-join -H --filter-file check_trim.tmp.tsv | \
tsv-join -H --filter-file check_fq.1.tmp.tsv \
> check.1.tsv

cat check.1.tsv | tsv-filter -H --ff-eq bam:trim_fq | wc -l
#204

cat check.1.tsv | tsv-filter -H --ff-eq bam:sum | wc -l
#1

# repeat processes above and get the results named check.2.tsv
cat check.2.tsv | tsv-filter -H --ff-eq bam:sum | wc -l
#204
```

After doing this, I found that the 1st alignment results were abnormal.

- Reads extraction

```bash
mkdir /mnt/e/project/srna/output/check/plant
cd /mnt/e/project/srna/output/check/plant
mkdir all mis unali
cd /mnt/e/project/srna/output/test
mkdir name

# check those files with high tRNA ratio
cat count/result.trna.tsv | tsv-filter -H --gt ratio:50 --eq group:1 > group1_trna_50.tsv

cat group1_trna_50.tsv | wc -l
#84

bash ext_all.sh &
bash ext_mis.sh &
bash ext_un.sh &

rm -rf name
```

- Why always tRNA

```bash
mkdir /mnt/e/project/srna/output/test/bac_mis
cd /mnt/e/project/srna/output/test/rna

for file in `cat ../group1_trna_50.tsv | cut -f 1 | sed '1d'`
do
samtools view -@ 5 ${file}_1mis.trna.sort.bam | \
tsv-filter --ne 2:4 | tsv-select -f 1,3,4,6,10 \
> ../bac_mis/${file}.trna.tsv
done

cd ..
bash bac_or_not.sh > mis_over_50.trna.tsv

cd bac_mis
parallel -j 6 " \
cat {}.trna.tsv | tsv-filter --iregex 4:X | \
tsv-select -f 2,3,4,5 | \
tsv-summarize --group-by 1,2,3,4 --count | \
sort -r -nk 5 > ../seq/{}.not-bac.trna.tsv \
" ::: $(ls *.trna.tsv | perl -p -e 's/\.trna\.tsv$//')

cd ..
cat seq/*.tsv >> seq.tsv
cat seq.tsv | tsv-summarize --group-by 1,2,3,4 --sum 5 > seq.num.tsv
cat seq.num.tsv | tsv-join --filter-file ../../rawname.tsv --key-fields 1 --append-fields 2 | tsv-select -f 6,5 | tsv-join --filter-file ../../name.tsv --key-fields 1 --append-fields 2 | tsv-summarize --group-by 3 --sum 2 | mlr --itsv --omd cat
```

| 1   | 30894 |
| --- | ----- |
| 2   | 23149 |
| 3   | 4090  |
| 4   | 1     |

So almost those non-bac reads were accumulated into group 1 and group 2.

Next question is: what were their characteristics?

- Characristics of non-bac reads in tRNA region

```bash
cd /mnt/e/project/srna/output/test
cat seq.num.tsv | tsv-select -f 4,3,5 > ../check/seq.raw.tsv
```

- Remove those top bacteria

```bash

```
