- [sRNA alignment using bowtie2](#srna-alignment-using-bowtie2)

	- [Aligning sRNA-seq data to plant genome](#aligning-srna-seq-data-to-plant-genome)

	- [Aligning different reads to bacterial genomes](#aligning-different-reads-to-bacterial-genomes)

	- [Convert sam to bam file for minimum storage stress](#convert-sam-to-bam-file-for-minimum-storage-stress)

- [Plant reads overview](#plant-reads-overview)

	- [Plant reads percentage](#plant-reads-percentage)

# sRNA alignment using bowtie2

Here we use *Arabidopsis thaliana* genome as an example to prepare potential tRFs screening protocol.

**Using Bowtie2 for reads alignment**

Because of the problem of downloading Bowtie. I decided to adopt Bowtie2, though Bowtie has advantages in mapping short reads.

It was relatively fast and convenient in using Bowtie2. Bowtie2 used reads seed region. The only 1 mismatch allowed within seed regions in genome alignment cannot be ignored though, especially analyzing sRNA.

##  Aligning sRNA-seq data to plant genome

First, I aligned sRNA-seq reads to the *A. tha* genome to validate whether there were reads that could not be aligned to the genome.

### Indexing

```bash
cd /mnt/e/project/srna/genome/plant/genome

bowtie-build --quiet --threads 12 Atha.fna Atha
```

* NOTICE

	> I used HPCC for such a number of files. So the bash script were 12 threads (my own computer), and I changed it into 24 threads when submit to  HPCC.
	>
	> Meanwhile, all the rsync step would be optional, depending on your choice.

```bash
mkdir -p /mnt/e/project/srna/output/bam/plant
mkdir -p /mnt/e/project/srna/output/fastq
cd /mnt/e/project/srna/trim
```

* Upload file to the HPCC

```bash
rsync -avP /mnt/e/project/srna/trim wangq@202.119.37.251:jyq/project/srna/bowtie2/
rsync -avP /mnt/e/project/srna/output wangq@202.119.37.251:jyq/project/srna/bowtie2/
rsync -avP /mnt/e/project/srna/genome wangq@202.119.37.251:jyq/project/srna/bowtie2/
rsync -avP /mnt/e/project/srna/annotation wangq@202.119.37.251:jyq/project/srna/bowtie2/
```

### Aligning

* This step could provide us -N 0 (no mismatch in seed region, default: 22) reads aligned.

**plantall.sh**

```bash
parallel -j 3 " \
bowtie2 -q {}_trimmed.fq.gz -N 0 \
-x ../genome/plant/Atha/Atha --al-gz ../output/fastq/{}_plantall.fq.gz \
--threads 4 -S ../output/bam/plant/{}_plantall.sam \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_trimmed.+gz$//')

# --al-gz: alignment fastq.gz extraction
# --no-unal: do not show unalignment reads in the .sam files
```

```bash
bsub -q mpi -n 24 -J pall -o . "bash plantall.sh"
```

Then we need another mapping round for the -N 1 reads (1 mismatch allowed in seed region, default: 22).

**plantall1mis.sh**

```bash
parallel -j 3 " \
bowtie2 -q {}_trimmed.fq.gz -N 1 \
-x ../genome/plant/Atha/Atha --al-gz ../output/fastq/{}_plantall1mis.fq.gz \
--un-gz ../output/fastq/{}_plantunali.fq.gz --threads 4 \
-S ../output/bam/plant/{}_plantall1mis.sam \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_trimmed.+gz$//')
```

```bash
bsub -q mpi -n 24 -J pall1mis -o . "bash plantall1mis.sh"
```

### Extract sequences of only 1 mismatch.

We aligend twice using -N 0 and -N 1 in bowtie2 plant alignment. Because we chose 240 sRNA-seq files, so the default seed region 22 could covered the most reads in our analysis. The reason why I needed to extract the only -N 1 reads was that -N 1 reads could be a pure subset of those reads which were not originally from *A. tha* (got mismatches with plant genome) and could possibly act through miRNA-like mechanism. It was more common that sRNA reads used strigent conmplementary rules to regulate plant gene expression in plants.

```bash
cd /mnt/e/project/srna/output/fastq

bash ../../script/extract_1mis.sh
# the reason why I avoided using parallel was to make sure that memory would not exceed limit.
# Meanwhile, it was not allowed by seqtk for long name during test, so I had to shorten the name.
# Check every extraction from list to fastq avoiding any problems.

: << EOF
parallel -j 2 " \
seqkit grep --quiet -n -v \
-f <(seqkit seq -n {}_plantaliall.fq.gz) \
{}_plantali1mis.fq.gz -o {}_plant1mis.fq.gz \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_.+\.gz$//' | uniq)
# -n, --by-name (default): match by full name instead of just ID
# -j 2 could prevent running out of memory

# aborted! the seqkit was too slow to complete this
EOF

ls SRR*_plant1mis.fq.gz | wc -l
# Better check file numbers you got
# I used to get less files than expectation (maybe run out of memory)
```

```bash
rsync -avP /mnt/e/project/srna/output/fastq/*_plant1mis.fq.gz \
wangq@202.119.37.251:jyq/project/srna/bowtie2/output/fastq
```

##  Aligning different reads to bacterial genomes

After aligned sRNA reads to the plant genome, we split all reads into two different parts, which was reads originally from plants and reads appeared out of nowhere. The goal was to explore where did those reads come from.

We selected bacterial genomes of 161 species from NCBI RefSeq database. From the previous step, we split reads to 3 types: matched seed regions without any mistake, 1 mismatch seed region allowed and unaligned reads.

###  Indexing

```bash
cd /mnt/e/project/srna/genome/bacteria

bowtie2-build --threads 12 --quiet bacteria.fna bacteria
```

###  Aligning

```bash
mkdir -p /mnt/e/project/srna/output/bam/bacteria
cd /mnt/e/project/srna/
```

Aligning unaligned reads to bacteria species.

**unali.sh**

```bash
cd ./output/fastq
parallel -j 4 " \
bowtie2 -q {}_plantunali.fq.gz -N 0 \
-x ../../genome/bacteria/bacteria --threads 6 -S ../bam/bacteria/{}_unali.sam \
" ::: $(ls SRR*_plantunali.fq.gz | perl -p -e 's/_plant.+gz$//')
```

```bash
bsub -q mpi -n 24 -J unali -o . "bash unali.sh"
```

Aligning 1 mismatch allowed reads to bacteria species.

**1mis.sh**

```bash
cd ./output/fastq
parallel -j 4 " \
bowtie2 -q {}_plant1mis.fq.gz -N 0 \
-x ../../genome/bacteria/bacteria --threads 6 -S ../bam/bacteria/{}_1mis.sam \
" ::: $(ls SRR*_plant1mis.fq.gz | perl -p -e 's/_plant.+gz$//')
```

```bash
bsub -q mpi -n 24 -J 1mis -o . "bash 1mis.sh"
```

Aligning perfectly matched reads to bacteria species.

**all.sh**

```bash
cd ./output/fastq
parallel -j 4 " \
bowtie2 -q {}_plantaliall.fq.gz -N 0 \
-x ../../genome/bacteria/bacteria --threads 6 -S ../bam/bacteria/{}_aliall.sam \
" ::: $(ls SRR*_plantaliall.fq.gz | perl -p -e 's/_plant.+gz$//')
```

```bash
bsub -q mpi -n 24 -J all -o . "bash all.sh"
```

## Convert sam to bam file for minimum storage stress.

```bash
cd /mnt/e/project/srna/output/bam/plant

parallel -j 3 " 
samtools sort -@ 4 {1}.sam > {1}.sort.bam 
samtools index {1}.sort.bam 
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')
```

```bash
cd /mnt/e/project/srna/output/bam/bacteria

parallel -j 3 " 
samtools sort -@ 4 {1}.sam > {1}.sort.bam 
samtools index {1}.sort.bam 
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')
```

```bash
rm *.sam
# clear all sam, bam files could be read by samtools view
```

# Plant reads overview

## Plant reads percentage

The sRNA-seq samples we gathered from NCBI were all completed in plant. So the original purpose was to detect the presence of sRNAs in plants. So the alignment first round in [sRNA_mapping_bowtie2.md](https://github.com/jdasfd/srna/blob/main/step_markdown/sRNA_mapping_bowtie2.md) provided two files: -N 0 reads aligning mode and -N 1 reads aligning mode.

Files originated from the -N 1 parameter contained all reads that could align to the plant genome without considering mismatches.


```bash
cd /mnt/e/project/srna/output/count

bash ../../script/plant_reads.sh > plant_reads.csv

cd ../count

Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
ct <- read.csv(args[1])
p <- ggplot(ct, aes(x = name, y = count, 
fill = factor(group, levels = c("unknown","bacteria","plant")))) +
geom_bar(stat = "identity", position = "fill") +
labs(x = "Seq files", y = "reads aligned ratio") +
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p <- p + scale_x_discrete(breaks = NULL) +
scale_fill_manual(name = "reads source",
labels = c("unknown", "bacteria", "plant"),
values = c("gray75", "darkgoldenrod", "seagreen3"))
ggsave(p, file = "../figure/all_file.pdf", width = 9, height = 4)
' all_file.csv
```

## Filter by ratio of aligning to plant

From the above results, it is clear that some sRNA-seq files cannot align to plant in relatively high ratio. A cut-off of more than 50% reads that cannot aligned to plant was used to filter those low quality seq files.

* sRNA-seq reads ratio

```bash
cd /mnt/e/project/srna/output/count

cat all_file.csv | mlr --icsv --otsv cat | \
tsv-summarize -H --group-by name --sum count > all_reads_count.tsv

cat all_file.csv | mlr --icsv --otsv cat | \
tsv-filter -H --str-eq group:plant | \
tsv-join -H --filter-file all_reads_count.tsv \
--key-fields name --append-fields count_sum | \
tsv-select -H -f name,count,count_sum | sed '1d' | \
perl -n -e 'chomp;@a = split/\t/,$_; $ratio=$a[1]/$a[2]*100;
printf("%s\t%.2f\n","$a[0]","$ratio");' | sed '1iname\tratio' > plant_ratio.tsv

cat plant_ratio.tsv | tsv-filter -H --ge ratio:50 > plant_50.tsv
# cut-off 50%

rm all_reads_count.tsv
```

```bash
cat all_file.csv | mlr --icsv --otsv cat | \
tsv-join -H --filter-file plant_50.tsv --key-fields name | \
mlr --itsv --ocsv cat > all_file_50.csv
```

```bash
Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
ct <- read.csv(args[1])
p <- ggplot(ct, aes(x = name, y = count, fill = factor(group, levels = c("unknown","bacteria","plant")))) +
geom_bar(stat = "identity", position = "fill") +
labs(x = "Seq files", y = "reads aligned ratio") +
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p <- p + scale_x_discrete(breaks = NULL) +
scale_fill_manual(name = "reads source",
labels = c("unknown", "bacteria", "plant"),
values = c("gray75", "darkgoldenrod", "seagreen3"))
ggsave(p, file = "../figure/all_file_50.pdf", width = 9, height = 4)
' all_file_50.csv
```

## Plant sRNA reads distribution

*A. tha* annotation is relatively abundant with full information. Using `.gff` file, it is better using gene to calculate col 3 rather than using directly RNA annotation, such as tRNA *et. al.*. It almost the same using two different methods, though there will be a few lines of difference, *e.g.* miRNA will provide you 5p and 3p, but gene will just give you a region. Because of the existence of transcript splicing, extracting annotation of gene could directly give out the different RNA region to reach the goal of deciding region that reads originated. 

* All 240 sRNA-seq files

```bash
cd /mnt/e/project/srna/annotation/plant/Atha

cat Atha.gff | grep -v '#' | tsv-filter --str-eq 3:gene --iregex 9:gene_biotype=tRNA > Atha_trna.gff
cat Atha.gff | grep -v '#' | tsv-filter --str-eq 3:gene --iregex 9:gene_biotype=rRNA > Atha_rrna.gff
cat Atha.gff | grep -v '#' | tsv-filter --str-eq 3:gene --iregex 9:gene_biotype=miRNA > Atha_mirna.gff
cat Atha.gff | grep -v '#' | tsv-filter --str-eq 3:gene --iregex 9:gene_biotype=snRNA > Atha_snrna.gff
cat Atha.gff | grep -v '#' | tsv-filter --str-eq 3:gene --iregex 9:gene_biotype=snoRNA > Atha_snorna.gff
cat Atha.gff | grep -v '#' | tsv-filter --str-eq 3:gene --iregex 9:gene_biotype=lncRNA > Atha_lncrna.gff
cat Atha.gff | grep -v '#' | tsv-filter --str-eq 3:gene --iregex 9:gene_biotype=ncRNA > Atha_ncrna.gff
cat Atha.gff | grep -v '#' | tsv-filter --str-eq 3:gene --not-iregex 9:gene_biotype=snoRNA \
--not-iregex 9:gene_biotype=snRNA \
--not-iregex 9:gene_biotype=miRNA \
--not-iregex 9:gene_biotype=rRNA \
--not-iregex 9:gene_biotype=tRNA \
--not-iregex 9:gene_biotype=lncRNA \
--not-iregex 9:gene_biotype=ncRNA \
> Atha_mrna.gff
```

```bash
parallel -j 4 " \
cat Atha_{}.gff | convert2bed --input=gff --output=bed > Atha_{}.bed \
" ::: $(ls *_*.gff | perl -p -e 's/^Atha_(.+)\.gff/$1/')
```

```bash
mkdir -p /mnt/e/project/srna/output/bam/plantrna
cd /mnt/e/project/srna/output/bam/plant

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/plant/Atha/Atha_trna.bed \
{}.sort.bam > ../plantrna/{}.trna.bam \
" ::: $(ls *_plantall.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/plant/Atha/Atha_rrna.bed \
{}.sort.bam > ../plantrna/{}.rrna.bam \
" ::: $(ls *_plantall.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/plant/Atha/Atha_mrna.bed \
{}.sort.bam > ../plantrna/{}.mrna.bam \
" ::: $(ls *_plantall.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/plant/Atha/Atha_mirna.bed \
{}.sort.bam > ../plantrna/{}.mirna.bam \
" ::: $(ls *_plantall.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/plant/Atha/Atha_snrna.bed \
{}.sort.bam > ../plantrna/{}.snrna.bam \
" ::: $(ls *_plantall.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/plant/Atha/Atha_snorna.bed \
{}.sort.bam > ../plantrna/{}.snorna.bam \
" ::: $(ls *_plantall.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/plant/Atha/Atha_lncrna.bed \
{}.sort.bam > ../plantrna/{}.lncrna.bam \
" ::: $(ls *_plantall.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/plant/Atha/Atha_ncrna.bed \
{}.sort.bam > ../plantrna/{}.ncrna.bam \
" ::: $(ls *_plantall.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
cd /mnt/e/project/srna/output/count
bash ../../script/plant_rna_ratio.sh > plant_rna_ratio.csv

Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
plant <- read.csv(args[1])
p <- ggplot(plant, aes(x = name, y = count, 
fill = factor(group, levels = c("mrna","rrna","ncrna","lncrna","snrna","snorna","trna","mirna")))) +
geom_bar(stat = "identity", position = "fill") +
labs(x = "Seq files", y = "reads aligned to plant ratio") +
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p <- p + scale_x_discrete(breaks = NULL) +
scale_fill_manual(name = "reads source",
labels = c("mrna","rrna","ncrna","lncrna","snrna","snorna","trna","mirna"),
values = c("#70ACAB","#FFF18F","#DAC847","#E3842C","#70795E","#3C3F38","#3A571F","#0B1F25"))
ggsave(p, file = "/mnt/e/project/srna/output/figure/plant_rna_ratio.pdf", width = 9, height = 4)
' plant_rna_ratio.csv
```

* After filtering by cut-off

```bash
cat plant_rna_ratio.csv | mlr --icsv --otsv cat | \
tsv-join -H --filter-file plant_50.tsv --key-fields name | \
mlr --itsv --ocsv cat > plant_rna_ratio_50.csv
```

```bash
Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
plant <- read.csv(args[1])
p <- ggplot(plant, aes(x = name, y = count, 
fill = factor(group, levels = c("mrna","rrna","ncrna","lncrna","snrna","snorna","trna","mirna")))) +
geom_bar(stat = "identity", position = "fill") +
labs(x = "Seq files", y = "reads aligned to plant ratio") +
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p <- p + scale_x_discrete(breaks = NULL) +
scale_fill_manual(name = "reads source",
labels = c("mrna","rrna","ncrna","lncrna","snrna","snorna","trna","mirna"),
values = c("#70ACAB","#FFF18F","#DAC847","#E3842C","#70795E","#3C3F38","#3A571F","#0B1F25"))
ggsave(p, file = "/mnt/e/project/srna/output/figure/plant_rna_ratio_50.pdf", width = 9, height = 4)
' plant_rna_ratio_50.csv
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

## Ratio of reads aligned to bacteria / all non-plant reads

Count all reads numbers from sort.bam files. The goal of this step is to acquire fraction of the reads aligned to bacteria from all reads. A shell script was written to reach the goal. **Attention: bowtie2 results contained -L 22 (default), it means that within 22 seed regions, there were no mismatches.**

* All 240 sRNA-seq files

```bash
mkdir -p /mnt/e/project/srna/output/count
cd /mnt/e/project/srna/output/count
mkdir trna rrna mrna all

cd /mnt/e/project/srna/output/bam/bacteria
bash ../../../script/read_count.sh > ../../count/read_count.csv

# bsub -q mpi -n 24 -o .. -J count "bash read_count.sh | tee ../../count/read_count.csv"
# tee will output results into *.out because of -o in bsub command
```

```bash
cd ../../count

Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
count <- read.csv(args[1])
s <- ggplot (data = count, aes(x = group, y = num)) +
geom_boxplot() +
geom_jitter(aes(color = name)) +
theme(legend.position = "none") +
labs(x = " ", y = "Bacterial reads / all reads")
ggsave(s, file = "../figure/read_count.pdf", width = 7, height = 4)
' read_count.csv
```

* Remove those seq files after filter

```bash
cat read_count.csv | mlr --icsv --otsv cat | \
tsv-join -H --filter-file plant_50.tsv --key-fields name | \
mlr --itsv --ocsv cat > read_count_50.csv
```

```bash
Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
count <- read.csv(args[1])
s <- ggplot (data = count, aes(x = group, y = num)) +
geom_boxplot() +
geom_jitter(aes(color = name)) +
theme(legend.position = "none") +
labs(x = " ", y = "Bacterial reads / all reads")
ggsave(s, file = "../figure/read_count_50.pdf", width = 7, height = 4)
' read_count_50.csv
```

* Give out the basic information among different groups.

```bash
cat read_count_50.csv | mlr --icsv --otsv cat | \
tsv-summarize -H --group-by group --mean 2 --median 2 | mlr --itsv --omd cat
```

| group | num_mean | num_median |
| --- | --- | --- |
| aliall | 6.174503 | 4.23615 |
| 1mis | 3.825783 | 1.5441 |
| unali | 16.5807745 | 9.45485 |