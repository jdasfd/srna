- [Aligning sRNA-seq to plant genome](#aligning-srna-seq-to-plant-genome)

	- [Indexing](#indexing)

	- [Aligning](#aligning)

	- [Plant reads percentage](#plant-reads-percentage)

- [Aligning different reads to bacterial genomes](#aligning-different-reads-to-bacterial-genomes)

	- [Extract sequences of only 1 mismatch](#extract-sequences-of-only-1-mismatch)

# Aligning sRNA-seq to plant genome

Here we use *Arabidopsis thaliana* genome as an example to prepare potential tRFs screening protocol.

**Using Bowtie2 for reads alignment**

Because of the problem of downloading Bowtie. I decided to adopt Bowtie2, though Bowtie has advantages in mapping short reads.

It was relatively fast and convenient in using Bowtie2. Bowtie2 used reads seed region. The only 1 mismatch allowed within seed regions in genome alignment cannot be ignored though, especially analyzing sRNA.


First, I aligned sRNA-seq reads to the *A. tha* genome to validate whether there were reads that could not be aligned to the genome.

## Indexing

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

## Aligning

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

## Convert sam to bam file for minimum storage stress

```bash
cd /mnt/e/project/srna/output/bam/plant

parallel -j 3 " 
samtools sort -@ 4 {1}.sam > {1}.sort.bam 
samtools index {1}.sort.bam 
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')
```

## Plant reads percentage

The sRNA-seq samples we gathered from NCBI were all completed in plant. So the original purpose was to detect the presence of sRNAs in plants. There were files orginated from two different parameters: -N 0 reads aligning mode and -N 1 reads aligning mode.

Files originated from the -N 1 parameter contained all reads that could align to the plant genome included all reads (even mismatches existed, details of -N 1 in bowtie2 could be seen in bowtie2 help information).


```bash
mkdir -p /mnt/e/project/srna/output/count
cd /mnt/e/project/srna/output/count
mkdir trna rrna mrna all
cd ..

bash ../../script/plant_reads.sh > plant_reads.csv

cat plant_reads.csv | perl -n -e 'chomp;
@a = split/,/,$_;
if($a[0] eq "name"){print "$_,num\n";}
else{if(not defined $name){$i = 1;
$name = $a[0];print"$_,$i\n";}
else{if($a[0] eq $name){print"$_,$i\n";}
else{$i++;$name=$a[0];print"$_,$i\n";}}}
' > all_file_plant_reads.csv
```

```bash
Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
ct <- read.csv(args[1])
p <- ggplot(ct, aes(x = num, y = count, 
fill = factor(group, levels = c("unknown","plant")))) +
geom_bar(stat = "identity", position = "fill") +
labs(x = "File numbers", y = "Ratio of reads source")
p <- p + scale_x_continuous(breaks = seq(0, 240, 10)) +
scale_fill_manual(name = "reads source",
labels = c("non-plant", "plant"),
values = c("gray75", "seagreen3")) +
theme(panel.background = element_blank(),
axis.text.x  = element_text(angle=45, vjust=0.5))
ggsave(p, file = "../figure/all_file_plant_reads.pdf", width = 9, height = 4)
' all_file_plant_reads.csv
```

## Filter by ratio of aligning to plant

From the above results, it is clear that some sRNA-seq files cannot align to plant in relatively high ratio. A cut-off of more than 50% reads that cannot aligned to plant was used to filter those low quality seq files.

* sRNA-seq reads ratio

```bash
cd /mnt/e/project/srna/output/count

cat all_file_plant_reads.csv | mlr --icsv --otsv cat | \
tsv-summarize -H --group-by name --sum count > all_file_plant_sum.tsv

cat all_file_plant_reads.csv | mlr --icsv --otsv cat | \
tsv-filter -H --str-eq group:plant | \
tsv-join -H --filter-file all_file_plant_sum.tsv \
--key-fields name --append-fields count_sum | \
tsv-select -H -f name,count,count_sum,num | sed '1d' | \
perl -n -e 'chomp;@a = split/\t/,$_; $ratio=$a[1]/$a[2]*100;
printf("%s\t%.2f\t%s\n","$a[0]","$ratio",$a[3]);' | sed '1iname\tratio\tnum' > plant_ratio.tsv

cat plant_ratio.tsv | tsv-filter -H --ge ratio:50 > all_file_plant_ratio_50.tsv
# cut-off 50%
```

```bash
cat all_file_plant_reads.csv | mlr --icsv --otsv cat | \
tsv-join -H --filter-file all_file_plant_ratio_50.tsv --key-fields name | \
mlr --itsv --ocsv cat > all_file_plant_reads_50.csv
```

```bash
Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
ct <- read.csv(args[1])
p <- ggplot(ct, aes(x = num, y = count, 
fill = factor(group, levels = c("unknown","plant")))) +
geom_bar(stat = "identity", position = "fill") +
labs(x = "File numbers", y = "Ratio of reads source")
p <- p + scale_x_continuous(breaks = seq(0, 240, 10)) +
scale_fill_manual(name = "reads source",
labels = c("non-plant", "plant"),
values = c("gray75", "seagreen3")) +
theme(panel.background = element_blank(),
axis.text.x  = element_text(angle=45, vjust=0.5))
ggsave(p, file = "../figure/all_file_plant_reads_50.pdf", width = 9, height = 4)
' all_file_plant_reads_50.csv
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

* Count reads from different plant region

Here we used `SRR*_plantall.sort.bam` because of we wanted to check the different sRNA ratio - whether they were all miRNA or not. 

```bash
cd /mnt/e/project/srna/output/count
bash ../../script/plant_rna_ratio.sh | tee plant_rna_ratio.csv

cat plant_rna_ratio.csv | perl -n -e 'chomp;
@a = split/,/,$_;
if($a[0] eq "name"){print "$_,num\n";}
else{if(not defined $name){$i = 1;
$name = $a[0];print"$_,$i\n";}
else{if($a[0] eq $name){print"$_,$i\n";}
else{$i++;$name=$a[0];print"$_,$i\n";}}}
' > all_file_plant_rna_ratio.csv
```

```bash
Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
plant <- read.csv(args[1])
p <- ggplot(plant, aes(x = num, y = count, 
fill = factor(group, levels = c("mrna","rrna","ncrna","lncrna","snrna","snorna","trna","mirna")))) +
geom_bar(stat = "identity", position = "fill") +
labs(x = "Seq files", y = "Ratio of reads from plant regions")
p <- p + scale_x_continuous(breaks = seq(0, 240, 10)) +
scale_fill_manual(name = "RNA region",
labels = c("mRNA","rRNA","ncRNA","lncRNA","snRNA","snoRNA","tRNA","miRNA"),
values = c("#70ACAB","#FFF18F","#DAC847","#E3842C","#70795E","#3C3F38","#3A571F","#0B1F25")) +
theme(panel.background = element_blank(),
axis.text.x  = element_text(angle=45, vjust=0.5))
ggsave(p, file = "/mnt/e/project/srna/output/figure/plant_rna_ratio.pdf", width = 9, height = 4)
' all_file_plant_rna_ratio.csv
```

* After filtering by cut-off

```bash
cat all_file_plant_rna_ratio.csv | mlr --icsv --otsv cat | \
tsv-join -H --filter-file all_file_plant_ratio_50.tsv --key-fields name | \
mlr --itsv --ocsv cat > all_file_plant_rna_ratio_50.csv
```

```bash
Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
plant <- read.csv(args[1])
p <- ggplot(plant, aes(x = num, y = count, 
fill = factor(group, levels = c("mrna","rrna","ncrna","lncrna","snrna","snorna","trna","mirna")))) +
geom_bar(stat = "identity", position = "fill") +
labs(x = "Seq files", y = "Ratio of reads from plant regions")
p <- p + scale_x_continuous(breaks = seq(0, 240, 10)) +
scale_fill_manual(name = "RNA region",
labels = c("mRNA","rRNA","ncRNA","lncRNA","snRNA","snoRNA","tRNA","miRNA"),
values = c("#70ACAB","#FFF18F","#DAC847","#E3842C","#70795E","#3C3F38","#3A571F","#0B1F25")) +
theme(panel.background = element_blank(),
axis.text.x  = element_text(angle=45, vjust=0.5))
ggsave(p, file = "/mnt/e/project/srna/output/figure/plant_rna_ratio_50.pdf", width = 9, height = 4)
' all_file_plant_rna_ratio_50.csv
```

## Summary of the plant alignment

* Give out the basic information among different groups.

```bash
cat plant_ratio.tsv | tsv-summarize -H --mean ratio --median ratio | mlr --itsv --omd cat
```

| ratio_mean | ratio_median |
| --- | --- |
| 77.6618333333 | 91.605 |


# Aligning different reads to bacterial genomes

After aligned sRNA reads to the plant genome, we split all reads into two different parts, which was reads originally from plants and reads appeared out of nowhere. The goal was to explore where did those reads come from.

We selected bacterial genomes of 161 species from NCBI RefSeq database. From the previous step, we split reads to 3 types: matched seed regions without any mistake, 1 mismatch seed region allowed and unaligned reads.

## Extract sequences of only 1 mismatch

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

##  Indexing

```bash
cd /mnt/e/project/srna/genome/bacteria

bowtie2-build --threads 12 --quiet bacteria.fna bacteria
```

##  Aligning

**Notice: path should be adjusted appropriately by yourself depending on your device**

```bash
mkdir -p /mnt/e/project/srna/output/bam/bacteria
cd /mnt/e/project/srna/output/fastq
```

Aligning unaligned reads to bacteria species.

**unali.sh**

```bash
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
parallel -j 4 " \
bowtie2 -q {}_plantaliall.fq.gz -N 0 \
-x ../../genome/bacteria/bacteria --threads 6 -S ../bam/bacteria/{}_aliall.sam \
" ::: $(ls SRR*_plantaliall.fq.gz | perl -p -e 's/_plant.+gz$//')
```

```bash
bsub -q mpi -n 24 -J all -o . "bash all.sh"
```

## Convert sam to bam file for minimum storage stress

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

## Ratio of reads aligned to bacteria / all non-plant reads

In the previous step, the reads originated from parameter -N 0 and -N 1 were all added up as the reads of plant genome. Our goal was to explore reads not from plant, so we aligned reads directly to bacteria (previous step).

Count all reads from sort.bam files. The goal of this step is to acquire fraction of the reads aligned to bacteria from all reads. A shell script was written to reach the goal. **Attention: bowtie2 results contained -L 22 (default), it means that within 22 seed regions, there were no mismatches.**

* All 240 sRNA-seq files

```bash
cd /mnt/e/project/srna/output/count
bash ../../script/read_count.sh | tee all_file_bac_reads.csv

# bsub -q mpi -n 24 -o .. -J count "bash read_count.sh | tee ../../count/read_count.csv"
# tee will output results into *.out because of -o in bsub command
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
ggsave(s, file = "../figure/all_file_bac_reads.pdf", width = 7, height = 4)
' all_file_bac_reads.csv
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

