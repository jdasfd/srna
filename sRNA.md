#  sRNA Mapping

Here we use *Arabidopsis thaliana* genome as an example to prepare potential tRFs screening protocol.

**Attention: I have mounted WSL2 to drive E, so all the ~ will be replaced by /mnt/e. Restricted by the hardware,  : |**

##  Using Bowtie2 for reads alignment

Because of the problem of downloading Bowtie. I decided to adopt Bowtie2, though Bowtie has advantages in mapping short reads.

It was relatively fast and convenient in using Bowtie2. The only 1 mismatch allowed in genome alignment cannot be ignored though, especially analyzing sRNA.

###  Aligning sRNA-seq data to plant genome

First, I aligned sRNA-seq reads to the *A. tha* genome to validate whether there were reads that could not be aligned to the genome.

#### Indexing

```bash
cd /mnt/e/project/srna/genome/plant/genome

bowtie-build --quiet --threads 12 Atha.fna Atha
```

#### Aligning

* NOTICE

	> I used HPCC for such a number of files. So the bash script were 12 threads (my own computer), and I changed it into 24 threads when submit to  HPCC.
	>
	> Meanwhile, I ignored the files transmission step using rsync. So the directory path were used as my own computer. You could substitute /mnt/e to ~.

```bash
mkdir -p /mnt/e/project/srna/output/bam/plant
```

Align to plant without mismatch.

```bash
cd ~/jyq/project/srna/trim

parallel -j 4 " \
bowtie -q {}_trimmed.fq.gz -v 0 -x ../genome/plant/Atha/Atha \
--al ../output/fastq/{}_plantaliall.fq --no-unal \
--threads 6 -S ../output/bam/plant/{}_plantall.sam \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_trimmed.+gz$//')
# -v: the max mismatch allowed
```

```bash
bsub -q fat_384 -n 80 -J plantaliall -o . "bash plantall.sh"
# plantall.sh use above command and bsub
```

Align to plant with 1 mismatch.

```bash
cd ~/jyq/project/srna/trim

parallel -j 4 " \
bowtie -q {}_trimmed.fq.gz -v 1 -x ../genome/plant/Atha/Atha \
--al ../output/fastq/{}_plantali1mis.fq --no-unal \
--threads 6 -S ../output/bam/plant/{}_plant1mis.sam \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_trimmed.+gz$//')
```

```bash
bsub -q fat_384 -n 80 -J plant1mis -o . "bash plant1mis.sh"
```

Align to plant with 2 mismatch.

```bash
cd ~/jyq/project/srna/trim

parallel -j 4 " \
bowtie -q {}_trimmed.fq.gz -v 2 -x ../genome/plant/Atha/Atha \
--al ../output/fastq/{}_plantali2mis.fq \
--un ../output/fastq/{}_plantunali.fq \
--threads 6 -S ../output/bam/plant/{}_plant2mis.sam \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_trimmed.+gz$//')
```

```bash
bsub -q fat_384 -n 80 -J plant2mis -o . "bash plant2mis.sh"
```

#### Extracting reads

Extract reads with only 1 mismatch to plant from fastq files

```bash
cd /mnt/e/project/srna/output/fastq

parallel -j 10 " \
gzip {} \
" ::: $(ls SRR*.fq)
# fq.gz will minimum the storage stress
```

```bash
cd /mnt/e/project/srna/output/fastq

parallel -j 3 " \
seqkit grep -j 2 --quiet -n -v \
-f <(seqkit seq -j 2 -n {}_plantaliall.fq.gz) \
{}_plantali1mis.fq.gz -o {}_plant1mis.fq.gz \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_.+\.gz$//' | uniq)
# -n, --by-name (default): match by full name instead of just ID
```

Extract reads with 2 mismatches to plant from fastq files

```bash
cd /mnt/e/project/srna/output/fastq

parallel -j 3 " \
seqkit grep -j 2 --quiet -n -v \
-f <(seqkit seq -j 2 -n {}_plantali1mis.fq.gz) \
{}_plantali2mis.fq.gz -o {}_plant2mis.fq.gz \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_.+\.gz$//' | uniq)
# -n, --by-name (default): match by full name instead of just ID

rm *_plantali1mis.fq
rm *_plantali2mis.fq
```

###  Aligning different reads to bacterial genomes

We chose 365 bacterial genomes from 191 species as our target bacteria. From the previous step, we split reads to 3 types: matched without any mistake, 1 mismatch allowed and unaligned reads (without any match on plant genome). 

####  Indexing

```bash
cd /mnt/e/project/srna/genome/bacteria

bowtie2-build --threads 12 --quiet bacteria.fna bacteria
```

####  Aligning

```bash
mkdir -p /mnt/e/project/srna/output/bam/bacteria
cd /mnt/e/project/srna/output/fastq
```

Aligning unaligned reads to bacteria species.

unali.sh:

```bash
parallel -j 4 " \
bowtie -q {}_plantunali.fq.gz -v 0 \
-x ../../genome/bacteria/bacteria --threads 6 -S ../bam/bacteria/{}_unali.sam \
" ::: $(ls SRR*_plantunali.fq.gz | perl -p -e 's/_plant.+gz$//')
```

```bash
bsub -q mpi -n 24 -J ali -o . "bash unali.sh"
```

Aligning 1 mismatch allowed reads to bacteria species.

1mis.sh:

```bash
parallel -j 4 " \
bowtie2 -q {}_plant1mis.fq.gz -v 0 \
-x ../../genome/bacteria/bacteria --threads 6 -S ../bam/bacteria/{}_1mis.sam \
" ::: $(ls SRR*_plant1mis.fq.gz | perl -p -e 's/_plant.+gz$//')
```

```bash
bsub -q mpi -n 24 -J 1mis -o . "bash 1mis.sh"
```

Aligning perfectly matched reads to bacteria species.

2mis.sh:

```bash
parallel -j 4 " \
bowtie2 -q {}_plant2mis.fq.gz -v 0 \
-x ../../genome/bacteria/bacteria --threads 6 -S ../bam/bacteria/{}_2mis.sam \
" ::: $(ls SRR*_plantaliall.fq.gz | perl -p -e 's/_plant.+gz$//')
```

```bash
bsub -q mpi -n 24 -J 2mis -o . "bash 2mis.sh"
```

aliall.sh:

```bash
parallel -j 4 " \
bowtie -q {}_plantaliall.fq.gz -v 0 \
-x ../../genome/bacteria/bacteria --threads 6 -S ../bam/bacteria/{}_aliall.sam \
" ::: $(ls SRR*_plantaliall.fq.gz | perl -p -e 's/_plant.+gz$//')
```



#### Convert sam to bam file for minimum storage stress.

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

#### Basic percentages of reads in plant or bacteria

```bash
cd /mnt/e/project/srna/output/bam
bash ../../script/all_file_count.sh > ../count/all_file.csv

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
ggsave(p, file = "../figure/all_file.pdf", width = 9, height = 4)
' all_file.csv
```

#### Filter by ratio of alignment

```bash
for file in `ls *.sort.bam | perl -p -e s/_.+bam$//`
do
unknown=`samtools view --count -@ 10 -f 4 ${file}_plantall.sort.bam`;
all=`samtools view --count -@ 10 ${file}_plantall.sort.bam`;
ratio=`echo "scale=4;$unknown*100/$all" | bc | awk '{printf "%.4f", $0}'`;
echo "${file},${ratio}";
done | tee ../../count/plant.csv

cat plant.csv | mlr --icsv --otsv cat | tsv-filter --le 2:50 > plant_50.tsv
```

```bash
cat all_file.csv | mlr --icsv --otsv cat | \
tsv-join -H --filter-file plant_50.tsv --key-fields 1 | \
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
