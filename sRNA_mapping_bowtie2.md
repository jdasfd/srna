# sRNA alignment using bowtie2

* NOTICE

	> I used HPCC for such a number of files. So the bash script were 12 threads (my own computer), and I changed it into 24 threads when submit to  HPCC.
	>
	> Meanwhile, I ignored the files transmission step using rsync. So the directory path were used as my own computer. You could substitute /mnt/e to ~.

```bash
mkdir -p /mnt/e/project/srna/output/bam/plant
cd /mnt/e/project/srna/trim
```

*alignall.sh:*

```bash
parallel -j 3 " \
bowtie2 -q {}_trimmed.fq.gz -N 0 \
-x ../genome/plant/Atha/Atha --al-gz ../output/fastq/{}_plantaliall.fq.gz \
--no-unal --threads 4 -S ../output/bam/plant/{}_plantall.sam \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_trimmed.+gz$//')

# --al-gz: alignment fastq.gz extraction
# --no-unal: do not show unalignment reads in the .sam files
```

```bash
bsub -q mpi -n 24 -J aliall -o . "bash alignall.sh"
```

Then we need another mapping round for the 1 mismatch allowed.

*align1mis.sh:*

```bash
parallel -j 3 " \
bowtie2 -q {}_trimmed.fq.gz -N 1 \
-x ../genome/plant/Atha/Atha --al-gz ../output/fastq/{}_plantali1mis.fq.gz \
--un-gz ../output/fastq/{}_plantunali.fq.gz --no-unal --threads 4 \
-S ../output/bam/plant/{}_plant1mis.sam \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_trimmed.+gz$//')
```

```bash
bsub -q mpi -n 24 -J ali1mis -o . "bash align1mis.sh"
```

#### Extract sequences of only 1 mismatch.

```bash
cd /mnt/e/project/srna/output/fastq

parallel -j 3 " \
seqkit grep -j 2 --quiet -n -v \
-f <(seqkit seq -j 2 -n {}_plantaliall.fq.gz) \
{}_plantali1mis.fq.gz -o {}_plant1mis.fq.gz \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_.+\.gz$//' | uniq)
# -n, --by-name (default): match by full name instead of just ID

rm *_plantali1mis.fq.gz
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

*unali.sh:*

```bash
parallel -j 3 " \
bowtie2 -q {}_plantunali.fq.gz \
-x ../../genome/bacteria/bacteria --threads 4 -S ../bam/bacteria/{}_unali.sam \
" ::: $(ls SRR*_plantunali.fq.gz | perl -p -e 's/_plant.+gz$//')
```

```bash
bsub -q mpi -n 24 -J ali -o .. "bash unali.sh"
```

Aligning 1 mismatch allowed reads to bacteria species.

*1mis.sh:*

```bash
parallel -j 3 " \
bowtie2 -q {}_plant1mis.fq.gz \
-x ../../genome/bacteria/bacteria --threads 4 -S ../bam/bacteria/{}_1mis.sam \
" ::: $(ls SRR*_plant1mis.fq.gz | perl -p -e 's/_plant.+gz$//')
```

```bash
bsub -q mpi -n 24 -J 1mis -o .. "bash 1mis.sh"
```

Aligning perfectly matched reads to bacteria species.

*all.sh:*

```bash
parallel -j 3 " \
bowtie2 -q {}_plantaliall.fq.gz \
-x ../../genome/bacteria/bacteria --threads 4 -S ../bam/bacteria/{}_aliall.sam \
" ::: $(ls SRR*_plantaliall.fq.gz | perl -p -e 's/_plant.+gz$//')
```

```bash
bsub -q mpi -n 24 -J unali -o .. "bash all.sh"
```

Convert sam to bam file for minimum storage stress.

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