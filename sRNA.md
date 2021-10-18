#  sRNA Mapping

Here we use *Arabidopsis thaliana* genome as an example to prepare potential tRFs screening protocol.

##  Prepare

###  *A. thaliana* genome

For this *A. thaliana* genome, we downloaded from ensembl. I found that some links may not work well. So I copied the ftp site below.

http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/

```bash
cd /mnt/e/project/srna
mkdir -p genome/plant/Atha

cd /mnt/e/project/srna/genome/plant/Atha

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz

gzip -d *
mv GCF_000001735.4_TAIR10.1_genomic.fna Atha.fna
# rename it for better using
```

###  *A. thaliana* gff

Download gff3 format annotation.

```bash
cd /mnt/e/project/srna/
mkdir -p annotation/plant/Atha
cd annotation/plant/Atha

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz

gzip -d GCF_000001735.4_TAIR10.1_genomic.gff.gz
mv GCF_000001735.4_TAIR10.1_genomic.gff Atha.gff
```

###  Biotools in protocol

using SRAtoolkit for to download SRA files.

```bash
brew install sratoolkit
```

But there were problems when using brew to install it.

Another way could be adopted.

```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.0/sratoolkit.2.11.0-ubuntu64.tar.gz
tar -xzvf sratoolkit.2.11.0-ubuntu64.tar.gz
mv sratoolkit.2.11.0-ubuntu64 sratoolkit.2.11.0
cd sratookit.2.11.0/bin

echo "export PATH="$(pwd):$PATH"" >> ~/.bashrc
# reopen a new linux shell
prefetch --help
```



##  Getting sRNA-seq data

###  Converting .sra to .fastq

Using fastq-dump in SRAtoolkit to convert sra files to fastq files.

```bash
cd /mnt/e/project/srna/sRNA/srr
parallel -j 12 "fastq-dump {}" ::: $(ls *.sra)
```

###  Quality control

Use trim galore for automatic adapters cut. Trim galore is a perl wrapper for Fastqc and cutadapt. It could be automatically or manually remove adapters.

```bash
cd /mnt/e/project/srna/srr
parallel -j 4 "
trim_galore --phred33 --small_rna --output_dir ../trim {}
" ::; $(ls *.fastq)
trim_galore --phred33 --small_rna --output_dir ../trim *.fastq
# --small_rna: automatically remove adapters for small rna sequencing
```



##  Using Bowtie

Bowtie is an ultrafast, memory-efficient short read aligner geared toward quickly aligning large sets of short DNA sequences (reads) to large genomes. It aligns 35-base-pair reads to the human genome at a rate of 25 million reads per hour on a typical workstation.

**Why use Bowtie rather than Bowtie2:** For relatively short reads (*e.g.* less than 50 bp), Bowtie 1 is sometimes faster and/or more sensitive. Bowtie 1 does not yet report gapped alignments; this is future work (Bowtie 2). For sRNA, usually we do not consider mismatches.

Because of Bowtie had been outdated and had problem in downloading it. I decided to adopt Bowtie2.

## Aligning sRNA-seq data to plant genome

First, I aligned sRNA-seq reads to the *A. tha* genome to validate whether there were reads that could not be align to the genome.

###  Indexing

```bash
cd /mnt/e/project/srna/genome/plant/Atha/genome

bowtie2-build --threads 12 --quiet Atha.fna Atha
```

###  Aligning

```bash
cd /mnt/e/project/srna/trim

nohup parallel -j 3 "
    bowtie -v 0 -a -m1 --best --strata --threads 4 --quiet \
    -x ../../plant/Atha/genome/Atha {} -S ../output/bowtie/{}.sam
" ::: $(ls SRR1004935*.fq)
# -v or -n: decide allowed max mismatches. Usually in sRNA it is 0
# -a: output all possible alignment
# -k: if there are multi-alignments, output all good aligments up to k
# -m: if there are multi-alignments, output all alignments more than m
# -S: report .sam file
# --best --strata: best means alignments sort from high score to low score (ties broken by quality). strata must use with best, means report the best part only
# bowtie can accept different file format input, included fastq and fasta

nohup parallel -j 3 "
    bowtie -v 0 -a -m1 --best --strata --threads 4 --quiet \
    -x ../../plant/Atha/genome/Atha {} --al ../output/bowtie/{}.ali.fastq \
" ::: $(ls SRR1004935*.fq)
# output aligned reads to the fastq file

nohup parallel -j 3 "
    bowtie -v 0 -a -m1 --best --strata --threads 4 --quiet \
    -x ../../plant/Atha/genome/Atha {} --un ../output/bowtie/{}.unali.fastq \
" ::: $(ls SRR1004935*.fq)
# output unaligned reads to the fastq file
```

Using Bowtie2

```bash
cd /mnt/e/project/srna/trim

nohup parallel -j 3 "
    bowtie2 -q {} --end-to-end -x ../genome/plant/Atha/Atha \
    --un ../output/{}.unali.fq --threads 4 -S ../output/{}.sam
" ::: $(ls SRR1004935*.fq | perl -p -e 's/_trimmed\.fq$//')

23131748 reads; of these:
  23131748 (100.00%) were unpaired; of these:
    1164936 (5.04%) aligned 0 times
    10693063 (46.23%) aligned exactly 1 time
    11273749 (48.74%) aligned >1 times
94.96% overall alignment rate
21565629 reads; of these:
  21565629 (100.00%) were unpaired; of these:
    1059845 (4.91%) aligned 0 times
    8465079 (39.25%) aligned exactly 1 time
    12040705 (55.83%) aligned >1 times
95.09% overall alignment rate
23119819 reads; of these:
  23119819 (100.00%) were unpaired; of these:
    956732 (4.14%) aligned 0 times
    11659236 (50.43%) aligned exactly 1 time
    10503851 (45.43%) aligned >1 times
95.86% overall alignment rate
22719863 reads; of these:
  22719863 (100.00%) were unpaired; of these:
    1307511 (5.75%) aligned 0 times
    9259237 (40.75%) aligned exactly 1 time
    12153115 (53.49%) aligned >1 times
94.25% overall alignment rate
22224601 reads; of these:
  22224601 (100.00%) were unpaired; of these:
    716329 (3.22%) aligned 0 times
    9437770 (42.47%) aligned exactly 1 time
    12070502 (54.31%) aligned >1 times
96.78% overall alignment rate
```

##  Aligning unaligned reads to bacterial genomes

###  Indexing

```bash
cd /mnt/e/project/srna/genome/bacteria

bowtie2-build --threads 12 --quiet bacteria.fna bacteria
```

###  Aligning

```bash
cd /mnt/e/project/srna/output

bowtie2 -x ../genome/bacteria/bacteria -q SRR4039757_trimmed.fq.unali.fq --end-to-end --threads 12 -S SRR4039757_trimmed.fq.unali.fq.sam

nohup parallel -j 3 "
    bowtie2 -q {} --end-to-end -x ../genome/bacteria/bacteria \
    --threads 4 -S ./{}.sam
" ::: $(ls SRR1004935*.unali.fq)
```



## Aligning unaligned reads to bacterial tRNAs

I chose tRNAs from 192 bacterial species as our target to check whether reads could align to bacterial tRNAs.

###  Indexing

```bash
cd /mnt/e/project/srna/bac/trna
bowtie-build --threads 12 --quiet bac_trna.fna bac
```

Get .bed file of selected bacterial tRNAs

```bash
cd /mnt/e/project/srna/annotation/bacteria

cat bac.gtf | grep -v '#' | tsv-filter --str-eq 2:RefSeq --str-eq 3:gene --regex '9:tRNA' > bac_trna.gtf
cat bac_trna.gtf | convert2bed --input=gtf > bac_trna.bed
```

```bash
cat bac.gtf | grep -v '#' | tsv-filter --str-eq 2:RefSeq --str-eq 3:gene --regex '9:rRNA' > bac_rrna.gtf
cat bac_rrna.gtf | convert2bed --input=gtf > bac_rrna.bed
```



```bash
cd /mnt/e/project/srna/output

parallel -j 3 "
	samtools sort -@ 4 {1}.sam > {1}.sort.bam
	samtools index {1}.sort.bam
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')
```

mosdepth for counting the coverage

```bash
parallel -j 3 "
	mosdepth -t 4 -b ../annotation/bacteria/bac_trna.bed {} {}
" ::: $(ls *unali.fq.sort.bam)
```

```bash
parallel -j 3 "
	mosdepth -t 4 -b ../annotation/bacteria/bac_rrna.bed {}rrna {}
" ::: $(ls *unali.fq.sort.bam)
```



###  Aligning

```bash
cd /mnt/e/project/srna/sRNA/output/bowtie

nohup parallel -j 3 "
    bowtie -v 0 -a -m1 --best --strata --threads 4 --quiet \
    -x ../../../bac/trna/bac {} --al ../bac_bowtie/{}.bac_ali.fastq \
" ::: $(ls *.unali.fastq)
```

```bash
nohup parallel -j 3 "
    bowtie -v 0 -a -m1 --best --strata --threads 4 --quiet \
    -x ../../../bac/trna/bac {} -S ../bac_bowtie/{}.bac_ali.sam \
" ::: $(ls *.unali.fastq)
```



```bash
cd /mnt/e/project/srna/sRNA/output/bac_bowtie

nohup parallel -j 3 "
    bowtie -v 0 -a -m1 --best --strata --threads 4 --quiet \
    -x ../../../Atha/genome/Atha_rna {} --al ../bac_bowtie/{}_mrna.ali.fastq \
" ::: $(ls *bac_ali.fastq)
```

```bash
nohup parallel -j 3 "
    bowtie -v 0 -a -m1 --best --strata --threads 4 --quiet \
    -x ../../../Atha/genome/Atha_rna {} -S ../bac_bowtie/{}_mrna.ali.sam \
" ::: $(ls *bac_ali.fastq)
```
