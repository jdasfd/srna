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
cd /mnt/e/project/srna/srr
parallel -j 12 "fastq-dump {}" ::: $(ls *.sra)
rm *.sra
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



##  Using Bowtie2

Because of Bowtie has been outdated and has problem in downloading it. I decided to adopt Bowtie2, though Bowtie has advantages in mapping short reads.

###  Aligning sRNA-seq data to plant genome

First, I aligned sRNA-seq reads to the *A. tha* genome to validate whether there were reads that could not be aligned to the genome.

#### Indexing

```bash
cd /mnt/e/project/srna/genome/plant/genome

bowtie2-build --threads 12 --quiet Atha.fna Atha
```

####  Aligning

```bash
cd /mnt/e/project/srna/trim

parallel -j 3 "
    bowtie2 -q {}_trimmed.fq --local -x ../genome/plant/Atha/Atha \
    --un ../output/{}_unali.fq --threads 4 -S ../output/{}_plant.sam
" ::: $(ls SRR1004935*.fq | perl -p -e 's/_trimmed\.fq//')
```

After run this, bowtie2 would show alignment results on screen

```bash
23131748 reads; of these:
  23131748 (100.00%) were unpaired; of these:
    9194122 (39.75%) aligned 0 times
    9421082 (40.73%) aligned exactly 1 time
    4516544 (19.53%) aligned >1 times
60.25% overall alignment rate
23119819 reads; of these:
  23119819 (100.00%) were unpaired; of these:
    8802807 (38.07%) aligned 0 times
    10349574 (44.76%) aligned exactly 1 time
    3967438 (17.16%) aligned >1 times
61.93% overall alignment rate
21565629 reads; of these:
  21565629 (100.00%) were unpaired; of these:
    8696299 (40.32%) aligned 0 times
    7765072 (36.01%) aligned exactly 1 time
    5104258 (23.67%) aligned >1 times
59.68% overall alignment rate
22224601 reads; of these:
  22224601 (100.00%) were unpaired; of these:
    8400742 (37.80%) aligned 0 times
    8933235 (40.20%) aligned exactly 1 time
    4890624 (22.01%) aligned >1 times
62.20% overall alignment rate
22719863 reads; of these:
  22719863 (100.00%) were unpaired; of these:
    8497773 (37.40%) aligned 0 times
    8664376 (38.14%) aligned exactly 1 time
    5557714 (24.46%) aligned >1 times
62.60% overall alignment rate
```



###  Aligning unaligned reads to bacterial genomes

We chose 365 bacterial genomes from 191 species as our target bacteria.

####  Indexing

```bash
cd /mnt/e/project/srna/genome/bacteria

bowtie2-build --threads 12 --quiet bacteria.fna bacteria
```

####  Aligning

```bash
cd /mnt/e/project/srna/output

parallel -j 3 "
    bowtie2 -q {}_unali.fq --local -x ../genome/bacteria/bacteria \
    --threads 4 -S ./{}_bac.sam
" ::: $(ls SRR1004935*_unali.fq | perl -p -e 's/_unali\.fq$//')
```

Alignment results:

```bash
8802807 reads; of these:
  8802807 (100.00%) were unpaired; of these:
    8727037 (99.14%) aligned 0 times
    6225 (0.07%) aligned exactly 1 time
    69545 (0.79%) aligned >1 times
0.86% overall alignment rate
8696299 reads; of these:
  8696299 (100.00%) were unpaired; of these:
    8669356 (99.69%) aligned 0 times
    7280 (0.08%) aligned exactly 1 time
    19663 (0.23%) aligned >1 times
0.31% overall alignment rate
9194122 reads; of these:
  9194122 (100.00%) were unpaired; of these:
    9163807 (99.67%) aligned 0 times
    10560 (0.11%) aligned exactly 1 time
    19755 (0.21%) aligned >1 times
0.33% overall alignment rate
8400742 reads; of these:
  8400742 (100.00%) were unpaired; of these:
    8390361 (99.88%) aligned 0 times
    1805 (0.02%) aligned exactly 1 time
    8576 (0.10%) aligned >1 times
0.12% overall alignment rate
8497773 reads; of these:
  8497773 (100.00%) were unpaired; of these:
    8448696 (99.42%) aligned 0 times
    18110 (0.21%) aligned exactly 1 time
    30967 (0.36%) aligned >1 times
0.58% overall alignment rate
```

Convert sam to bam file for minimum storage stress.

```bash
cd /mnt/e/project/srna/output

parallel -j 3 "
	samtools sort -@ 4 {1}.sam > {1}.sort.bam
	samtools index {1}.sort.bam
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')

rm *.sam
```

### mosdepth for counting the coverage

```bash
parallel -j 3 "
	mosdepth -t 4 -b ../annotation/bacteria/bac_trna.bed {}_trna {}_bac.sort.bam
" ::: $(ls *_bac.sort.bam | perl -p -e 's/_bac\.sort\.bam//')
# -t: threads, it has been said that 4 could reach the max speed 
```

```bash
cat SRR10049355_trna.mosdepth.summary.txt | perl ../species/chi.pl > ../species/result.tsv
```

```bash
cd /mnt/e/project/srna/species

cat result.tsv |
    parallel --colsep '\t' -j 1 -k '
        echo "==> {1}"
        Rscript -e "
            x <- matrix(c({2}, {4}, {3}, {5}), nrow=2)
            x
            chisq.test(x)
        "
    '
```

