# Preparations for the sRNA analytic process

In this process, I prepared to analyse small RNA-seq files in *A. thaliana*. There were few things should be specified. I recorded below.

##  Prepare

### *Arabidopsis thaliana*

#### *A. thaliana* genome

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

#### *A. thaliana* annotation

Download gff3 format annotation.

```bash
cd /mnt/e/project/srna/
mkdir -p annotation/plant/Atha
cd annotation/plant/Atha

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz

gzip -d GCF_000001735.4_TAIR10.1_genomic.gff.gz
mv GCF_000001735.4_TAIR10.1_genomic.gff Atha.gff
```

### Bacteria

I used ASSEMBLY in [bacteria_ar.md](https://github.com/wang-q/withncbi/blob/master/pop/bacteria_ar.md). Please go check the markdown to get the database.

#### Bacteria annotation

There is a demand that separating RNA GFF from annotation. tRNA, rRNA and mRNA regions were three different types.

```bash
cat bacteria.gff | grep -v '#' | tsv-filter --str-eq 3:gene --iregex 9:tRNA > bac_trna.gff
cat bacteria.gff | grep -v '#' | tsv-filter --str-eq 3:gene --iregex 9:rRNA > bac_rrna.gff
cat bacteria.gff | grep -v '#' | tsv-filter --str-eq 3:CDS --not-iregex 9:tRNA --not-iregex 9:rRNA > bac_mrna.gff

# All CDS contain tRNAs and rRNAs, which should be removed 
```

```bash
cat bac_trna.gff | convert2bed --input=gff --output=bed > trna.bed
cat bac_rrna.gff | convert2bed --input=gff --output=bed > rrna.bed
cat bac_mrna.gff | convert2bed --input=gff --output=bed > mrna.bed
```

###  Biotools in protocol

#### sratoolkit (not recommended)

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

####  Samtools

The basic SAM/BAM file manipulating tools.

```bash
brew install samtools
```

#### Fastqc

Fastqc for checking sequencing quality.

```bash
brew install fastqc
```

#### intspan and anchr

[intspan](https://github.com/wang-q/intspan) and [anchr](https://github.com/wang-q/anchr) were written by my lab professor [Qiang Wang](https://www.github.com/wang-q).

Go check and install.

#### Bowtie2

```bash
brew install bowtie2
```

#### Seqkit

Using seqkit to extract reads.

```bash
brew install seqkit
```

#### mosdepth

```bash
brew install mosdepth
```

#### tsv-utils

tsv-utils is an really powerful toolbox for manipulating .tsv format files. It is highly recommended that you mastered it.

Make sure you completed software installation [wang-q/dotfiles](https://github.com/wang-q/dotfiles) and `brew tap wang-q/tap` could be adopted. 

```bash
brew install wang-q/tap/tsv-utils
```

#### convert2bed in bedops

```bash
brew install bedops
```


##  Getting sRNA-seq data

###  Download fastq from NCBI

Using anchr for downloading fastq, which could avoid using fastq-dump (time-wasted).

These information needed for downloading: Experiment (SRX number for each sample), Sample_Name (represent each sample name), Bases (size of the sequencing files). We should put it together in an tsv file format and then anchr could get the ftp automatically.

SraRunTable.txt is metadata with all sequencing file information for one project. We could use SraRunSelector for getting this metadata. Download every SraRunTable manually and then extract information from it.

Only ecotype of Col-0 could be kept in this protocol, according to our reference genome was *A. tha*.

```bash
mkdir -p /mnt/e/project/srna/ena/thale
cd /mnt/e/project/srna/ena
# <download metadata file here>

cat SraRunTable.txt | mlr --icsv --otsv cat | \
tsv-select -H -f Experiment,"Sample\ Name",Bases > ./thale/SraRunTable.tsv
# rm and get next

cat SraRunTable.txt | mlr --icsv --otsv cat | \
tsv-select -H -f Experiment,"Sample\ Name",Bases | sed '1d' >> ./thale/SraRunTable.tsv
# continue this until all SraRunTable had been manipulated
# if there are RNA-seq and ncRNA/sRNA-seq, you could use tsv-filter

cd thale
cat SraRunTable.tsv | sed '1 s/^/#/' | \
keep-header -- tsv-sort -k2,2 -k3,3nr | \
tsv-uniq -H -f "Sample\ Name" --max 1 | mlr --itsv --ocsv cat > source.csv

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv
# viewing md format of each sequencing file

aria2c -j 4 -x 4 -s 2 --file-allocation=none -c -i ena_info.ftp.txt
# aria2c could be used, although it was better using aspera (shell need to be modified if using aspera)

md5sum --check ena_info.md5.txt
# check if there was any mistake during downloading
```

I just classified different sequencing file according to plant in manual. Thale, rice, maize ... (needed to be completed).

All sRNA-seq files using *A. thaliana* were added up to 244 files.

###  Quality control

Use fastqc to check the quality of sequencing for each fastq file.

 ```bash
mkdir -p /mnt/e/project/srna/fastqc/raw
cd /mnt/e/project/srna/ena/thale

fastqc -t 12 --quiet -o ../../fastqc/raw *.gz
# -t: threads number
 ```

Use trim galore for automatic adapters cut. Trim galore is a perl wrapper for fastqc and cutadapt. It could be automatically or manually remove adapters.

```bash
cd /mnt/e/project/srna/ena/thale

parallel -j 4 " \
trim_galore --phred33 -j 3 \
--length 17 --output_dir ../../trim {} \
" ::: $(ls *.fastq.gz)
# --small_rna: automatically remove adapters for small rna sequencing. If using this, trim_galore would only trim small RNA adapters
```

Check the quality of sequencing files after trim.

```bash
mkdir -p /mnt/e/project/srna/fastqc/after_trim
cd /mnt/e/project/srna/trim

fastqc -t 12 --quiet -o ../fastqc/after_trim *.gz
```

After I used fastqc, there were few fail in per base sequence content. It was common in this item because of fastqc not applicable for the sRNA-seq task in this step. But we could see that adapters were all removed perfectly after trimming. Now we could use sequencing files for analyze.

### Convert fastq to fasta

Using faops to convert all fastq to fasta format because of the non-essential use of quality lines.

```bash
parallel -j 20 " \
faops filter -l 0 {}.fq.gz {}.fa \
" ::: $(ls *.fq.gz | perl -p -e 's/\.fq\.gz$//')
```