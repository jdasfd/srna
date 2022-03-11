#  sRNA Mapping

Here we use *Arabidopsis thaliana* genome as an example to prepare potential tRFs screening protocol.

**Attention: I have mounted WSL2 to drive E, so all the ~ will be replaced by /mnt/e. Restricted by the hardware,  : |**



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

#### convert2bed in bedops

```bash
brew install bedops
```



##  Getting sRNA-seq data

###  Download fastq from NCBI

Using anchr for downloading fastq, which could avoid using fastq-dump (time-wasted).

These information needed for downloading: Experiment (SRX number for each sample), Sample_Name (represent each sample name), Bases (size of the sequencing files). We should put it together in an tsv file format and then anchr could get the ftp automatically.

SraRunTable.txt is metadata with all sequencing file information for one project. We could use SraRunSelector for getting this metadata. Download every SraRunTable manually and then extract information from it.

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

 

##  Using Bowtie2 for reads alignment

Because of Bowtie has been outdated and has problem in downloading it. I decided to adopt Bowtie2, though Bowtie has advantages in mapping short reads.

###  Aligning sRNA-seq data to plant genome

First, I aligned sRNA-seq reads to the *A. tha* genome to validate whether there were reads that could not be aligned to the genome.

#### Indexing

```bash
cd /mnt/e/project/srna/genome/plant/genome

bowtie2-build --threads 12 --quiet Atha.fna Atha
```

#### Aligning

* NOTICE

	> I used HPCC for such a number of files. So the bash script were 12 threads (my own computer), and I changed it into 24 threads when submit to  HPCC.
	>
	> Meanwhile, I ignored the files transmission step using rsync. So the directory path were used as my own computer. You could substitute /mnt/e to ~.

```bash
mkdir -p /mnt/e/project/srna/output/bam/plant
cd /mnt/e/project/srna/trim
```

alignall.sh:

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
bsub -q mpi -n 24 -J aliall -o .. "bash alignall.sh"
```

Then we need another mapping round for the 1 mismatch allowed.

align1mis.sh:

```bash
parallel -j 3 " \
bowtie2 -q {}_trimmed.fq.gz -N 1 \
-x ../genome/plant/Atha/Atha --al-gz ../output/fastq/{}_plantali1mis.fq.gz \
--un-gz ../output/fastq/{}_plantunali.fq.gz --no-unal --threads 4 \
-S ../output/bam/plant/{}_plant1mis.sam \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_trimmed.+gz$//')
```

```bash
bsub -q mpi -n 24 -J ali1mis -o .. "bash align1mis.sh"
```

Extract sequences of only 1 mismatch.

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

unali.sh:

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

1mis.sh:

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

all.sh:

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



## Count the ratio of bacterial sources in all samples

### Ratio of reads aligned to bacteria / all non-plant reads

Count all reads numbers from sort.bam files. The goal of this step is to acquire fraction of the reads aligned to bacteria from all reads. I wrote a shell script to reach the goal.

```bash
cd /mnt/e/project/srna/output/bam/bacteria

bash ../../../script/read_count.sh | tee ../../count/read_count.csv
```

```bash
bsub -q mpi -n 24 -o .. -J count "bash read_count.sh | tee ../../count/read_count.csv"
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
labs(x = " ", y = "bacterial reads / all reads")
ggsave(s, file = "read_count.pdf", width = 7, height = 4)
' read_count.csv
```

### Ratio of tRNA reads / aligned reads for different groups

I classified bacteria into five different types according to bacteria-plant(host) relations. Four categories are: endo/epiphyte, environment, gut and marine.

```bash
cd /mnt/e/project/srna/output/count
mkdir trna rrna mrna all
```

Using bed of rna to extract mapping reads from different RNA regions.

```bash
cd /mnt/e/project/srna/output/bam/bacteria

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
cd /mnt/e/project/srna/output/bam/rna

parallel -j 3 " 
samtools sort -@ 4 {1}.bam > {1}.sort.bam 
samtools index {1}.sort.bam
" ::: $(ls *.bam | perl -p -e 's/\.bam$//')

parallel -j 3 " \
rm {}.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

Use idxstats to count the different chromosome cover. Because there were 190 bacteria included, so the chromosome numbers should be joined to its own name.

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
cd /mnt/e/project/srna/output/bam/bacteria

parallel -j 4 " \
samtools idxstats {}.sort.bam | \
tsv-select -f 1,3 | grep -v '*' | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 \
> ../../count/all/{}.all.tsv \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

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
-f result.trna.tsv -t tRNA -o trna_percent_group.pdf -a trna_percent.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.rrna.tsv -t rRNA -o rrna_percent_group.pdf -a rrna_percent.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.mrna.tsv -t mRNA -o mrna_percent_group.pdf -a mrna_percent.pdf
```



## Frequencies of bacteria among different groups

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
Rscript /mnt/e/project/srna/script/rna_plot.r -f {}.num.tsv -o {}_freq.pdf -t {} -y frequencies \
" ::: $(ls *.tsv | perl -p -e 's/\.num\.tsv//')
```



## Reads TPM among different groups

The coverage of each bed region and the mean depth using coverage/RNA region length were calculated and compared among groups.

```bash
mkdir -p /mnt/e/project/srna/output/cover
cd /mnt/e/project/srna/output/cover
mkdir mrna rrna trna
```

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

Rscript /mnt/e/project/srna/script/rna_plot.r -f ttpm.tsv -o ttpm.pdf -t trna -y TPM
Rscript /mnt/e/project/srna/script/rna_plot.r -f rtpm.tsv -o rtpm.pdf -t rrna -y TPM
Rscript /mnt/e/project/srna/script/rna_plot.r -f mtpm.tsv -o mtpm.pdf -t mrna -y TPM

tsv-select -H -f file,group,catgry,ttpm tpm.tsv | sed '1d' | sed '1ifile\tgroup\tcatgry\tnum' > ttpm.tsv
tsv-select -H -f file,group,catgry,rtpm tpm.tsv | sed '1d' | sed '1ifile\tgroup\tcatgry\tnum' > rtpm.tsv
tsv-select -H -f file,group,catgry,mtpm tpm.tsv | sed '1d' | sed '1ifile\tgroup\tcatgry\tnum' > mtpm.tsv
```



Because of the mosdepth cannot decide output directory. So it was better using mosdepth in the directory you want to save the results.

```bash
cd /mnt/e/project/srna/output/cover/trna

parallel -j 3 " \
mosdepth -t 4 -b ../../../annotation/bacteria/trna.bed \
{} ../../bam/bacteria/{}.sort.bam \
" ::: $(ls ../../bam/bacteria/*.sort.bam | \
perl -p -e 's/\.sort\.bam//' | perl -p -e 's/^.+?\/S/S/')
```

```bash
cd /mnt/e/project/srna/output/cover/rrna

parallel -j 3 " \
mosdepth -t 4 -b ../../../annotation/bacteria/rrna.bed \
{} ../../bam/bacteria/{}.sort.bam \
" ::: $(ls ../../bam/bacteria/*.sort.bam | \
perl -p -e 's/\.sort\.bam//' | perl -p -e 's/^.+?\/S/S/')
```

```bash
cd /mnt/e/project/srna/output/cover/mrna

parallel -j 3 " \
mosdepth -t 4 -b ../../../annotation/bacteria/mrna.bed \
{} ../../bam/bacteria/{}.sort.bam \
" ::: $(ls ../../bam/bacteria/*.sort.bam | \
perl -p -e 's/\.sort\.bam//' | perl -p -e 's/^.+?\/S/S/')
```

```bash
cd /mnt/e/project/srna/output/cover

for dir in `ls`
do
if [[ -d "$dir" ]]
then
cd $dir;
for file in `ls *.per-base.bed.gz`
do
gzip -d ${file}
done
cd ..
fi
done
```



##  Reads cover depth and  position in different RNA region

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

### Mosdepth

Mosdepth for counting the coverage.

```bash
cd /mnt/e/project/srna/output/bam/bacteria

parallel -j 3 " \
mosdepth -t 4 {}_aliall {}_aliall.sort.bam \
" ::: $(ls *_aliall.sort.bam | perl -p -e 's/_aliall\.sort\.bam//')
# -t: threads, it has been said that 4 could reach the max speed

mv *_aliall.per-base.bed* *.txt ../../mosdepth
```

```bash
parallel -j 3 " \
mosdepth -t 4 {}_1mis {}_1mis.sort.bam \
" ::: $(ls *_1mis.sort.bam | perl -p -e 's/_1mis\.sort\.bam//')

mv *_1mis.per-base.bed* *.txt ../../mosdepth
```

```bash
parallel -j 3 " \
mosdepth -t 4 {}_unali {}_unali.sort.bam \
" ::: $(ls *_unali.sort.bam | perl -p -e 's/_unali\.sort\.bam//')

mv *_unali.per-base.bed* *.txt ../../mosdepth
```

Unzip all per-base.bed.gz.

```bash
cd /mnt/e/project/srna/output/mosdepth

gzip -d *.per-base.bed.gz
```

Convert bed format file to runlist files for better manipulating using a perl script. Each covered base will be counted.

```bash
parallel -j 10 " \
cat {}.per-base.bed | perl ../../script/bed2yml.pl > {}.yml \
" ::: $(ls *.per-base.bed | perl -p -e 's/\.per-base\.bed//')
```

### spanr  for tRNA coverage length

Use spanr from wang-q [intspan](https://github.com/wang-q/intspan).

Get all used bacteria genome size in .chr.sizes format.

```bash
mkdir -p /mnt/e/project/srna/output/opt
mkdir -p /mnt/e/project/srna/output/result
cd /mnt/e/project/srna/output/opt

faops size ../../genome/bacteria/bacteria.fna > bacteria.chr.sizes
```

Calculate the coverage.

The .csv file contains 4 columns, chr, chrLength, size and coverage. We need the column 2 ‘chrLength’ (representing genome length) and the column 3 ‘size’ (representing genome covered length). Column 4 coverage: column 3 (coverage length) / column 2 (chrLength), so it does not mean the reads cover.

```bash
cd /mnt/e/project/srna/output/mosdepth

parallel -j 10 " \
spanr stat ../opt/bacteria.chr.sizes {}.yml -o ../opt/{}.csv \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')
```

Extracte tRNA from .gff file.

```bash
cd /mnt/e/project/srna/annotation/bacteria

spanr gff bacteria.gff --tag tRNA > tRNA.yml
# --tag: selected gene name

spanr gff bacteria.gff --tag rRNA > rRNA.yml

spanr gff bacteria.gff --tag CDS > mRNA_all.yml
spanr compare --op diff mRNA_all.yml tRNA.yml -o mRNA1.yml
spanr compare --op diff mRNA1.yml rRNA.yml -o mRNA.yml
rm mRNA_all.yml mRNA1.yml
# mRNA region without tRNA and rRNA
# spanr compare could manipulate aggregation, including intersect (default), union, diff or xor
```

Compare two runlist files for their intersect part.

```bash
cd /mnt/e/project/srna/output/mosdepth

parallel -j 10 " \
spanr compare {}.yml ../../annotation/bacteria/tRNA.yml \
-o ../opt/{}.tRNA.intersect.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')

parallel -j 10 " \
spanr compare --op diff {}.yml ../../annotation/bacteria/tRNA.yml \
-o ../opt/{}.tRNA.diff.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')

parallel -j 10 " \
spanr compare {}.yml ../../annotation/bacteria/rRNA.yml \
-o ../opt/{}.rRNA.intersect.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')

parallel -j 10 " \
spanr compare --op diff {}.yml ../../annotation/bacteria/rRNA.yml \
-o ../opt/{}.rRNA.diff.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')

parallel -j 10 " \
spanr compare {}.yml ../../annotation/bacteria/mRNA.yml \
-o ../opt/{}.mRNA.intersect.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')

parallel -j 10 " \
spanr compare --op diff {}.yml ../../annotation/bacteria/mRNA.yml \
-o ../opt/{}.mRNA.diff.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')
```

Get the same .csv results including tRNA length and tRNA covered length. rRNA and mRNA regions were used as control.

```bash
cd /mnt/e/project/srna/output/opt

spanr stat bacteria.chr.sizes ../../annotation/bacteria/tRNA.yml -o tRNA.csv
spanr stat bacteria.chr.sizes ../../annotation/bacteria/rRNA.yml -o rRNA.csv
spanr stat bacteria.chr.sizes ../../annotation/bacteria/mRNA.yml -o mRNA.csv
```

```bash
cd /mnt/e/project/srna/output/opt

parallel -j 10 " \
spanr stat bacteria.chr.sizes {}.tRNA.intersect.yml \
-o {}.tRNA.intersect.csv \
" ::: $(ls *.tRNA.intersect.yml | perl -p -e 's/\.tRNA.+yml$//')

parallel -j 10 " \
spanr stat bacteria.chr.sizes {}.tRNA.diff.yml \
-o {}.tRNA.diff.csv \
" ::: $(ls *.tRNA.diff.yml | perl -p -e 's/\.tRNA.+yml$//')

parallel -j 10 " \
spanr stat bacteria.chr.sizes {}.rRNA.intersect.yml \
-o {}.rRNA.intersect.csv \
" ::: $(ls *.rRNA.intersect.yml | perl -p -e 's/\.rRNA.+yml$//')

parallel -j 10 " \
spanr stat bacteria.chr.sizes {}.rRNA.diff.yml \
-o {}.rRNA.diff.csv \
" ::: $(ls *.rRNA.diff.yml | perl -p -e 's/\.rRNA.+yml$//')

parallel -j 10 " \
spanr stat bacteria.chr.sizes {}.mRNA.intersect.yml \
-o {}.mRNA.intersect.csv \
" ::: $(ls *.mRNA.intersect.yml | perl -p -e 's/\.mRNA.+yml$//')

parallel -j 10 " \
spanr stat bacteria.chr.sizes {}.mRNA.diff.yml \
-o {}.mRNA.diff.csv \
" ::: $(ls *.mRNA.diff.yml | perl -p -e 's/\.mRNA.+yml$//')
```

Convert .csv to .tsv format.

```bash
parallel -j 10 " \
cat {}.csv | csv2tsv -H > {}.tsv \
" ::: $(ls *.csv | perl -p -e 's/\.csv$//')

rm *.csv *.yml
```

Use tsv-utils join tsv together.

```bash
cd /mnt/e/project/srna/output/opt

parallel -j 3 " \
cat {}.tRNA.diff.tsv | cut -f 1,2,3 | \
tsv-join -H --filter-file tRNA.tsv --key-fields 1 --append-fields 3 | \
tsv-join -H --filter-file {}.tRNA.intersect.tsv --key-fields 1 --append-fields 3 | \
tsv-join --filter-file ../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 6 --sum 2 --sum 3 --sum 4 --sum 5 | \
tsv-join --filter-file ../../name.tsv --key-fields 1 --append-fields 2 \
> ../result/{}.tRNA.tsv \
" ::: $(ls *.tRNA.*.tsv | perl -p -e 's/\..*tsv$//' | uniq)

parallel -j 3 " \
cat {}.rRNA.diff.tsv | cut -f 1,2,3 | \
tsv-join -H --filter-file rRNA.tsv --key-fields 1 --append-fields 3 | \
tsv-join -H --filter-file {}.rRNA.intersect.tsv --key-fields 1 --append-fields 3 | \
tsv-join --filter-file ../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 6 --sum 2 --sum 3 --sum 4 --sum 5 | \
tsv-join --filter-file ../../name.tsv --key-fields 1 --append-fields 2 \
> ../result/{}.rRNA.tsv \
" ::: $(ls *.rRNA.*.tsv | perl -p -e 's/\..*tsv$//' | uniq)

parallel -j 3 " \
cat {}.mRNA.diff.tsv | cut -f 1,2,3 | \
tsv-join -H --filter-file mRNA.tsv --key-fields 1 --append-fields 3 | \
tsv-join -H --filter-file {}.mRNA.intersect.tsv --key-fields 1 --append-fields 3 | \
tsv-join --filter-file ../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 6 --sum 2 --sum 3 --sum 4 --sum 5 | \
tsv-join --filter-file ../../name.tsv --key-fields 1 --append-fields 2 \
> ../result/{}.mRNA.tsv \
" ::: $(ls *.mRNA.*.tsv | perl -p -e 's/\..*tsv$//' | uniq)
```

###  Chi-square test

Use chi-square tests performing chi-squared contingency  interspecies.

```bash
mkdir -p /mnt/e/project/srna/output/chi/inter
cd /mnt/e/project/srna/output/result

for tsv in `ls *.tsv`
do
cat $tsv | perl ../../script/interspecie.pl | \
parallel --colsep '\t' -j 1 -k '
    echo "==> {1}"
    Rscript -e "
        x <- matrix(c({2},{4},{3},{5}), nrow=2)
        x
        chisq.test(x)
        "
' > ../chi/inter/$tsv.inter.txt
done
```

Use a simple script for better looking.

```bash
cd /mnt/e/project/srna/output/chi/inter

for file in `ls *.inter.txt | perl -p -e 's/\.tsv.+txt$//'`
do
echo "${file}" >> ../result.tsv
cat ${file}.tsv.inter.txt | perl ../../../script/square.pl >> ../result.tsv
done
```

###  Extract tRNA reads

```bash
mkdir -p /mnt/e/project/srna/output/target/trna
mkdir -p /mnt/e/project/srna/output/target/rrna
```

```bash
cd /mnt/e/project/srna/output/bam/rna

parallel -j 3 " \
samtools fasta -@ 4 -F 4 {}.trna.sort.bam > ../../target/trna/{}.fasta \
" ::: $(ls *_1mis.trna.sort.bam | perl -p -e 's/\.trna.+bam$//')
```

```bash
parallel -j 3 " \
samtools fasta -@ 4 -F 4 {}.rrna.sort.bam > ../../target/rrna/{}.fasta \
" ::: $(ls *_1mis.rrna.sort.bam | perl -p -e 's/\.rrna.+bam$//')
```

###  tRNA reads statistical information

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



###  psRNATarget

Use psRNATarget websites for miRNA target prediction. Download files from the website.



###  ClusterProfiler

```bash
cd /mnt/e/project/srna/output/kegg

cat psRNATargetJob-trna.txt | sed 1,2d | grep 'Cleavage' | cut -f 2 | cut -d '.' -f 1 | sed 1iTAIR > gene_trna.tsv
cat psRNATargetJob-rrna.txt | sed 1,2d | grep 'Cleavage' | cut -f 2 | cut -d '.' -f 1 | sed 1iTAIR > gene_rrna.tsv
```

```R
library(readr)
library(clusterProfiler)
library(org.At.tair.db) #Arabidopsis thaliana annotation for clusterProfiler
trna <- read_tsv("gene_trna.tsv")
rrna <- read_tsv("gene_rrna.tsv")

tgo <- enrichGO(gene = tgene$go_id,
                OrgDb = org.At.tair.db,
                ont = "ALL",
                keyType = 'GO',
                pvalueCutoff = 0.05,
                readable = TRUE)
rgo <- enrichGO(gene = rrna$TAIR, OrgDb = org.At.tair.db, ont = "ALL", keyType = 'TAIR', pvalueCutoff = 1)

rrna <- grrna$TAIR
# extract values from data.frame.

library(biomaRt)
listMarts(host = "plants.ensembl.org")
# the default listMarts() connected to ensembl.org for animals, so we used host to connect plant.ensembl.org for our databases
#             biomart                      version
# 1       plants_mart      Ensembl Plants Genes 51
# 2 plants_variations Ensembl Plants Variations 51

mart <- useMart("plants_mart", host = "plants.ensembl.org")
listDatasets(mart)
# the first is biomart, we used plants_mart in the previous step (row 1)

mart <- useMart("plants_mart", "athaliana_eg_gene", host = "plants.ensembl.org")
listAttributes(mart)
# view all output ID formats

tgene <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "go_id", "tair_locus"),
               filters = 'tair_locus',
               values = locus,
               mart = mart)
rgene <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "go_id"), 
               filters = 'tair_locus',
               values = rrna$TAIR,
               mart = mart)
# attributes: output ID formats, filters: our ID formats, values: our data, mart: chosen mart
```



----------------------------------------------------

Use goodness-of-fit tests performing genome covered length and tRNA covered length consistent with expectation.

```bash
mkdir -p /mnt/e/project/srna/output/chi/intra
cd /mnt/e/project/srna/output/result

for tsv in `ls *.tsv`
do
cat $tsv | parallel --colsep '\t' -j 1 -k '
        echo "==> {1}"
        Rscript -e "
            x <- c({5}, {3})
            a <- {4}/{2}
            b <- 1-({4}/{2})
            chisq.test(x, p = c(a, b))
        "
' > ../chi/intra/$tsv.intra.txt
done
```

Use a simple script for getting chi-square results.

```bash
cd /mnt/e/project/srna/output/chi/inter

parallel -j 3 " \
cat {}.tsv.chi-square.txt | perl ../../../script/square.pl > {}.chi.tsv \
" ::: $(ls *.tsv.chi-square.txt | perl -p -e 's/\.tsv.*txt$//')
```

Join tsv together for better analysis.

```bash
parallel -j 3 " \
cat {}.chi.tsv | tsv-join --filter-file name.tsv -k 1 --append-fields 2,3 > {}.tsv \
" ::: $(ls *.chi.tsv | perl -p -e 's/\.chi\.tsv//')

rm *.chi.tsv
```

```R
bar <- ggplot(data = name, mapping = aes(x = 'category', y = count, fill = cate)) +
geom_bar(stat = 'identity', position = 'stack') +
coord_polar(theta = 'y') +
labs(x = '', y = '', title = '') +

```



```bash

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

