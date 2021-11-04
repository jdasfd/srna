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
parallel -j 4 " \
trim_galore --phred33 --small_rna --output_dir ../trim {} \
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

parallel -j 3 " \
bowtie2 -q {}_trimmed.fq -N 0 -x ../genome/plant/Atha/Atha \
--al ../output/fastq/{}_aliall.fq --no-unal --threads 4 -S ../output/bam/{}_plantall.sam \
" ::: $(ls SRR1004935*.fq | perl -p -e 's/_trimmed\.fq//')
```

After run the above command, we would get results (none mismatch allowed).

```bash
21565629 reads; of these:
  21565629 (100.00%) were unpaired; of these:
    1059845 (4.91%) aligned 0 times
    8465079 (39.25%) aligned exactly 1 time
    12040705 (55.83%) aligned >1 times
95.09% overall alignment rate
23131748 reads; of these:
  23131748 (100.00%) were unpaired; of these:
    1164936 (5.04%) aligned 0 times
    10693063 (46.23%) aligned exactly 1 time
    11273749 (48.74%) aligned >1 times
94.96% overall alignment rate
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

Then we need another mapping round for the 1 mismatch allowed.

```bash
parallel -j 3 " \
bowtie2 -q {}_trimmed.fq -N 1 -x ../genome/plant/Atha/Atha \
--al ../output/fastq/{}_ali1mis.fq --un ../output/fastq/{}_unali.fq --no-unal --threads 4 -S ../output/bam/{}_plant1mis.sam \
" ::: $(ls SRR1004935*.fq | perl -p -e 's/_trimmed\.fq//')
```

After run this, bowtie2 would show alignment results on screen.

```bash
23131748 reads; of these:
  23131748 (100.00%) were unpaired; of these:
    1067707 (4.62%) aligned 0 times
    10500848 (45.40%) aligned exactly 1 time
    11563193 (49.99%) aligned >1 times
95.38% overall alignment rate
23119819 reads; of these:
  23119819 (100.00%) were unpaired; of these:
    878194 (3.80%) aligned 0 times
    11431813 (49.45%) aligned exactly 1 time
    10809812 (46.76%) aligned >1 times
96.20% overall alignment rate
21565629 reads; of these:
  21565629 (100.00%) were unpaired; of these:
    969920 (4.50%) aligned 0 times
    8299226 (38.48%) aligned exactly 1 time
    12296483 (57.02%) aligned >1 times
95.50% overall alignment rate
22719863 reads; of these:
  22719863 (100.00%) were unpaired; of these:
    1207658 (5.32%) aligned 0 times
    9077386 (39.95%) aligned exactly 1 time
    12434819 (54.73%) aligned >1 times
94.68% overall alignment rate
22224601 reads; of these:
  22224601 (100.00%) were unpaired; of these:
    633246 (2.85%) aligned 0 times
    9225018 (41.51%) aligned exactly 1 time
    12366337 (55.64%) aligned >1 times
97.15% overall alignment rate
```

#### Extract sequences of only 1 mismatch.

```bash
cd /mnt/e/project/srna/output/fastq

parallel -j 3 " \
seqkit grep --quiet -n -v -f <(seqkit seq -n {}_aliall.fq) {}_ali1mis.fq -o {}_1mis.fq \
" ::: $(ls SRR1004935*.fq | perl -p -e 's/_.+\.fq$//' | uniq)

rm *_ali1mis.fq
```



###  Aligning different reads to bacterial genomes

We chose 365 bacterial genomes from 191 species as our target bacteria. From the previous step, we split reads to 3 types: matched without any mistake, 1 mismatch allowed and unaligned reads (without any match on plant genome). 

####  Indexing

```bash
cd /mnt/e/project/srna/genome/bacteria

bowtie2-build --threads 12 --quiet bacteria.fna bacteria
```

####  Aligning

Aligning unaligned reads to bacteria species.

```bash
cd /mnt/e/project/srna/output/fastq

parallel -j 3 " \
bowtie2 -q {}_unali.fq -x ../../genome/bacteria/bacteria --threads 4 -S ../bam/{}_unali.sam \
" ::: $(ls SRR1004935*_unali.fq | perl -p -e 's/_unali\.fq$//')
```

Alignment results:

```bash
878194 reads; of these:
  878194 (100.00%) were unpaired; of these:
    656447 (74.75%) aligned 0 times
    36265 (4.13%) aligned exactly 1 time
    185482 (21.12%) aligned >1 times
25.25% overall alignment rate
969920 reads; of these:
  969920 (100.00%) were unpaired; of these:
    850850 (87.72%) aligned 0 times
    47548 (4.90%) aligned exactly 1 time
    71522 (7.37%) aligned >1 times
12.28% overall alignment rate
1067707 reads; of these:
  1067707 (100.00%) were unpaired; of these:
    928017 (86.92%) aligned 0 times
    52590 (4.93%) aligned exactly 1 time
    87100 (8.16%) aligned >1 times
13.08% overall alignment rate
633246 reads; of these:
  633246 (100.00%) were unpaired; of these:
    556391 (87.86%) aligned 0 times
    31810 (5.02%) aligned exactly 1 time
    45045 (7.11%) aligned >1 times
12.14% overall alignment rate
1207658 reads; of these:
  1207658 (100.00%) were unpaired; of these:
    1037552 (85.91%) aligned 0 times
    69296 (5.74%) aligned exactly 1 time
    100810 (8.35%) aligned >1 times
14.09% overall alignment rate
```

Aligning 1 mismatch allowed reads to bacteria species.

```bash
parallel -j 3 " \
bowtie2 -q {}_1mis.fq -x ../../genome/bacteria/bacteria --threads 4 -S ../bam/{}_1mis.sam \
" ::: $(ls SRR1004935*_1mis.fq | perl -p -e 's/_1mis\.fq$//')
```

```bash
78576 reads; of these:
  78576 (100.00%) were unpaired; of these:
    75892 (96.58%) aligned 0 times
    478 (0.61%) aligned exactly 1 time
    2206 (2.81%) aligned >1 times
3.42% overall alignment rate
89955 reads; of these:
  89955 (100.00%) were unpaired; of these:
    87678 (97.47%) aligned 0 times
    980 (1.09%) aligned exactly 1 time
    1297 (1.44%) aligned >1 times
2.53% overall alignment rate
97262 reads; of these:
  97262 (100.00%) were unpaired; of these:
    94600 (97.26%) aligned 0 times
    995 (1.02%) aligned exactly 1 time
    1667 (1.71%) aligned >1 times
2.74% overall alignment rate
83126 reads; of these:
  83126 (100.00%) were unpaired; of these:
    81217 (97.70%) aligned 0 times
    833 (1.00%) aligned exactly 1 time
    1076 (1.29%) aligned >1 times
2.30% overall alignment rate
99894 reads; of these:
  99894 (100.00%) were unpaired; of these:
    97139 (97.24%) aligned 0 times
    1181 (1.18%) aligned exactly 1 time
    1574 (1.58%) aligned >1 times
2.76% overall alignment rate
```

Aligning perfectly matched reads to bacteria species.

```bash
parallel -j 3 " \
bowtie2 -q {}_aliall.fq -x ../../genome/bacteria/bacteria --threads 4 -S ../bam/{}_aliall.sam \
" ::: $(ls SRR1004935*_aliall.fq | perl -p -e 's/_aliall\.fq$//')
```

```bash
20505784 reads; of these:
  20505784 (100.00%) were unpaired; of these:
    19612390 (95.64%) aligned 0 times
    429307 (2.09%) aligned exactly 1 time
    464087 (2.26%) aligned >1 times
4.36% overall alignment rate
22163087 reads; of these:
  22163087 (100.00%) were unpaired; of these:
    21471794 (96.88%) aligned 0 times
    404653 (1.83%) aligned exactly 1 time
    286640 (1.29%) aligned >1 times
3.12% overall alignment rate
21966812 reads; of these:
  21966812 (100.00%) were unpaired; of these:
    21269151 (96.82%) aligned 0 times
    366359 (1.67%) aligned exactly 1 time
    331302 (1.51%) aligned >1 times
3.18% overall alignment rate
21412352 reads; of these:
  21412352 (100.00%) were unpaired; of these:
    20676114 (96.56%) aligned 0 times
    361101 (1.69%) aligned exactly 1 time
    375137 (1.75%) aligned >1 times
3.44% overall alignment rate
21508272 reads; of these:
  21508272 (100.00%) were unpaired; of these:
    20748319 (96.47%) aligned 0 times
    376299 (1.75%) aligned exactly 1 time
    383654 (1.78%) aligned >1 times
3.53% overall alignment rate
```

Convert sam to bam file for minimum storage stress.

```bash
cd /mnt/e/project/srna/output/bam

parallel -j 3 " 
samtools sort -@ 4 {1}.sam > {1}.sort.bam 
samtools index {1}.sort.bam
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')

rm *.sam
```



### Get .bed file of selected bacterial tRNAs

```bash
cd /mnt/e/project/srna/annotation/bacteria

cat bac.gtf | grep -v '#' | tsv-filter --str-eq 2:RefSeq --str-eq 3:gene --regex '9:tRNA' > bac_trna.gtf
cat bac_trna.gtf | convert2bed --input=gtf > bac_trna.bed

cat bac.gtf | grep -v '#' | tsv-filter --str-eq 2:RefSeq --str-eq 3:gene --ff-lt 4:5 > bac_gene.gtf
cat bac_gene.gtf | convert2bed --input=gtf > bac_gene.bed
```



### Counting the coverage

####  Mosdepth

Mosdepth for counting the coverage.

```bash
cd /mnt/e/project/srna/output/bam

parallel -j 3 " mosdepth -t 4 {}_aliall {}_aliall.sort.bam " ::: $(ls *_aliall.sort.bam | perl -p -e 's/_aliall\.sort\.bam//')
# -t: threads, it has been said that 4 could reach the max speed

mv *_aliall.per-base.bed* *.txt ../mosdepth
```

```bash
parallel -j 3 " mosdepth -t 4 {}_1mis {}_1mis.sort.bam " ::: $(ls *_1mis.sort.bam | perl -p -e 's/_1mis\.sort\.bam//')

mv *_1mis.per-base.bed* *.txt ../mosdepth
```

```bash
parallel -j 3 " mosdepth -t 4 {}_unali {}_unali.sort.bam " ::: $(ls *_unali.sort.bam | perl -p -e 's/_unali\.sort\.bam//')

mv *_unali.per-base.bed* *.txt ../mosdepth
```

Unzip all per-base.bed.gz.

```bash
gzip -d *.per-base.bed.gz
```

Convert bed format file to runlist files for better manipulating using a perl script.

```bash
parallel -j 3 " \
cat {}.per-base.bed | perl ../../script/bed2yml.pl > {}.yml \
" ::: $(ls *.per-base.bed | perl -p -e 's/\.per-base\.bed//')
```

####  spanr

Use spanr from wang-q [intspan](https://github.com/wang-q/intspan).

Get all used bacteria genome size in .chr.sizes format.

```bash
cd /mnt/e/project/srna/output/mosdepth

faops size ../../genome/bacteria/bacteria.fna > bacteria.chr.sizes
```

Calculate the coverage.

The .csv file contains 4 columns, chr, chrLength, size and coverage. We need the column 2 ‘chrLength’ (representing genome length) and the column 3 ‘size’ (representing genome covered length).

```bash
parallel -j 3 " \
spanr stat bacteria.chr.sizes {}.yml -o ../result/{}.csv \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')
```

Extracte tRNA from .gff file.

```bash
cd /mnt/e/project/srna/annotation/bacteria

spanr gff bacteria.gff --tag tRNA > tRNA.yml
# --tag: selected gene name

spanr gff bacteria.gff --tag rRNA > rRNA.yml
```

Compare two runlist files for their intersect part.

```bash
cd /mnt/e/project/srna/output/mosdepth

parallel -j 3 " \
spanr compare {}.yml ../../annotation/bacteria/tRNA.yml -o {}.intersect.yml \
" ::: $(ls *.yml | perl -p -e 's/\.yml//')
# spanr compare could manipulate aggregation, including intersect (default), union, diff or xor
```

Get the same .csv results including tRNA length and tRNA covered length

```bash
spanr stat bacteria.chr.sizes ../../annotation/bacteria/tRNA.yml -o ../result/tRNA.csv
```

```bash
cd /mnt/e/project/srna/output/mosdepth

parallel -j 3 " \
spanr stat bacteria.chr.sizes {}.intersect.yml -o ../result/{}.intersect.csv \
" ::: $(ls *.intersect.yml | perl -p -e 's/\.inter.+yml$//')
```

Convert .csv to .tsv format.

```bash
cd /mnt/e/project/srna/output/result

parallel -j 3 " \
cat {}.csv | csv2tsv -H > {}.tsv \
" ::: $(ls *.csv | perl -p -e 's/\.csv$//')
```

Use tsv-utils join tsv together.

```bash
parallel -j 3 " \
cat {}.tsv | cut -f 1,2,3 | \
tsv-join -H --filter-file tRNA.tsv --key-fields 1 --append-fields 3 | \
tsv-join -H --filter-file {}.intersect.tsv --key-fields 1 --append-fields 3 | \
sed '1d' > {}.result.tsv \
" ::: $(ls *.tsv | perl -p -e 's/\..*tsv$//' | uniq | grep -v 'tRNA')
```

###  Chi-square test

Use chi-square test for calculating whether there are differences between genome covered length and tRNA covered length. 

```bash
mkdir /mnt/e/project/srna/output/chi
cd /mnt/e/project/srna/output/result

for tsv in `ls *.result.tsv`
do
cat $tsv |
    parallel --colsep '\t' -j 1 -k '
        echo "==> {1}"
        Rscript -e "
            x <- matrix(c({2}, {4}, {3}, {5}), nrow=2)
            x
            chisq.test(x)
        "
    ' > ../chi/$tsv.chi-square.txt
done
```

Get the name list.

```bash
cd /mnt/e/project/srna

perl script/name.pl name.txt > output/chi/name.tsv
```

Use a simple script for getting chi-square results.

```bash
cd /mnt/e/project/srna/output/chi

parallel -j 3 " \
cat {}.result.tsv.chi-square.txt | perl ../../script/square.pl > {}.chi.tsv \
" ::: $(ls *.result.tsv.chi-square.txt | perl -p -e 's/\.result.+txt$//')
```

Join tsv together for better analysis.

```bash
parallel -j 3 " \
cat {}.chi.tsv | tsv-join --filter-file name.tsv -k 1 --append-fields 2,3 > {}.result.tsv \
" ::: $(ls *.chi.tsv | perl -p -e 's/\.chi\.tsv//')
```



###  rRNA

```bash
parallel -j 3 " \
spanr compare {}.yml ../../annotation/bacteria/rRNA.yml -o {}.rRNA.intersect.yml \
" ::: $(ls *.yml | grep -v 'intersect' | perl -p -e 's/\.yml//')
```

```bash
spanr stat bacteria.chr.sizes ../../annotation/bacteria/rRNA.yml -o ../result/rRNA.csv
```

```bash
cd /mnt/e/project/srna/output/mosdepth

parallel -j 3 " \
spanr stat bacteria.chr.sizes {}.intersect.yml -o ../result/{}.intersect.csv \
" ::: $(ls *.rRNA.intersect.yml | perl -p -e 's/\.inter.+yml$//')
```

```bash
cd /mnt/e/project/srna/output/result

parallel -j 3 " \
cat {}.csv | csv2tsv -H > {}.tsv \
" ::: $(ls *.csv | perl -p -e 's/\.csv$//')
```

```bash
cd /mnt/e/project/srna/output/result

parallel -j 3 " \
cat {}.tsv | cut -f 1,2,3 | \
tsv-join -H --filter-file rRNA.tsv --key-fields 1 --append-fields 3 | \
tsv-join -H --filter-file {}.rRNA.intersect.tsv --key-fields 1 --append-fields 3 | \
sed '1d' > {}.rRNA.result.tsv \
" ::: $(ls *.rRNA.intersect.tsv | perl -p -e 's/\.rRNA.*tsv$//')
```

```bash
cd /mnt/e/project/srna/output/result

for tsv in `ls *.rRNA.result.tsv`
do
cat $tsv |
    parallel --colsep '\t' -j 1 -k '
        echo "==> {1}"
        Rscript -e "
            x <- matrix(c({2}, {4}, {3}, {5}), nrow=2)
            x
            chisq.test(x)
        "
    ' > ../chi/$tsv.chi-square.txt
done
```

```bash
cd /mnt/e/project/srna/output/chi

parallel -j 3 " \
cat {}.result.tsv.chi-square.txt | perl ../../script/square.pl > {}.chi.tsv \
" ::: $(ls *.rRNA.*.txt | perl -p -e 's/\.result.+txt$//')
```

```bash
parallel -j 3 " \
cat {}.chi.tsv | tsv-join --filter-file name.tsv -k 1 --append-fields 2,3 > {}.result.tsv \
" ::: $(ls *.rRNA.chi.tsv | perl -p -e 's/\.chi\.tsv//')
```



