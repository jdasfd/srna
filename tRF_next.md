# Specified tRF manipulating process

After I finished my sRNA analysis in sRNA.md. I realized that there are some specific phenomena I could use more detailed process dealing with them. So this process will be a continued manipulating for sRNA in plant sRNA-seq.

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

```bash
parallel -j 3 " \
bowtie2 -q {}_trimmed.fq.gz -N 0 \
-x ../genome/plant/Atha/Atha \
--al-gz ../output/fastq/{}_plantaliall.fq.gz \
--un-gz ../output/fastq/{}_plantunali.fq.gz \
--threads 4 -S ../output/bam/plant/{}_plantall.sam \
" ::: $(ls SRR*.fq.gz | perl -p -e 's/_trimmed.+gz$//')
```

```bash
bsub -q mpi -n 24 -J ali2plant -o . "bash alignall.sh"
```

Convert sam to bam for reducing storage stress

```bash
cd /mnt/e/project/srna/output/bam/plant

parallel -j 3 " 
samtools sort -@ 4 {1}.sam > {1}.sort.bam 
samtools index {1}.sort.bam 
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')
```

### Aligning those reads from non-plant sources to bacteria

We chose 365 bacterial genomes from 190 species as our target bacteria.

####  Indexing

```bash
cd /mnt/e/project/srna/genome/bacteria

bowtie2-build --threads 12 --quiet bacteria.fna bacteria
```

### Aligning

```bash
mkdir -p /mnt/e/project/srna/output/bam/bacteria
```

There will be 2 different files: plant sRNA reads and non-plant sRNA reads. We aligned them to bacteria separately. Their biological information were plant and bacteria homologous reads and non-plant reads directly from bacteria.

```bash
cd /mnt/e/project/srna/output/fastq

parallel -j 3 " \
bowtie2 -q {}_plantunali.fq.gz \
-x ../../genome/bacteria/bacteria --threads 4 -S ../bam/bacteria/{}_unali.sam \
" ::: $(ls SRR*_plantunali.fq.gz | perl -p -e 's/_plant.+gz$//')
```

```bash
bsub -q mpi -n 24 -J unali2bac -o . "bash unali2bac.sh"
```

```bash
cd /mnt/e/project/srna/output/fastq

parallel -j 3 " \
bowtie2 -q {}_plantaliall.fq.gz \
-x ../../genome/bacteria/bacteria --threads 4 -S ../bam/bacteria/{}_ali.sam \
" ::: $(ls SRR*_plantaliall.fq.gz | perl -p -e 's/_plant.+gz$//')
```

```bash
bsub -q mpi -n 24 -J ali2bac -o . "bash ali2bac.sh"
```



## Reads from bacteria percentages

```bash
echo "name,ratio,group";
cd /mnt/e/project/trf/output/bam/plant
for file in `ls SRR*.sort.bam | perl -p -e 's/_.+bam$//'`
do
cd ..;
unalip=`samtools view -c -f 4 -@ 10 ./plant/${file}_plantall.sort.bam`;
allp=`samtools view -c -@ 10 ./plant/${file}_plantall.sort.bam`;
alib=`samtools view -c -F 4 -@ 10 ./bacteria/${file}_unali.sort.bam`;
un=`echo "scale=2;$unalip*100/$allp" | bc | awk '{printf "%.2f", $0}'`;
alitob=`echo "scale=2;$alib*100/$unalip" | bc | awk '{printf "%.2f", $0}'`;
echo "${file},${un},unali";
echo "${file},${alitob},alitobac";
cd ./plant;
done
```

```bash
bash ../../script/percentage.sh | tee /mnt/e/project/srna/output/bam/percentage.csv
cd /mnt/e/project/srna/output/bam
cp percentage.csv /mnt/c/Users/59717/Documents/
```



```R
library(ggplot2)
library(readr)

count <- read.csv("read_count.csv")
s <- ggplot (data = count, aes(x = group, y = num)) +
geom_boxplot() + 
geom_jitter(aes(color = name)) +
theme(legend.position = 'none')

p <- ggplot(data = per, aes(x = name, y = ratio, group = group, fill = group)) +
geom_bar(stat = "identity", position = "dodge") +
labs(x = NULL) +
theme(axis.text.x = element_blank())
```


