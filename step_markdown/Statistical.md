- [Bacterial related reads statistical information](#bacterial-related-reads-statistical-information)

    - [Basic percentages of reads aligned to plant or bacteria](#basic-percentages-of-reads-aligned-to-plant-or-bacteria)

    - [Filter by ratio of aligning to plant](#filter-by-ratio-of-aligning-to-plant)

# Bacterial related reads statistical information

In this markdown, I recorded all reads aligned to bacteria and their characristics.


## Basic percentages of reads aligned to plant

The sRNA-seq samples we gathered from NCBI were all completed in plant. So the original purpose was to detect the presence of sRNAs in plants. So the alignment first round in [sRNA_mapping_bowtie2.md](https://github.com/jdasfd/srna/blob/main/step_markdown/sRNA_mapping_bowtie2.md) provided two files: -N 0 reads aligning mode and -N 1 reads aligning mode.

Files originated from the -N 1 parameter contained all reads that could align to the plant genome without considering mismatches.


```bash
cd /mnt/e/project/srna/output/bam/plant

bash ../../script/all_file_count.sh > ../count/all_file.csv

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


## Ratio of reads aligned to bacteria / all non-plant reads among categories

* All 240 sRNA-seq files

```bash
cd /mnt/e/project/srna/output/bam/bacteria

parallel -j 4 " \
samtools idxstats -@ 3 {}.sort.bam | grep '*' | \
tsv-select -f 4 > ../../count/all/{}.unknown.tsv \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools idxstats -@ 3 {}.sort.bam | grep -v '*' | \
tsv-select -f 1,3 | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 > ../../count/all/{}.bac.tsv \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
cd /mnt/e/project/srna/output/count/all

parallel -j 6 " \
cat {}.bac.tsv | tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 > {}.all.tsv \
" ::: $(ls *.bac.tsv | perl -p -e 's/\.bac\.tsv$//')

parallel -j 12 " \
perl ../../../script/count.pl -a {}.all.tsv \
-r {}.unknown.tsv -o {}.tsv \
" ::: $(ls *.all.tsv | perl -p -e 's/\.all\.tsv$//')

rm *.all.tsv *.unknown.tsv
```

```bash
cd /mnt/e/project/srna/output/count/all

rm ../bacreads_group.tsv
# because next step will use >> so clear first

for file in `ls *.bac.tsv | perl -p -e 's/\.bac\.tsv$//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.tsv | awk -v name=$name -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
>> ../bacreads_group.tsv
done

cd ..
sed -i '1i\name\tgroup\tratio\tcatgry' bacreads_group.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
 -n 5 -f bacreads_group.tsv -t RNA_in_group -y "Bac-reads" -o ../figure/bacreads_group.pdf
# add the scale_y_continuous limits
```

#### After filtering

```bash
cat bacreads_group.tsv | \
tsv-join -H --filter-file plant_50.tsv --key-fields name \
> bacreads_group_50.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-n 5 -f bacreads_group_50.tsv -t RNA_in_group -y "Bac-reads" -o ../figure/bacreads_group_50.pdf
```

```bash
cat bacreads_group_50.tsv | \
tsv-summarize -H --group-by group,catgry --mean ratio | \
mlr --itsv --omd cat
```

| group | catgry | ratio_mean |
| --- | --- | --- |
| 3 | 1mis | 0.634710151378 |
| 1 | 1mis | 1.02873775248 |
| 4 | 1mis | 0.00620741048405 |
| 2 | 1mis | 1.50359081252 |
| 1 | aliall | 1.48569685888 |
| 4 | aliall | 0.208781526894 |
| 3 | aliall | 1.16217552165 |
| 2 | aliall | 2.14279758871 |
| 1 | unali | 3.40977514859 |
| 4 | unali | 0.0586593955001 |
| 2 | unali | 4.55701684177 |
| 3 | unali | 7.21581767507 |

```bash
cat bacreads_group_50.tsv | \
tsv-summarize -H --group-by group,catgry --median ratio | \
mlr --itsv --omd cat
```

| group | catgry | ratio_median |
| --- | --- | --- |
| 3 | 1mis | 0.212392565465 |
| 1 | 1mis | 0.425748152575 |
| 4 | 1mis | 0.0033455682172 |
| 2 | 1mis | 0.208726262452 |
| 1 | aliall | 0.937687844361 |
| 4 | aliall | 0.112857058589 |
| 3 | aliall | 0.750917485459 |
| 2 | aliall | 1.45728697456 |
| 1 | unali | 2.02597994901 |
| 4 | unali | 0.0434075553298 |
| 2 | unali | 3.1837129376 |
| 3 | unali | 1.20582662298 |

