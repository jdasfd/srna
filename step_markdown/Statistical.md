## Ratio of reads aligned to bacteria / all non-plant reads

Count all reads numbers from sort.bam files. The goal of this step is to acquire fraction of the reads aligned to bacteria from all reads. I wrote a shell script to reach the goal.

#### All 240 seq files

```bash
mkdir -p /mnt/e/project/srna/output/count
cd /mnt/e/project/srna/output/count
mkdir trna rrna mrna all

cd /mnt/e/project/srna/output/bam/bacteria
bash ../../../script/read_count.sh > ../../count/read_count.csv
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
labs(x = " ", y = "Bacterial reads / all reads")
ggsave(s, file = "../figure/read_count.pdf", width = 7, height = 4)
' read_count.csv
```

#### Remove those seq files after filter

```bash
cat read_count.csv | mlr --icsv --otsv cat | \
tsv-join -H --filter-file plant_50.tsv --key-fields 1 | \
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



## Plant sRNA reads distribution

*A. tha* annotation is relatively abundant with full information. Using `.gff` file, it is better using gene to calculate col 3 rather than using directly RNA annotation, such as tRNA *et. al.*. It almost the same using two different methods, though there will be a few lines of difference, *e.g.* miRNA will provide you 5p and 3p, but gene will just give you a region. Because of the existence of transcript splicing, using gene could directly give out the mRNA region to meet my expectations

#### All seq files

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
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/plant/Atha/Atha_rrna.bed \
{}.sort.bam > ../plantrna/{}.rrna.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/plant/Atha/Atha_mrna.bed \
{}.sort.bam > ../plantrna/{}.mrna.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/plant/Atha/Atha_mirna.bed \
{}.sort.bam > ../plantrna/{}.mirna.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/plant/Atha/Atha_snrna.bed \
{}.sort.bam > ../plantrna/{}.snrna.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/plant/Atha/Atha_snorna.bed \
{}.sort.bam > ../plantrna/{}.snorna.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/plant/Atha/Atha_lncrna.bed \
{}.sort.bam > ../plantrna/{}.lncrna.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/plant/Atha/Atha_ncrna.bed \
{}.sort.bam > ../plantrna/{}.ncrna.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
cd /mnt/e/project/srna/output/count
bash ../../script/plantali_count.sh > plantali_count.csv

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
ggsave(p, file = "/mnt/e/project/srna/output/figure/plantali.pdf", width = 9, height = 4)
' plantali_count.csv
```

#### After filtering

```bash
cat plantali_count.csv | mlr --icsv --otsv cat | \
tsv-join -H --filter-file plant_50.tsv --key-fields 1 | \
mlr --itsv --ocsv cat > plantali_count_50.csv
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
ggsave(p, file = "/mnt/e/project/srna/output/figure/plantali_50.pdf", width = 9, height = 4)
' plantali_count_50.csv
```



## Ratio of reads aligned to bacteria / all non-plant reads among categories

#### All seq files

```bash
cd /mnt/e/project/srna/output/bam/bacteria

parallel -j 4 " \
samtools idxstats -@ 3 {}.sort.bam | grep '*' | \
tsv-select -f 4 > ../../count/bacteria/{}.rest.tsv \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
cd /mnt/e/project/srna/output/count/all

parallel -j 6 " \
cat {}.all.tsv | tsv-join --filter-file ../../../name.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 > ../bacteria/{}.all.tsv \
" ::: $(ls *.all.tsv | perl -p -e 's/\.all\.tsv$//')
```

```bash
cd /mnt/e/project/srna/output/count/bacteria

parallel -j 12 " \
perl ../../../script/count.pl -a {}.all.tsv \
-r {}.rest.tsv -o {}.tsv \
" ::: $(ls *.all.tsv | perl -p -e 's/\.all\.tsv$//')

rm *.all.tsv *.rest.tsv

rm ../result.tsv
# because next step will use >> so clear first

for file in `ls | perl -p -e 's/\.tsv//'`
do
name=${file%%_*};
catgry=${file#*_};
cat ${file}.tsv | awk -v name=$name -v catgry=$catgry '{print name"\t"$1"\t"$2"\t"catgry}' \
>> ../result.tsv
done

cd ..
sed -i '1i\name\tgroup\tratio\tcatgry' result.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
 -n 5 -f result.tsv -t RNA_in_group -y "Bac-reads" -o ../figure/all_RNA_group.pdf
# add the scale_y_continuous limits
```

#### After filtering

```bash
cat result.tsv | tsv-join -H --filter-file plant_50.tsv --key-fields 1 > result_50.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-n 5 -f result_50.tsv -t RNA_in_group -y "Bac-reads" -o ../figure/all_RNA_group_50.pdf
```



## Ratio of sRNA reads / aligned reads from different RNA regions in groups

Using bed of rna to extract mapping reads from different RNA regions.

#### All seq files

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
-f result.trna.tsv -t tRNA_region -y "Bac-reads in tRNA" -o ../figure/trna_reads.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.rrna.tsv -t rRNA_region -y "Bac-reads in rRNA" -o ../figure/rrna_reads.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.mrna.tsv -t mRNA_region -y "Bac-reads in mRNA" -o ../figure/mrna_reads.pdf
```

#### After filtering

```bash
cat result.trna.tsv | tsv-join -H --filter-file plant_50.tsv --key-fields 1 > result_50.trna.tsv
cat result.rrna.tsv | tsv-join -H --filter-file plant_50.tsv --key-fields 1 > result_50.rrna.tsv
cat result.mrna.tsv | tsv-join -H --filter-file plant_50.tsv --key-fields 1 > result_50.mrna.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.trna.tsv -t tRNA_region -y "Bac-reads in tRNA" -o ../figure/trna_reads_50.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.rrna.tsv -t rRNA_region -y "Bac-reads in rRNA" -o ../figure/rrna_reads_50.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.mrna.tsv -t mRNA_region -y "Bac-reads in mRNA" -o ../figure/mrna_reads_50.pdf
```



## Bacteria occurred frequencies among different groups (waiting for update)

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
Rscript /mnt/e/project/srna/script/rna_plot.r -f {}.num.tsv \
-o ../figure/{}_freq.pdf -t {} -y frequencies \
" ::: $(ls *.tsv | perl -p -e 's/\.num\.tsv//')
```



## Sequence among all files

So the main reason I do this step is to select those frequently occurred among all sequence files, that is, the most likely sRNA appeared among *A. tha* sRNA-seq files.

#### All seq files

```bash
mkdir -p /mnt/e/project/srna/output/tier/trna
cd /mnt/e/project/srna/output/bam/rna

parallel -j 4 " \
samtools view -@ 2 {}.trna.sort.bam | \
tsv-select -f 3,10 > ../../tier/trna/{}.trna.tsv \
" ::: $(ls *.trna.sort.bam | perl -p -e 's/\.trna.+bam$//')
```

```bash
mkdir -p /mnt/e/project/srna/output/tier/among
mkdir -p /mnt/e/project/srna/output/tier/file
cd /mnt/e/project/srna/output/tier/trna

parallel -j 6 " \
cat {}_aliall.trna.tsv {}_1mis.trna.tsv {}_unali.trna.tsv | \
tsv-summarize --group-by 2 --count > ../file/{}.trna.tsv \
" ::: $(ls *.trna.tsv | perl -p -e 's/_.+tsv$//' | uniq)

cd ../file

rm ../among/all_seq.tsv # ensure >> will not just output after the old file
cat *.tsv | tsv-select -f 1 >> ../among/all_seq.tsv
cd ../among
cat all_seq.tsv | tsv-summarize --group-by 1 --count > all_seq.count.tsv

cat all_seq.count.tsv | tsv-summarize --group-by 2 --count | sed '1inum\tcount'> seq_num.tsv
```

```bash
Rscript -e '
library(ggplot2)
library(readr)
library(ggforce)
args <- commandArgs(T)
count <- read_tsv(args[1])
s <- ggplot (data = count, aes(x = num, y = count)) +
geom_bar(stat = "identity") +
theme(axis.ticks.x = element_blank()) +
labs(x = "file num", y = "frequency") +
scale_x_continuous(limits = c(0,220)) +
facet_zoom(ylim = c(0, 2000)) +
geom_col()
ggsave(s, file = "../../figure/seq_distri.pdf", width = 10, height = 4)
' seq_num.tsv
```

```bash
cat all_seq.count.tsv | tsv-filter --ge 2:120 > tier1.tsv
cat all_seq.count.tsv | tsv-filter --ge 2:60 --lt 2:120 > tier2.tsv
cat all_seq.count.tsv | tsv-filter --lt 2:60 > tier3.tsv
```

```bash
mkdir tier1 tier2 tier3
cd /mnt/e/project/srna/output/tier/trna

parallel -j 10 " \
perl /mnt/e/project/srna/script/select_seq.pl -i {}.trna.tsv \
-t ../among/tier1.tsv -o ../tier1/{}.tier1.tsv \
" ::: $(ls *.trna.tsv | perl -p -e 's/\.trna\.tsv$//')

parallel -j 10 " \
perl /mnt/e/project/srna/script/select_seq.pl -i {}.trna.tsv \
-t ../among/tier2.tsv -o ../tier2/{}.tier2.tsv \
" ::: $(ls *.trna.tsv | perl -p -e 's/\.trna\.tsv$//')

parallel -j 10 " \
perl /mnt/e/project/srna/script/select_seq.pl -i {}.trna.tsv \
-t ../among/tier3.tsv -o ../tier3/{}.tier3.tsv \
" ::: $(ls *.trna.tsv | perl -p -e 's/\.trna\.tsv$//')
```

```bash
cd /mnt/e/project/srna/output/tier/tier1

parallel -j 6 " \
cat {}.tier1.tsv | \
tsv-summarize --group-by 1 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-select -f 3,2 | perl ../../../script/filter.pl | \
tsv-join --filter-file ../../count/all/{}.all.tsv --key-fields 1 --append-fields 2 \
> {}.tsv \
" ::: $(ls *.tier1.tsv | perl -p -e 's/\.tier.+tsv$//')

rm *.tier1.tsv

bash ../../../script/group_rna_count.sh > ../name_count.tier1.tsv
```

```bash
cd /mnt/e/project/srna/output/tier/tier2

parallel -j 6 " \
cat {}.tier2.tsv | \
tsv-summarize --group-by 1 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-select -f 3,2 | perl ../../../script/filter.pl | \
tsv-join --filter-file ../../count/all/{}.all.tsv --key-fields 1 --append-fields 2 \
> {}.tsv \
" ::: $(ls *.tier2.tsv | perl -p -e 's/\.tier.+tsv$//')

rm *.tier2.tsv

bash ../../../script/group_rna_count.sh > ../name_count.tier2.tsv
```

```bash
cd /mnt/e/project/srna/output/tier/tier3

parallel -j 6 " \
cat {}.tier3.tsv | \
tsv-summarize --group-by 1 --count | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-select -f 3,2 | perl ../../../script/filter.pl | \
tsv-join --filter-file ../../count/all/{}.all.tsv --key-fields 1 --append-fields 2 \
> {}.tsv \
" ::: $(ls *.tier3.tsv | perl -p -e 's/\.tier.+tsv$//')

rm *.tier3.tsv

bash ../../../script/group_rna_count.sh > ../name_count.tier3.tsv
```

```bash
cd ..

cat name_count.tier1.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[3];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.tier1.tsv

cat name_count.tier2.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[3];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.tier2.tsv

cat name_count.tier3.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[3];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.tier3.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.tier1.tsv -t ">120_files" -y "Bac-reads in tRNA (T1)" -o ../figure/tier1_percent.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.tier2.tsv -t "60-120_files" -y "Bac-reads in tRNA (T2)" -o ../figure/tier2_percent.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.tier3.tsv -t "<60_files" -y "Bac-reads in tRNA (T3)" -o ../figure/tier3_percent.pdf
```

#### After filtering

```bash
cat result.tier1.tsv | tsv-join -H --filter-file ../count/plant_50.tsv --key-fields 1 > result_50.tier1.tsv
cat result.tier2.tsv | tsv-join -H --filter-file ../count/plant_50.tsv --key-fields 1 > result_50.tier2.tsv
cat result.tier3.tsv | tsv-join -H --filter-file ../count/plant_50.tsv --key-fields 1 > result_50.tier3.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.tier1.tsv -t ">120_files" -y "Bac-reads in tRNA (T1)" -o ../figure/tier1_percent_50.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.tier2.tsv -t "60-120_files" -y "Bac-reads in tRNA (T2)" -o ../figure/tier2_percent_50.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.tier3.tsv -t "<60_files" -y "Bac-reads in tRNA (T3)" -o ../figure/tier3_percent_50.pdf
```



## tRF region

tRF3/5 regions and other_tRNA regions were extracted from the tRNA.bed file according to tRF characteristics.

#### All seq files

```bash
cd /mnt/e/project/srna/annotation/bacteria
```

```bash
cat trna.bed | perl -e 'while(<>){
    chomp;
    @a = split/\t/,$_;
    $end1 = $a[1] + 25;
    $end2 = $a[2] - 25;
    print"$a[0]\t$a[1]\t$end1\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\n";
    print"$a[0]\t$end2\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\n";
}' > trf3_5.bed
```

```bash
cat trna.bed | perl -e 'while(<>){
    chomp;
    @a = split/\t/,$_;
    $start = $a[1] + 26;
    $end = $a[2] - 26;
    print"$a[0]\t$start\t$end\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\n";
}' > other_trf.bed
```

```bash
mkdir -p /mnt/e/project/srna/output/bam/trf
cd /mnt/e/project/srna/output/bam/rna

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/bacteria/trf3_5.bed \
{}.trna.sort.bam > ../trf/{}.trf3_5.bam \
" ::: $(ls *.trna.sort.bam | perl -p -e 's/\.trna.+bam$//')

parallel -j 4 " \
samtools view -@ 3 -bh -L ../../../annotation/bacteria/other_trf.bed \
{}.trna.sort.bam > ../trf/{}.other_trf.bam \
" ::: $(ls *.trna.sort.bam | perl -p -e 's/\.trna.+bam$//')
```

``` bash
cd /mnt/e/project/srna/output/bam/trf

parallel -j 3 " 
samtools sort -@ 4 {1}.bam > {1}.sort.bam 
samtools index {1}.sort.bam
" ::: $(ls *.bam | perl -p -e 's/\.bam$//')

parallel -j 3 " \
rm {}.bam \
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
cd /mnt/e/project/srna/output/count
mkdir trf3_5 other_trf
cd /mnt/e/project/srna/output/bam/trf

parallel -j 6 " \
samtools idxstats {}.trf3_5.sort.bam | \
tsv-select -f 1,3 | grep -v '*' | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 \
> ../../count/trf3_5/{}.trf3_5.tsv \
" ::: $(ls *.trf3_5.sort.bam | perl -p -e 's/\.trf.+bam$//')

parallel -j 6 " \
samtools idxstats {}.other_trf.sort.bam | \
tsv-select -f 1,3 | grep -v '*' | \
tsv-join --filter-file ../../../rawname.tsv --key-fields 1 --append-fields 2 | \
tsv-summarize --group-by 3 --sum 2 \
> ../../count/other_trf/{}.other_trf.tsv \
" ::: $(ls *.other_trf.sort.bam | perl -p -e 's/\.other.+bam$//')
```

```bash
cd /mnt/e/project/srna/output/count/trf3_5

for file in `ls *.trf3_5.tsv | perl -p -e 's/\.trf.+tsv//'`
do
cat ${file}.trf3_5.tsv | \
tsv-join --filter-file ../all/${file}.all.tsv --key-fields 1 --append-fields 2 \
> ${file}.tsv
done

rm *.trf3_5.tsv

bash ../../../script/group_rna_count.sh > ../name_count.trf3_5.tsv
```

```bash
cd /mnt/e/project/srna/output/count/other_trf

for file in `ls *.other_trf.tsv | perl -p -e 's/\.other.+tsv//'`
do
cat ${file}.other_trf.tsv | \
tsv-join --filter-file ../all/${file}.all.tsv --key-fields 1 --append-fields 2 \
> ${file}.tsv
done

rm *.other_trf.tsv

bash ../../../script/group_rna_count.sh > ../name_count.other_trf.tsv
```

```bash
cd /mnt/e/project/srna/output/count

cat name_count.trf3_5.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[3];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.trf3_5.tsv

cat name_count.other_trf.tsv | perl -n -e 'while(<>){chomp;
@a=split/\t/,$_;
$b=$a[2]*100/$a[3];
print"$a[0]\t$a[1]\t$b\t$a[4]\n";
}' | sed -e '1i\name\tgroup\tratio\tcatgry' > result.other_trf.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.trf3_5.tsv -t tRF-3/5_region -y "tRF3/5 Bac-reads" -o ../figure/trf3_5_percent.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result.other_trf.tsv -t other-tRNA_region -y "other tRNA Bac-reads" -o ../figure/other_trf_percent.pdf
```

#### After filtering

```bash
cat result.trf3_5.tsv | tsv-join -H --filter-file plant_50.tsv --key-fields 1 > result_50.trf3_5.tsv
cat result.other_trf.tsv | tsv-join -H --filter-file plant_50.tsv --key-fields 1 > result_50.other_trf.tsv
```

```bash
Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.trf3_5.tsv -t tRF-3/5_region -y "tRF3/5 Bac-reads" -o ../figure/trf3_5_percent_50.pdf

Rscript /mnt/e/project/srna/script/rna_percent.r \
-f result_50.other_trf.tsv -t other-tRNA_region -n 25 -y "other tRNA Bac-reads" -o ../figure/other_trf_percent_50.pdf
```






## Extract tRF sequence from different bacteria

### Get the *A. thaliana* transcripts

```bash
mkdir -p /mnt/e/project/srna/genome/plant_CDS/Atha
cd /mnt/e/project/srna/genome/plant_CDS/Atha
gffread ../../../annotation/plant/Atha/Atha.gff -g ../../plant/Atha/Atha.fna \
-w Atha_transcript.fa -x Atha_CDS.fa
```

### Extract sequence from the tRF3/5 region and high frequency sequence

```bash
mkdir -p /mnt/e/project/srna/output/seq/trf
cd /mnt/e/project/srna/output/seq/trf
mkdir 1 2 3
cd /mnt/e/project/srna/output/bam/trf

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:1 \
> ../../seq/trf/1/{}.tsv \
" ::: $(ls *.trf3_5.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:2 \
> ../../seq/trf/2/{}.tsv \
" ::: $(ls *.trf3_5.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:3 \
> ../../seq/trf/3/{}.tsv \
" ::: $(ls *.trf3_5.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
mkdir -p /mnt/e/project/srna/output/seq/other
cd /mnt/e/project/srna/output/seq/other
mkdir 1 2 3
cd /mnt/e/project/srna/output/bam/trf

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:1 \
> ../../seq/other/1/{}.tsv \
" ::: $(ls *.other_trf.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:2 \
> ../../seq/other/2/{}.tsv \
" ::: $(ls *.other_trf.sort.bam | perl -p -e 's/\.sort\.bam$//')

parallel -j 3 " \
samtools view -@ 4 {}.sort.bam | \
tsv-select -f 3,1,4,10 | \
tsv-join --key-fields 1 --filter-file ../../../rawname.tsv --append-fields 2 | \
tsv-select -f 5,1,2,3,4 | \
tsv-join --key-fields 1 --filter-file ../../../name.tsv --append-fields 2 | \
tsv-filter --eq 6:3 \
> ../../seq/other/3/{}.tsv \
" ::: $(ls *.other_trf.sort.bam | perl -p -e 's/\.sort\.bam$//')
```

```bash
cd /mnt/e/project/srna/output/seq/trf/1

cat *_1mis.trf3_5.tsv | tsv-select -f 2,4,5 | \
perl -n -e 'chomp;@a=split/\t/,$_;$l=length($a[2]);print"$a[0]\t$a[1]\t$l\t$a[2]\n";' | \
tsv-summarize --group-by 1,2,3,4 --count | \
perl -n -e 'chomp;@a = split/\t/,$_;print">$a[0]_$a[1]_$a[2]\n";print"$a[3]\n"' \
> ../1mis.1.fasta
```



## Gene list  extraction

```bash
cat Atha_gene.gff | perl -n -e 'chomp;
@a = split/\t/,$_;
@b = split/;/,$a[8];
$b[0] =~ s/^ID=gene-//;
print "$a[0]\t$b[0]\t$a[3]\t$a[4]\t$a[6]\n";
'
```

```bash
cd /mnt/e/project/srna/output/bam/trf
```



## ClusterProfiler for GO and KEGG analysis

```bash
cd /mnt/e/project/srna/output/kegg


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
# 1       plants_mart      Ensembl Plants Genes 52
# 2 plants_variations Ensembl Plants Variations 52

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



---

# Waiting for update



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



----------------------------------------------------



```R
bar <- ggplot(data = name, mapping = aes(x = 'category', y = count, fill = cate)) +
geom_bar(stat = 'identity', position = 'stack') +
coord_polar(theta = 'y') +
labs(x = '', y = '', title = '') +

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



```bash
Rscript -e '
library(readr)
library(ggplot2)
library(dplyr)
library(forcats)

count <- read_tsv("all_seq.count.tsv")
count$num <- as.character(count$num)

tplot <- count %>%
mutate(seq = fct_reorder(seq, desc(count))) %>%
ggplot(aes(x = seq, y = num)) +
geom_bar(stat="identity") +
theme(axis.text.x = NULL)
ggsave(tplot, file = "freq.pdf", width = 7, height = 4)
'
```



## Reads TPM among different groups (waiting for update)

The goal is to filter out those reads among one file.

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

tsv-select -H -f file,group,catgry,ttpm tpm.tsv | sed '1d' | sed '1ifile\tgroup\tcatgry\tnum' > ttpm.tsv
tsv-select -H -f file,group,catgry,rtpm tpm.tsv | sed '1d' | sed '1ifile\tgroup\tcatgry\tnum' > rtpm.tsv
tsv-select -H -f file,group,catgry,mtpm tpm.tsv | sed '1d' | sed '1ifile\tgroup\tcatgry\tnum' > mtpm.tsv

Rscript /mnt/e/project/srna/script/rna_plot.r -f ttpm.tsv -o ../figure/ttpm.pdf -t trna -y TPM
Rscript /mnt/e/project/srna/script/rna_plot.r -f rtpm.tsv -o ../figure/rtpm.pdf -t rrna -y TPM
Rscript /mnt/e/project/srna/script/rna_plot.r -f mtpm.tsv -o ../figure/mtpm.pdf -t mrna -y TPM
```



## Using RIsearch2 to  predict RNA-RNA interaction

### Create suffix array

```bash
mkdir -p /mnt/e/project/srna/output/risearch
cd /mnt/e/project/srna/genome/plant_CDS/Atha
risearch2.x -c Atha_transcript.fa -t 12 -o Atha_transcript.suf
```

### Interaction prediction

```bash
mkdir -p /mnt/e/project/srna/output/risearch
cd /mnt/e/project/srna/output/risearch
risearch2.x -q ../seq/trf/1mis.1.fasta -i ../../genome/plant_CDS/Atha/Atha_transcript.suf -s 7 -e -10
```
